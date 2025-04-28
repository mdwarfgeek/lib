#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lfa.h"

/* Implementation of GP model following Foreman-Mackey et al. (2017)
   "Fast and Scalable Gaussian Process Modeling with Applications to
   Astronomical Time Series", 2017, AJ, 154, 220.  Compared to their
   implementation this is less optimized / efficient but has minimal
   dependencies. */

/* Convenience function to generate SHO kernel, result array should
   have 8 elements allocated in general to accommodate q < 0.5 unless
   this is known not to be the case.  The number of kernels actually
   used is returned. */

int fsgp_kern_sho (double s0, double w0, double q, double *kern) {
  double tmp;
  int nkern;

  if(q < 0.5) {
    tmp = sqrt(1.0 - 4.0*q*q);

    kern[0] = 0.5 * s0 * w0 * q * (1.0 + 1.0/tmp);
    kern[1] = 0;
    kern[2] = 0.5 * w0 * (1.0 - tmp) / q;
    kern[3] = 0;

    kern[4] = 0.5 * s0 * w0 * q * (1.0 - 1.0/tmp);
    kern[5] = 0;
    kern[6] = 0.5 * w0 * (1.0 + tmp) / q;
    kern[7] = 0;

    nkern = 2;
  }
  else {
    if(q == 0.5)
      tmp = 1.0e-3;
    else
      tmp = sqrt(4.0*q*q - 1.0);
    
    kern[0] = s0 * w0 * q;
    kern[1] = kern[0] / tmp;
    kern[2] = 0.5 * w0 / q;
    kern[3] = kern[2] * tmp;

    nkern = 1;
  }
  
  return(nkern);
}

/* Convenience function to generate approximate Matern 3/2 kernel
   k(tau) = var * (1 + sqrt(3) * tau / rho) exp(-sqrt(3) * tau / rho) */

int fsgp_kern_matern (double var, double rho, double eps, double *kern) {
  double w0;
  int nkern;

  w0 = sqrt(3) / rho;
  
  /* var * exp(-w0*tau) (cos(eps*tau) + w0*sin(eps*tau)/eps) */
  kern[0] = var;
  kern[1] = var * w0 / eps;
  kern[2] = w0;
  kern[3] = eps;

  nkern = 1;

  return(nkern);
}

/* Check valididy of kernel where coefficients can take arbitrary values.
   Simple version where we require all nkern sub-kernels to be individually
   positive definite.  Returns count of valid sub-kernels.  Shouldn't be
   used for SHO or Matern kernels from above, which are always valid but
   correspond to the condition of equality in the validity test criterion
   below so are prone to numerical issues. */

int fsgp_kern_valid (double *kern, int nkern) {
  int jkern, nvalid;
  double aj, bj, cj, dj;
  
  nvalid = 0;
  for(jkern = 0; jkern < nkern; jkern++) {
    aj = kern[4*jkern];
    bj = kern[4*jkern+1];
    cj = kern[4*jkern+2];
    dj = kern[4*jkern+3];

    if(fabs(bj * dj) <= aj * cj)
      nvalid++;
  }

  return(nvalid);
}

/* Compute Cholesky factorization of kernel K = L D L^T where
   L = I + tril(U W^T)
   with the reparameterization suggested in Sect. 5.2 applied.

   The kernel coefficients are given as an array of nkern*4
   elements where the most rapidly varying direction are the
   four coefficients a_j, b_j, c_j, d_j for the ith kernel
   element in the summation.  The kernel is then:

   K(n,m)   = sigma_n^2 delta(n,m) + sum(j=1,nkern) K_j(n,m)
   K_j(n,m) = exp(-c_j dt) (a_j cos(d_j dt) + b_j sin(d_j dt))

   where dt is the time difference between points n,m.

   Data points must be in time order. */

int fsgp_compute (struct fsgp_fac *fac,
                  double *kern, int nkern,
                  double *t, double *yerr, int ndp) {
  double *kerncopy = (double *) NULL;
  double *tcopy = (double *) NULL;
  double *yvar = (double *) NULL;
  double *phi = (double *) NULL;
  double *u = (double *) NULL;
  double *v = (double *) NULL;
  double *w = (double *) NULL;
  double *d = (double *) NULL;
  double *sqrtd = (double *) NULL;
  double *s = (double *) NULL;
  int nckern;
  int idp, offp, off, jkern, jckern, kckern;
  double aj, bj, cj, dj, sinarg, cosarg, exparg, sumaj;
  double sum;
  
  /* Return error if there are no data points */
  if(ndp < 1)
    return(-1);

  nckern = 2*nkern;

  /* Allocate arrays */
  kerncopy = (double *) malloc(4*nkern * sizeof(double));
  tcopy = (double *) malloc(ndp * sizeof(double));
  yvar = (double *) malloc(ndp * sizeof(double));
  phi = (double *) malloc(ndp*nckern * sizeof(double));
  u = (double *) malloc(ndp*nckern * sizeof(double));
  v = (double *) malloc(ndp*nckern * sizeof(double));
  w = (double *) malloc(ndp*nckern * sizeof(double));
  d = (double *) malloc(ndp * sizeof(double));
  sqrtd = (double *) malloc(ndp * sizeof(double));
  s = (double *) malloc(nckern*nckern * sizeof(double));
  if(!kerncopy || !tcopy || !yvar || !phi || !u || !v || !w ||
     !sqrtd || !d || !s)
    goto error;

  /* Copies of inputs needed for "predict" operation */
  memcpy(kerncopy, kern, 4*nkern * sizeof(double));
  memcpy(tcopy, t, ndp * sizeof(double));
  
  /* Sum of aj */
  sumaj = 0;

  for(jkern = 0; jkern < nkern; jkern++) {
    aj = kern[4*jkern];
    sumaj += aj;
  }
  
  /* First term */
  yvar[0] = yerr[0]*yerr[0];
  d[0] = yvar[0] + sumaj;
  sqrtd[0] = sqrt(d[0]);

  for(jkern = 0; jkern < nkern; jkern++) {
    aj = kern[4*jkern];
    bj = kern[4*jkern+1];
    cj = kern[4*jkern+2];
    dj = kern[4*jkern+3];
    
    inline_sincos(dj * t[0], sinarg, cosarg);
    
    u[2*jkern]     = aj * cosarg + bj * sinarg;
    u[2*jkern+1]   = aj * sinarg - bj * cosarg;
    v[2*jkern]     = cosarg;
    v[2*jkern+1]   = sinarg;
    w[2*jkern]     = cosarg / d[0];  /* directly calculate W */
    w[2*jkern+1]   = sinarg / d[0];
    phi[2*jkern]   = 0;
    phi[2*jkern+1] = 0;
  }

  for(jckern = 0; jckern < nckern; jckern++)
    for(kckern = 0; kckern <= jckern; kckern++)
      s[jckern*nckern+kckern] = 0;

  /* Subsequent terms */
  for(idp = 1; idp < ndp; idp++) {
    offp = (idp-1) * nckern;
    off = idp * nckern;
    
    for(jkern = 0; jkern < nkern; jkern++) {
      aj = kern[4*jkern];
      bj = kern[4*jkern+1];
      cj = kern[4*jkern+2];
      dj = kern[4*jkern+3];
      
      inline_sincos(dj * t[idp], sinarg, cosarg);

      u[off + 2*jkern]     = aj * cosarg + bj * sinarg;
      u[off + 2*jkern+1]   = aj * sinarg - bj * cosarg;
      v[off + 2*jkern]     = cosarg;
      v[off + 2*jkern+1]   = sinarg;

      exparg = exp(-cj * (t[idp]-t[idp-1]));
      
      phi[off + 2*jkern]   = exparg;
      phi[off + 2*jkern+1] = exparg;
    }

    /* Update S and accumulate sums for D */
    sum = 0;

    for(jckern = 0; jckern < nckern; jckern++) {
      /* Off-diagonal elements */
      for(kckern = 0; kckern < jckern; kckern++) {
        s[jckern*nckern+kckern] = phi[off + jckern] * phi[off + kckern] * (s[jckern*nckern+kckern] + d[idp-1] * w[offp + jckern] * w[offp + kckern]);
        sum += 2 * s[jckern*nckern+kckern] * u[off+jckern] * u[off+kckern];
      }

      /* Diagonal elements */
      s[jckern*nckern+jckern] = phi[off + jckern] * phi[off + jckern] * (s[jckern*nckern+jckern] + d[idp-1] * w[offp + jckern] * w[offp + jckern]);
      sum += s[jckern*nckern+jckern] * u[off+jckern] * u[off+jckern];
    }

    /* Calculate D */
    yvar[idp] = yerr[idp]*yerr[idp];
    d[idp] = yvar[idp] + sumaj - sum;
    sqrtd[idp] = sqrt(d[idp]);

    /* Calculate W */
    for(jckern = 0; jckern < nckern; jckern++) {
      sum = 0;
      
      for(kckern = 0; kckern <= jckern; kckern++)
        sum += s[jckern*nckern+kckern] * u[off+kckern];

      for(kckern = jckern+1; kckern < nckern; kckern++)
        sum += s[kckern*nckern+jckern] * u[off+kckern];

      w[off + jckern] = (v[off + jckern] - sum) / d[idp];
    }
  }

  free((void *) s);
  s = (double *) NULL;

  fac->kern = kerncopy;
  fac->sumaj = sumaj;
  fac->t = tcopy;
  fac->yvar = yvar;
  fac->phi = phi;
  fac->u = u;
  fac->v = v;
  fac->w = w;
  fac->d = d;
  fac->sqrtd = sqrtd;
  fac->nkern = nkern;
  fac->nckern = nckern;
  fac->ndp = ndp;
  
  return(0);
  
 error:
  if(kerncopy)
    free((void *) kerncopy);
  if(tcopy)
    free((void *) tcopy);
  if(yvar)
    free((void *) yvar);
  if(phi)
    free((void *) phi);
  if(u)
    free((void *) u);
  if(v)
    free((void *) v);
  if(w)
    free((void *) w);
  if(d)
    free((void *) d);
  if(sqrtd)
    free((void *) sqrtd);
  if(s)
    free((void *) s);

  return(-1);
}

/* Apply inverse of kernel */

void fsgp_apply (struct fsgp_fac *fac, double *fg,
                 double *y, double *z) {
  int idp, offp, off, jckern;
  double sum;

  /* Forward substitution */
  for(jckern = 0; jckern < fac->nckern; jckern++)
    fg[jckern] = 0;
    
  z[0] = y[0];
    
  for(idp = 1; idp < fac->ndp; idp++) {
    offp = (idp-1) * fac->nckern;
    off = idp * fac->nckern;
      
    sum = 0;
      
    for(jckern = 0; jckern < fac->nckern; jckern++) {
      fg[jckern] = fac->phi[off+jckern] * (fg[jckern] + fac->w[offp+jckern] * z[idp-1]);
      sum += fac->u[off+jckern] * fg[jckern];
    }
      
    z[idp] = y[idp] - sum;
  }
    
  /* Back substitution */
  for(jckern = 0; jckern < fac->nckern; jckern++)
    fg[jckern] = 0;
    
  z[fac->ndp-1] /= fac->d[fac->ndp-1];
    
  for(idp = fac->ndp-1; idp > 0; idp--) {
    offp = (idp-1) * fac->nckern;
    off = idp * fac->nckern;
      
    sum = 0;
      
    for(jckern = 0; jckern < fac->nckern; jckern++) {
      fg[jckern] = fac->phi[off+jckern] * (fg[jckern] + fac->u[off+jckern] * z[idp]);
      sum += fac->w[offp+jckern] * fg[jckern];
    }
      
    z[idp-1] = z[idp-1] / fac->d[idp-1] - sum;
  }
}

/* Apply inverse of kernel to diagonal matrix */

static double fsgp_apply_diag (struct fsgp_fac *fac, double *fg, double *z,
                               double yd, int irhs) {
  int idp, offp, off, jckern;
  double sum;

  /* Generate this right-hand side */
  for(idp = 0; idp < fac->ndp; idp++)
    z[idp] = 0;
  
  z[irhs] = yd;
  
  /* Forward substitution */
  for(jckern = 0; jckern < fac->nckern; jckern++)
    fg[jckern] = 0;
  
  for(idp = 1; idp < fac->ndp; idp++) {
    offp = (idp-1) * fac->nckern;
    off = idp * fac->nckern;
    
    sum = 0;
    
    for(jckern = 0; jckern < fac->nckern; jckern++) {
      fg[jckern] = fac->phi[off+jckern] * (fg[jckern] + fac->w[offp+jckern] * z[idp-1]);
      sum += fac->u[off+jckern] * fg[jckern];
    }
    
    z[idp] -= sum;
  }
  
  /* Back substitution, only need to go until we've computed the
     element of interest, so loop terminates early. */
  for(jckern = 0; jckern < fac->nckern; jckern++)
    fg[jckern] = 0;
  
  z[fac->ndp-1] /= fac->d[fac->ndp-1];
  
  for(idp = fac->ndp-1; idp > irhs; idp--) {
    offp = (idp-1) * fac->nckern;
    off = idp * fac->nckern;
    
    sum = 0;
    
    for(jckern = 0; jckern < fac->nckern; jckern++) {
      fg[jckern] = fac->phi[off+jckern] * (fg[jckern] + fac->u[off+jckern] * z[idp]);
      sum += fac->w[offp+jckern] * fg[jckern];
    }
    
    z[idp-1] = z[idp-1] / fac->d[idp-1] - sum;
  }

  return(z[irhs]);
}

/* "Predict" (interpolation/extrapolation) function based on result
   from "apply" above.  Data points must be in time order. */

int fsgp_predict (struct fsgp_fac *fac, double *y,
                  double *tpred, double *ypred, double *varpred, int npred) {
  double *z = (double *) NULL;
  double *fg = (double *) NULL;
  
  double *sinarg = (double *) NULL;
  double *cosarg = (double *) NULL;
  double ttmp;
  int ipred, offpred;
  
  double *q = (double *) NULL;
  int idp, off, offn, jkern, jckern;
  double aj, bj, cj, dj, exparg, u, v, mu;

  double *ktt = (double *) NULL;
  double *ktp = (double *) NULL;
  double *ztmp = (double *) NULL;
  double sum, dt, sgn, kv;

  /* Get out of here if there's no work to do */
  if(npred < 1)
    return(0);

  /* Calculate z */
  z = (double *) malloc(fac->ndp * sizeof(double));
  fg = (double *) malloc(fac->nckern * sizeof(double));
  if(!z || !fg)
    goto error;
  
  fsgp_apply(fac, fg, y, z);
  
  if(tpred) {    
    /* Precalculate sin, cos terms and initialize outputs */
    sinarg = (double *) malloc(npred*fac->nkern * sizeof(double));
    cosarg = (double *) malloc(npred*fac->nkern * sizeof(double));
    if(!sinarg || !cosarg)
      goto error;

    for(ipred = 0; ipred < npred; ipred++) {
      off = ipred*fac->nkern;
      ttmp = tpred[ipred];

      for(jkern = 0; jkern < fac->nkern; jkern++) {
        dj = fac->kern[4*jkern+3];
        inline_sincos(dj * ttmp, sinarg[off+jkern], cosarg[off+jkern]);
      }

      ypred[ipred] = 0;
    }
  
    /* Allocate workspace */
    q = (double *) malloc(fac->nckern * sizeof(double));
    if(!q)
      goto error;

    /* Forward pass: init Q */
    for(jckern = 0; jckern < fac->nckern; jckern++)
      q[jckern] = 0;

    /* Extrapolation prior to start, all means are zero so just skip them */
    ipred = 0;
    while(ipred < npred && tpred[ipred] < fac->t[0])
      ipred++;

    /* Interpolation */
    for(idp = 0; idp < fac->ndp-1; idp++) {
      off = idp * fac->nckern;
      offn = (idp+1) * fac->nckern;

      for(jckern = 0; jckern < fac->nckern; jckern++)
        q[jckern] = (q[jckern] + z[idp] * fac->v[off+jckern]) * fac->phi[offn+jckern];

      while(ipred < npred && tpred[ipred] < fac->t[idp+1]) {
        offpred = ipred * fac->nkern;

        mu = 0;
        for(jkern = 0; jkern < fac->nkern; jkern++) {
          aj = fac->kern[4*jkern];
          bj = fac->kern[4*jkern+1];
          cj = fac->kern[4*jkern+2];
        
          exparg = exp(-cj * (tpred[ipred]-fac->t[idp+1]));
        
          u = aj * cosarg[offpred+jkern] + bj * sinarg[offpred+jkern];
          mu += q[2*jkern] * u * exparg;
        
          u = aj * sinarg[offpred+jkern] - bj * cosarg[offpred+jkern];
          mu += q[2*jkern+1] * u * exparg;
        }
      
        ypred[ipred] += mu;
      
        ipred++;
      }
    }

    /* Extrapolation after the end */
    off = (fac->ndp-1) * fac->nckern;

    for(jckern = 0; jckern < fac->nckern; jckern++)
      q[jckern] = q[jckern] + z[fac->ndp-1] * fac->v[off+jckern];

    while(ipred < npred) {
      offpred = ipred * fac->nkern;
    
      mu = 0;
      for(jkern = 0; jkern < fac->nkern; jkern++) {
        aj = fac->kern[4*jkern];
        bj = fac->kern[4*jkern+1];
        cj = fac->kern[4*jkern+2];
      
        exparg = exp(-cj * (tpred[ipred]-fac->t[fac->ndp-1]));
      
        u = aj * cosarg[offpred+jkern] + bj * sinarg[offpred+jkern];
        mu += q[2*jkern] * u * exparg;
      
        u = aj * sinarg[offpred+jkern] - bj * cosarg[offpred+jkern];
        mu += q[2*jkern+1] * u * exparg;
      }
    
      ypred[ipred] += mu;
    
      ipred++;
    }

    /* Reverse pass: init Q */
    for(jckern = 0; jckern < fac->nckern; jckern++)
      q[jckern] = 0;
  
    /* Extrapolation off end, all means are zero so just skip them */
    ipred = npred - 1;
    while(ipred >= 0 && tpred[ipred] >= fac->t[fac->ndp-1])
      ipred--;

    /* Interpolation */
    for(idp = fac->ndp-1; idp > 0; idp--) {
      off = idp * fac->nckern;

      for(jckern = 0; jckern < fac->nckern; jckern++)
        q[jckern] = (q[jckern] + z[idp] * fac->u[off+jckern]) * fac->phi[off+jckern];

      while(ipred >= 0 && tpred[ipred] >= fac->t[idp-1]) {
        offpred = ipred * fac->nkern;

        mu = 0;
        for(jkern = 0; jkern < fac->nkern; jkern++) {
          cj = fac->kern[4*jkern+2];
        
          exparg = exp(-cj * (fac->t[idp-1]-tpred[ipred]));
        
          v = cosarg[offpred+jkern];
          mu += q[2*jkern] * v * exparg;
        
          v = sinarg[offpred+jkern];
          mu += q[2*jkern+1] * v * exparg;
        }
      
        ypred[ipred] += mu;
      
        ipred--;
      }
    }

    /* Extrapolation off the bottom */
    for(jckern = 0; jckern < fac->nckern; jckern++)
      q[jckern] = q[jckern] + z[0] * fac->u[jckern];

    while(ipred >= 0) {
      offpred = ipred * fac->nkern;
    
      mu = 0;
      for(jkern = 0; jkern < fac->nkern; jkern++) {
        cj = fac->kern[4*jkern+2];
      
        exparg = exp(-cj * (fac->t[0]-tpred[ipred]));
      
        v = cosarg[offpred+jkern];
        mu += q[2*jkern] * v * exparg;
      
        v = sinarg[offpred+jkern];
        mu += q[2*jkern+1] * v * exparg;
      }
    
      ypred[ipred] += mu;
    
      ipred--;
    }

    /* Covariance matrix, if requested (expensive) */
    if(varpred) {
      /* Allocate workspace for K(t,t*) */
      ktt = (double *) malloc(fac->ndp * sizeof(double));
      ktp = (double *) malloc(fac->ndp * sizeof(double));
      if(!ktt || !ktp)
        goto error;

      /* Form K(t,t*) directly using precomputed sin,cos tables.
         Exponentials are computed for each element to avoid numerical
         overflow rather than using precomputed ratios from "phi". */
      for(ipred = 0; ipred < npred; ipred++) {
        offpred = ipred * fac->nkern;
      
        for(idp = 0; idp < fac->ndp; idp++) {
          off = idp * fac->nckern;
        
          sum = 0;

          for(jkern = 0; jkern < fac->nkern; jkern++) {
            aj = fac->kern[4*jkern];
            bj = fac->kern[4*jkern+1];
            cj = fac->kern[4*jkern+2];

            /* tpred - t and sign */
            dt = tpred[ipred] - fac->t[idp];
            sgn = copysign(1.0, dt);
            dt = fabs(dt);
          
            /* exp(-cj * |tpred - t|) */
            exparg = exp(-cj * dt);

            /* cos(dj * |tpred - t|) = cos(dj * (tpred - t))
               = cos(dj * tpred) * cos(dj * t) + sin(dj * tpred) * sin(dj * t) */
            sum += aj * exparg * (cosarg[offpred+jkern] * fac->v[off+2*jkern] + sinarg[offpred+jkern] * fac->v[off+2*jkern+1]);

            /* sin(dj * |tpred - t|) = sgn(tpred - t) * sin(dj * (tpred - t))
               noting sin(dj * (tpred - t))
               = sin(dj * tpred) * cos(dj * t) - cos(dj * tpred) * sin(dj * t) */
            sum += bj * exparg * sgn * (sinarg[offpred+jkern] * fac->v[off+2*jkern] - cosarg[offpred+jkern] * fac->v[off+2*jkern+1]);
          }

          ktt[idp] = sum;
        }

        /* Compute K^-1(t,t) K(t,t*) */
        fsgp_apply(fac, fg, ktt, ktp);

        /* Finally compute diagonal elements of covariance
           = K(t*,t*) - K(t,t*) ktp 
           = sumaj - K(t*,t)^T ktp */
        sum = 0;
      
        for(idp = 0; idp < fac->ndp; idp++)
          sum += ktt[idp] * ktp[idp];

        varpred[ipred] = fac->sumaj - sum;
      }
    
      free((void *) ktt);
      ktt = (double *) NULL;
      free((void *) ktp);
      ktp = (double *) NULL;
    }
    
    free((void *) sinarg);
    sinarg = (double *) NULL;
    free((void *) cosarg);
    cosarg = (double *) NULL;
  
    free((void *) q);
    q = (double *) NULL;
  }
  else {  /* special case predicting at the points used in the GP */
    /* Check there aren't too many data points requested.
       User is allowed to request less. */
    if(npred > fac->ndp)
      goto error;
    
    /* mu = y - yvar * z */
    for(ipred = 0; ipred < npred; ipred++)
      ypred[ipred] = y[ipred] - fac->yvar[ipred] * z[ipred];

    /* Variance, if requested */
    if(varpred) {
      ztmp = (double *) malloc(fac->ndp * sizeof(double));
      if(!ztmp)
        goto error;

      for(ipred = 0; ipred < npred; ipred++) {
        /* K^-1 I yvar */
        kv = fsgp_apply_diag(fac, fg, ztmp, fac->yvar[ipred], ipred);

        /* var = yvar - yvar (K^-1 I yvar) */
        varpred[ipred] = fac->yvar[ipred] * (1.0 - kv);
      }
      
      free((void *) ztmp);
      ztmp = (double *) NULL;
    }
  }
  
  free((void *) z);
  free((void *) fg);
  
  return(0);

 error:
  if(z)
    free((void *) z);
  if(fg)
    free((void *) fg);
  if(sinarg)
    free((void *) sinarg);
  if(cosarg)
    free((void *) cosarg);
  if(q)
    free((void *) q);
  if(ktt)
    free((void *) ktt);
  if(ktp)
    free((void *) ktp);
  if(ztmp)
    free((void *) ztmp);
  
  return(-1);
}

/* Return weighted residuals.
   Apply Cholesky decomposition of kernel, such that
   K = C C^T, or C = L D^1/2
   then
   y = C z
   and we can solve for z by forward substitution.  The quantity z is
   what is needed when using non-linear least-squares, the weighted
   residual, where the "model part" of the log likelihood is half of
   the sum of squares of this quantity.  This separate routine can be
   used to implement non-linear least-squares fitting with a GP if the
   hyperparameters are fixed. */

int fsgp_residual (struct fsgp_fac *fac,
                   double *y, double *z, int nrhs) {
  double *fg = (double *) NULL;

  int irhs, idp, offp, off, jckern;
  double *thisy, *thisz, sum;

  /* Get out of here if there's no work to do */
  if(nrhs < 1)
    return(0);

  /* Allocate workspace */
  fg = (double *) malloc(fac->nckern * sizeof(double));
  if(!fg)
    goto error;

  for(irhs = 0; irhs < nrhs; irhs++) {
    thisy = y + irhs * fac->ndp;
    thisz = z + irhs * fac->ndp;

    /* Forward substitution */
    for(jckern = 0; jckern < fac->nckern; jckern++)
      fg[jckern] = 0;
    
    thisz[0] = thisy[0];
    
    for(idp = 1; idp < fac->ndp; idp++) {
      offp = (idp-1) * fac->nckern;
      off = idp * fac->nckern;
      
      sum = 0;
      
      for(jckern = 0; jckern < fac->nckern; jckern++) {
        fg[jckern] = fac->phi[off+jckern] * (fg[jckern] + fac->w[offp+jckern] * thisz[idp-1]);
        sum += fac->u[off+jckern] * fg[jckern];
      }
      
      thisz[idp] = (thisy[idp] - sum);
    }

    /* Divide by square root of diagonal */
    for(idp = 0; idp < fac->ndp; idp++)
      thisz[idp] /= fac->sqrtd[idp];
  }

  free((void *) fg);

  return(0);

 error:
  if(fg)
    free((void *) fg);

  return(-1);
}

/* Sample from GP given vector of standard Gaussian deviates q,
   returns result in y for zero mean. */

int fsgp_sample (struct fsgp_fac *fac, double *q, double *y) {
  double *fg = (double *) NULL;
  int idp, offp, off, jckern;
  double sqrtd, sum;

  /* Allocate workspace */
  fg = (double *) malloc(fac->nckern * sizeof(double));
  if(!fg)
    goto error;
  
  for(jckern = 0; jckern < fac->nckern; jckern++)
    fg[jckern] = 0;
  
  sqrtd = fac->sqrtd[0];

  y[0] = sqrtd * q[0];
  
  for(idp = 1; idp < fac->ndp; idp++) {
    offp = (idp-1) * fac->nckern;
    off = idp * fac->nckern;
    
    sum = 0;

    for(jckern = 0; jckern < fac->nckern; jckern++) {
      fg[jckern] = fac->phi[off+jckern] * (fg[jckern] + fac->w[offp+jckern] * sqrtd * q[idp-1]);
      sum += fac->u[off+jckern] * fg[jckern];
    }
    
    sqrtd = fac->sqrtd[idp];

    y[idp] = y[idp] * sqrtd + sum;
  }

  free((void *) fg);

  return(0);

 error:
  if(fg)
    free((void *) fg);

  return(-1);
}

/* Compute log determinant of kernel */

double fsgp_logdet (struct fsgp_fac *fac) {
  double sum;
  int idp;

  sum = 0;
  
  for(idp = 0; idp < fac->ndp; idp++)
    sum += log(fac->d[idp]);
  
  return sum;
}

/* Calculate contribution of GP to log likelihood */

int fsgp_loglike (struct fsgp_fac *fac, double *y, double *loglike) {
  double *ztmp = (double *) NULL;
  double *fg = (double *) NULL;
  double sum, logdet, logtwopi;
  int idp;

  logtwopi = log(TWOPI);
  
  ztmp = (double *) malloc(fac->ndp * sizeof(double));
  fg = (double *) malloc(fac->nckern * sizeof(double));
  if(!ztmp || !fg)
    goto error;
  
  fsgp_apply(fac, fg, y, ztmp);
  
  logdet = fsgp_logdet(fac);
  
  sum = 0;
  for(idp = 0; idp < fac->ndp; idp++)
    sum += y[idp] * ztmp[idp];

  *loglike = -0.5*(sum + logdet + fac->ndp * logtwopi);

  free((void *) ztmp);
  free((void *) fg);
  
  return(0);
  
 error:
  if(ztmp)
    free((void *) ztmp);
  if(fg)
    free((void *) fg);

  return(-1);
}
  
void fsgp_free (struct fsgp_fac *fac) {
  free((void *) fac->kern);
  fac->kern = (double *) NULL;
  free((void *) fac->t);
  fac->t = (double *) NULL;
  free((void *) fac->yvar);
  fac->yvar = (double *) NULL;
  
  free((void *) fac->phi);
  fac->phi = (double *) NULL;
  free((void *) fac->u);
  fac->u = (double *) NULL;
  free((void *) fac->v);
  fac->v = (double *) NULL;
  free((void *) fac->w);
  fac->w = (double *) NULL;
  free((void *) fac->d);
  fac->d = (double *) NULL;
  free((void *) fac->sqrtd);
  fac->sqrtd = (double *) NULL;
}
