/* lautil.c:  Various linear algebra related utility routines that depend
              on LAPACK. */

#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "lfa.h"
#include "lautil.h"

/* Prototypes for LAPACK Fortran routines */
void dgelss_ (int *, int *, int *,
              double *, int *, double *, int *, double *, double *, int *,
              double *, int *, int *);
void dgesvd_ (char *, char *,
	      int *, int *, double *, int *,
	      double *, double *, int *, double *, int *,
	      double *, int *, int *, unsigned int, unsigned int);
void dpotrf_ (char *, int *, double *, int *, int *, unsigned int);

/* Computes minimum norm solution only for a linear least squares
   problem.  Simple n=1 driver for the LAPACK dgelss routine. */

int svdsolve (double *a, double *b, int m, double scale) {
  int lwork = 6*m;
  double rcond = -1.0;
  int nrhs = 1, rank, info;

  VLAONSTACK(double, s, m);
  VLAONSTACK(double, work, lwork);

  dgelss_(&m, &m, &nrhs,
          a, &m, b, &m, s, &rcond, &rank,
          work, &lwork, &info);

  return(info);
}

/* Computes minimum norm solution and covariance matrix for a linear
   least squares problem. */

int svdsolcov (double *a, double *b, int m, double scale) {
  int lwork = 6*m;
  int mthr, info;

  int i, j, k;
  double thresh, sum, sol;

  VLAONSTACK(double, s, m);
  VLAONSTACK(double, u, m*m);
  VLAONSTACK(double, vt, m*m);
  VLAONSTACK(double, work, lwork);

  /* Perform decomposition */
  dgesvd_("A", "A",
	  &m, &m,
	  a, &m, s, u, &m, vt, &m,
	  work, &lwork, &info, 1, 1);

  if(scale < 0)
    /* Machine precision */
    scale = FLT_RADIX * DBL_EPSILON;

  /* Threshold for singular values */
  thresh = scale*s[0];
  if(thresh < DBL_MIN)  /* smallest positive number */
    thresh = DBL_MIN;

  /* Take reciprocal of elements of s above threshold */
  for(k = 0; k < m && s[k] > thresh; k++)
    s[k] = 1.0 / s[k];

  mthr = m;

  /* Compute inverse and solution vector */
  for(i = 0; i < m; i++) {
    sol = 0;

    for(j = 0; j < m; j++) {
      /* i'th row of v multiplied by j'th column of u transpose */
      sum = 0;
      for(k = 0; k < mthr; k++)
	sum += vt[i*m+k]*u[k*m+j] * s[k];

      a[i*m+j] = sum;

      sol += b[j]*sum;
    }

    work[i] = sol;
  }    

  for(i = 0; i < m; i++)
    b[i] = work[i];

  return(info);
}

/* Multivariate Gaussian deviates.  User-passed workspace can be used to
   store Cholesky decomposition of covariance matrix if multiple sets of
   deviates need to be generated for the same covariance. */

int mvgaussdev (struct rng_state *s,
                double *mean, double *cov,
                double *work,
                double *ans, int n) {
  int i, j, info;
  double zz;

  VLAONSTACK(double, z, n);

  if(cov) {
    /* Copy into workspace */
    for(j = 0; j < n; j++)
      for(i = 0; i < n; i++)
        work[j*n+i] = cov[j*n+i];
    
    /* Cholesky decomposition of the covariance matrix */
    dpotrf_("L", &n, work, &n, &info, 1);
  }
  /* otherwise, assume Cholesky is already done */

  if(ans) {
    /* ans = mean + work * Z where Z is a vector of Gaussian deviates */
    
    /* Start off with the means */
    for(j = 0; j < n; j++)
      ans[j] = mean[j];
    
    /* Generate Gaussian deviates */
    rng_fetch_gauss(s, z, n);
    
    /* Loop over columns of the lower triangular matrix */
    for(i = 0; i < n; i++) {
      zz = z[i];
      
      /* Now add into the answer */
      for(j = i; j < n; j++)
        ans[j] += zz*work[i*n+j];
    }
  }

  return(info);
}
