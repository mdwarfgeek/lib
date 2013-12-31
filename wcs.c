#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "lfa.h"

#include "cvtunit.h"
#include "util.h"

#define TINY (FLT_RADIX*DBL_EPSILON)
#define SAFE 1.0e-10

/* Notes: see Calabretta & Greisen (2002) Sect. 2.3 and 5.1 for
   equations.  We only handle the default LATPOLE and LONPOLE
   values, of LATPOLE = theta_0 = pi/2, LONPOLE = phi_0 = pi
   so sin(phi - phi_p) = -sin(phi), cos(phi - phi_p) = -cos(phi)
   ZPN with zero-order terms in r is not supported, it blows up
   at the origin and would need a special implementation. */

void wcs_vec2xy (struct wcs_info *wcs, double *vec, double *x, double *y) {
  double camtcd, ctsp, ctcp, rt, st, fx = 0, fy = 0, rfac, denom;
  double ctsq, tt, fac;
  int i;

  /* cos(alpha-alpha_0)*cos(delta) */
  camtcd = vec[0]*wcs->cosa + vec[1]*wcs->sina;

  /* cos(theta)*sin(phi-phi_p) = cos(delta)*sin(alpha-alpha_0) */
  ctsp = vec[1]*wcs->cosa - vec[0]*wcs->sina;

  /* cos(theta)*cos(phi-phi_p) = cos(delta)*sin(delta_0)*cos(alpha-alpha_0)
                                 - sin(delta)*cos(delta_0) */
  ctcp = camtcd*wcs->sind - vec[2]*wcs->cosd;

  /* sin(theta) = sin(delta)*sin(delta_0)
                  + cos(delta)*cos(delta_0)*cos(alpha-alpha_0) */
  st = vec[2]*wcs->sind + camtcd*wcs->cosd;

  /* Convert to appropriate focal plane angular coordinates for proj. */
  switch(wcs->proj) {
  case PROJ_TAN:
    /* Gnomonic projection: R_theta = cot theta */
    if(st < TINY)
      tt = 0;  /* arbitrary, but no real image should do this! */
    else
      tt = 1.0 / st;

    fx = ctsp*tt;
    fy = -ctcp*tt;
    break;
  case PROJ_SIN:
    /* Orthographic projection: R_theta = cos theta */
    fx = ctsp;
    fy = -ctcp;
    break;
  case PROJ_ARC:
    /* ARC projection: R_theta = pi/2 - theta => cos R_theta = sin theta */

    /* Roundoff error can sometimes push st slightly out of range, so supply
       the appropriate values at the limits in case this happens. */
    if(st < -1.0)
      rt = M_PI;
    else if(st > 1.0)
      rt = 0;
    else
      rt = acos(st);

    /* cos^2(theta) = sin^2(rt) */
    ctsq = ctsp*ctsp + ctcp*ctcp;

    /* Is there any distortion? */
    if(wcs->mpc > 0) {
      /* Evaluate polynomial */
      rfac = wcs->pc[wcs->mpc];
      for(i = wcs->mpc-1; i > 0; i--)
	rfac = rfac * rt + wcs->pc[i];
    }
    else
      rfac = 1.0;

    /* Handle case of origin */
    if(ctsq < SAFE) {  /* use series rt / sin(rt) to O(rt^4) */
      tt = rt*rt;
      fac = ((7.0/360.0) * tt + (1.0/6.0)) * tt + 1.0;
    }
    else  /* actual expression: rt / cos(theta) */
      fac = rt / sqrt(ctsq);

    fx = rfac*fac*ctsp;
    fy = -rfac*fac*ctcp;
    break;
  }

  /* Now transform to pixels */
  denom = wcs->a * wcs->e - wcs->d * wcs->b;

  *x = (fx * wcs->e - fy * wcs->b) / denom + wcs->c;
  *y = (fy * wcs->a - fx * wcs->d) / denom + wcs->f;
}

void wcs_xy2vec (struct wcs_info *wcs, double x, double y, double *vec) {
  double fx, fy, st = 0, ctsp = 0, ctcp = 0, rtsq, tt, ttsq, rfac, drfac, ct, camtcd;
  double fac, rt, dtt;
  int i, iter;

  /* First transform from pixels to focal plane angular coords. */
  x -= wcs->c;
  y -= wcs->f;

  fx = wcs->a * x + wcs->b * y;
  fy = wcs->d * x + wcs->e * y;

  /* Convert to the quantities cos(theta)*sin(phi), cos(theta)*cos(phi), and
   * sin(theta) needed to do the inversion to alpha,delta.  Case for each proj.
   */
  switch(wcs->proj) {
  case PROJ_TAN:
    /* Gnomonic projection: R_theta = cot theta */
    st = sqrt(1.0/(1.0+fx*fx+fy*fy));
    ctsp = fx*st;
    ctcp = -fy*st;
    break;
  case PROJ_SIN:
    /* Orthographic projection: R_theta = cos theta */
    ctsp = fx;
    ctcp = -fy;
    st = sqrt(1.0-fx*fx-fy*fy);
    break;
  case PROJ_ARC:
    /* ARC projection: R_theta = pi/2 - theta => cos R_theta = sin theta */
    rtsq = fx*fx+fy*fy;
    rt = sqrt(rtsq);

    /* Is there any distortion? */
    if(wcs->mpc > 0) {
      /* Invert polynomial to get tt = 90-theta */
      tt = rt;

      for(iter = 0; iter < 10; iter++) {
	rfac = wcs->pc[wcs->mpc];
	drfac = wcs->mpc * wcs->pc[wcs->mpc];
	for(i = wcs->mpc-1; i > 0; i--) {
	  rfac = rfac * tt + wcs->pc[i];
	  drfac = drfac * tt + i * wcs->pc[i];
	}

	dtt = (tt*rfac - rt) / drfac;

	if(fabs(dtt) < TINY)
	  break;

	tt -= dtt;
      }

      dsincos(tt, &ct, &st);  
      ttsq = tt*tt;

      if(ttsq < SAFE)  /* tt/rt * series sin(tt)/tt to O(rt^4) */
	fac = (((1.0/120.0) * ttsq - (1.0/6.0)) * ttsq + 1.0) / rfac;
      else  /* actual expression: cos(theta) / rt */
	fac = ct / rt;
    }
    else {
      dsincos(rt, &ct, &st);  
      
      if(rtsq < SAFE)  /* use series for sin(rt) / rt to O(rt^4) */
	fac = ((1.0/120.0) * rtsq - (1.0/6.0)) * rtsq + 1.0;
      else  /* actual expression: cos(theta) / rt */
	fac = ct / rt;
    }

    ctsp = fx*fac;
    ctcp = -fy*fac;

    break;
  }

  /* sin(alpha-alpha_0)*cos(delta) = -cos(theta)*sin(phi-phi_0) = ctsp */
  /* cos(alpha-alpha_0)*cos(delta) = sin(theta)*cos(delta_0) 
                                     - cos(theta)*sin(delta_0)*cos(phi-phi_p) */
  camtcd = st*wcs->cosd + ctcp*wcs->sind;  /* ctcp = -cos(theta)*cos(phi-phi_p) */
  
  /* cos(alpha)*cos(delta) = cos((alpha-alpha_0)+alpha_0)*cos(delta) */
  vec[0] = camtcd*wcs->cosa - ctsp*wcs->sina;

  /* sin(alpha)*cos(delta) = sin((alpha-alpha_0)+alpha_0)*cos(delta) */
  vec[1] = ctsp*wcs->cosa + camtcd*wcs->sina;

  /* sin(delta) = sin(theta)*sin(delta_0) + cos(theta)*cos(delta_0)*cos(phi-phi_p) */
  vec[2] = st*wcs->sind - ctcp*wcs->cosd;  /* ctcp = -cos(theta)*cos(phi-phi_p) */
}

void wcs_ad2xy (struct wcs_info *wcs, double a, double d,
		double *x, double *y) {
  double sina, cosa, sind, cosd, vec[3];

  dsincos(a, &sina, &cosa);
  dsincos(d, &sind, &cosd);

  vec[0] = cosa*cosd;
  vec[1] = sina*cosd;
  vec[2] = sind;

  wcs_vec2xy(wcs, vec, x, y);
}

void wcs_xy2ad (struct wcs_info *wcs, double x, double y, double *a, double *d) {
  double vec[3];

  wcs_xy2vec(wcs, x, y, vec);

  /* Finally, back to RA/DEC */
  *a = atan2(vec[1], vec[0]);
  if((*a) < 0.0)
    *a += TWOPI;

  *d = atan2(vec[2], sqrt(vec[0]*vec[0]+vec[1]*vec[1]));
}

void wcs_xy2xy (struct wcs_info *wcs1, struct wcs_info *wcs2,
		double x1, double y1, double *x2, double *y2) {
  double vec[3];

  wcs_xy2vec(wcs1, x1, y1, vec);
  wcs_vec2xy(wcs2, vec, x2, y2);
}

void wcs_tp2xy (struct wcs_info *wcs, double fx, double fy,
		double *x, double *y) {
  double cotsq, cot, rt, st, rfac, denom;
  double tt, fac;
  int i;

  /* Convert to appropriate focal plane angular coordinates for proj. */
  switch(wcs->proj) {
  case PROJ_TAN:
    /* Already done */
    break;
  case PROJ_SIN:
    /* Gnomonic to Orthographic: multiply by sin theta */
    st = sqrt(1.0/(1.0+fx*fx+fy*fy));
    fx *= st;
    fy *= st;
    break;
  case PROJ_ARC:
    /* cot^2 theta */
    cotsq = fx*fx+fy*fy;
    cot = sqrt(cotsq);

    /* R_theta = pi/2 - theta => tan R_theta = cot theta */
    rt = atan(cot);

    /* Is there any distortion? */
    if(wcs->mpc > 0) {
      /* Evaluate polynomial */
      rfac = wcs->pc[wcs->mpc];
      for(i = wcs->mpc-1; i > 0; i--)
	rfac = rfac * rt + wcs->pc[i];

      fx *= rfac;
      fy *= rfac;
    }

    /* Handle case of origin */
    if(cotsq < SAFE) {  /* use series rt / tan(rt) to O(rt^6) */
      tt = rt*rt;
      fac = (((-2.0/945.0) * tt + (1.0/45.0)) * tt - (1.0/3.0)) * tt + 1.0;
    }
    else  /* actual expression: rt / cot(theta) */
      fac = rt / cot;

    fx *= fac;
    fy *= fac;
    break;
  }

  /* Now transform to pixels */
  denom = wcs->a * wcs->e - wcs->d * wcs->b;

  *x = (fx * wcs->e - fy * wcs->b) / denom + wcs->c;
  *y = (fy * wcs->a - fx * wcs->d) / denom + wcs->f;
}

void wcs_xy2tp (struct wcs_info *wcs, double x, double y,
		double *fx, double *fy) {
  double rtsq, csct, tt, ttsq, rfac, drfac;
  double fac, rt, dtt;
  int i, iter;

  /* First transform from pixels to focal plane angular coords. */
  x -= wcs->c;
  y -= wcs->f;

  (*fx) = wcs->a * x + wcs->b * y;
  (*fy) = wcs->d * x + wcs->e * y;

  /* Convert to Gnomonic */
  switch(wcs->proj) {
  case PROJ_TAN:
    /* Already done */
    break;
  case PROJ_SIN:
    /* Orthographic to Gnomonic: divide by sin theta */
    csct = sqrt(1.0/(1.0-(*fx)*(*fx)-(*fy)*(*fy)));
    (*fx) *= csct;
    (*fy) *= csct;
    break;
  case PROJ_ARC:
    /* Extract R_theta */
    rtsq = (*fx)*(*fx)+(*fy)*(*fy);
    rt = sqrt(rtsq);

    /* Is there any distortion? */
    if(wcs->mpc > 0) {
      /* Invert polynomial to get tt = 90-theta */
      tt = rt;

      for(iter = 0; iter < 10; iter++) {
	rfac = wcs->pc[wcs->mpc];
	drfac = wcs->mpc * wcs->pc[wcs->mpc];
	for(i = wcs->mpc-1; i > 0; i--) {
	  rfac = rfac * tt + wcs->pc[i];
	  drfac = drfac * tt + i * wcs->pc[i];
	}

	dtt = (tt*rfac - rt) / drfac;

	if(fabs(dtt) < TINY)
	  break;

	tt -= dtt;
      }

      ttsq = tt*tt;

      if(ttsq < SAFE)  /* tt/rt * series tan(tt)/tt to O(rt^6) */
	fac = ((((17.0/315.0) * ttsq + (2.0/15.0)) * ttsq + (1.0/3.0)) * ttsq + 1.0) / rfac;
      else  /* actual expression: cot(theta) / rt = tan(tt) / rt */
	fac = tan(tt) / rt;
    }
    else {
      if(rtsq < SAFE)  /* use series for tan(rt) / rt to O(rt^6) */
	fac = (((17.0/315.0) * rtsq + (2.0/15.0)) * rtsq + (1.0/3.0)) * rtsq + 1.0;
      else  /* actual expression: cot(theta) / rt = tan(rt) / rt */
	fac = tan(rt) / rt;
    }

    (*fx) *= fac;
    (*fy) *= fac;

    break;
  }
}
