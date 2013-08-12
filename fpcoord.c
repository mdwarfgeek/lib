#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "astro.h"
#include "cvtunit.h"
#include "util.h"

#define MIN_SCL (FLT_RADIX*DBL_EPSILON)

int vec2tp (double s[3], double tp[3], double *x, double *y) {
  double cdtpsq, ccc, st, fac;

  /* cos(delta_tp)**2 */
  cdtpsq = tp[0]*tp[0] + tp[1]*tp[1];

  if(cdtpsq > DBL_MIN) {
    /* cos(delta)*cos(delta_tp)*cos(alpha-alpha_tp) */
    ccc = s[0]*tp[0] + s[1]*tp[1];
    
    /* Compute sin(theta) */
    st = ccc + s[2]*tp[2];

    /* Bounds check */
    if(st < MIN_SCL)
      return(-1);  /* get out now */
    
    /* Gnomonic projection: R_theta = cot theta */
    fac = 1.0 / (sqrt(cdtpsq) * st);

    (*x) = (s[1]*tp[0] - s[0]*tp[1]) * fac;
    (*y) = (s[2]*cdtpsq - tp[2]*ccc) * fac;
  }
  else {
    /* Pole, assume RA of tangent point is zero */
    st = s[2]*tp[2];

    if(st < MIN_SCL)
      return(-1);

    (*x) =  s[1] / st;
    (*y) = -s[0] / s[2];
  }

  return(0);
}

void tp2vec (double x, double y, double tp[3], double s[3]) {
  double cdtpsq, cdtp, secdtp, fac;

  /* cos(delta_tp)^2 */
  cdtpsq = tp[0]*tp[0] + tp[1]*tp[1];

  /* sin(theta) */
  fac = 1.0 / sqrt(1.0 + x*x + y*y);

  if(cdtpsq > DBL_MIN) {
    cdtp = sqrt(cdtpsq);
    secdtp = 1.0 / cdtp;

    s[0] = (tp[0] - (x*tp[1] + y*tp[0]*tp[2])*secdtp) * fac;
    s[1] = (tp[1] + (x*tp[0] - y*tp[1]*tp[2])*secdtp) * fac;
    s[2] = (tp[2] + y*cdtp) * fac;
  }
  else {
    /* Pole, assume RA of tangent point is zero */
    s[0] = -y * tp[2] * fac;
    s[1] =          x * fac;
    s[2] =      tp[2] * fac;
  }
}

