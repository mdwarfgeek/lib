#include <stdio.h>
#include <math.h>

#include "lfa.h"

#define KEPLER_PREC    1.0e-13
#define KEPLER_MAXITER 100

double kepler (double ma, double ecc) {
  double ea, tmp, alpha, beta, z, ste;
  double se, ce, f, df, delta;
  int i;

  /* Reduce mean anomaly to [-pi, pi] */
  ma = remainder(ma, TWOPI);

  /* Eccentric or circular? */
  if(ecc > 0) {
    /* For eccentric orbits, use cubic approximation from Mikkola (1987)
       to provide initial guess of eccentric anomaly. */
    
    /* Eq. 9a */
    tmp = 1.0 / (4 * ecc + 0.5);
    
    alpha = (1.0 - ecc) * tmp;
    beta = 0.5 * ma * tmp;
    
    /* Eq. 9b */
    z = cbrt(beta + copysign(sqrt(beta*beta + alpha*alpha*alpha), beta));
    
    /* Eq. 9c: initial value of sin(E/3) */
    ste = z - alpha / z;
    
    /* Eq. 7: 5th order correction term */
    ste -= 0.078 * ste*ste*ste*ste*ste / (1.0 + ecc);
    
    /* Eq. 8: eccentric anomaly */
    ea = ma + ecc * ste * (3.0 - 4.0 * ste*ste);
    
    /* Refine solution of Kepler's equation using Newton's method */
    for(i = 0; i < KEPLER_MAXITER; i++) {
      inline_bare_sincos(ea, se, ce);
      
      f = ea - ecc * se - ma;
      df = 1.0 - ecc * ce;
      delta = f / df;
      
      ea -= delta;
      
      if(fabs(delta) < KEPLER_PREC) {
        /* I think that's enough... */
        break;
      }
    }
    
    if(i >= KEPLER_MAXITER)
      fprintf(stderr, "kepler: iteration limit reached\n");
  }
  else
    ea = ma;

  return(ea);
}
