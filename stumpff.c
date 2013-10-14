#include <math.h>

#include "lfa.h"

#define THRESH 1.0e-3

static void pdsincosh (double x, double *s, double *c);

/* Evaluates s^k c_k(x) for k=0,1,2,3 where x = alpha*s^2 */

void stumpff (double s,
	      double alpha,
	      double sqrtalpha,
	      double *c) {
  /* Coefficients in series for c_2 and c_3.  Series truncated so error
     is less than 1.0e-15 when using x in [0, 1.0e-3].  These are used
     rather than c_0 and c_1 to avoid division. */
  static const double coeff[] = {
    /* c_2                 c_3 */
    1.0/3628800,        1.0/39916800,
    1.0/40320,          1.0/362880,
    1.0/720,            1.0/5040,
    1.0/24,             1.0/120,
    1.0/2,              1.0/6
  };

  double x, srx, ssx, csx, ss, ralpha;

  /* Compute x */
  x = alpha * s*s;
  
  if(x > THRESH) {
    /* Use sin, cos of sqrt(x) */
    srx = sqrtalpha * fabs(s);

    dsincos(srx, &ssx, &csx);

    c[0] = csx;
    c[1] = ssx / sqrtalpha;
    if(s < 0)
      c[1] *= -1;

    /* Use recursion relation for the other two:
       x c_k+2(x) = 1/k! - c_k(x) */
    ralpha = 1.0 / alpha;

    c[2] = (1.0 - c[0]) * ralpha;
    c[3] = (s - c[1]) * ralpha;
  }
  else if(x < -THRESH) {
    /* Use sinh, cosh of sqrt(-x) */
    srx = sqrtalpha * fabs(s);

    pdsincosh(srx, &ssx, &csx);

    c[0] = csx;
    c[1] = ssx / sqrtalpha;
    if(s < 0)
      c[1] *= -1;

    /* Use recursion relation for the other two:
       x c_k+2(x) = 1/k! - c_k(x) */
    ralpha = 1.0 / alpha;

    c[2] = (1.0 - c[0]) * ralpha;
    c[3] = (s - c[1]) * ralpha;
  }
  else {
    /* Use series for small arguments */
    sum_poly(x, coeff, 2, 5, &(c[2]));

    ss = s*s;

    c[2] *= ss;
    c[3] *= s*ss;
    
    /* Use recursion relation for the other two:
       c_k(x) = 1/k! - x c_k+2(x) */
    c[0] = 1.0 - c[2] * alpha;
    c[1] = s - c[3] * alpha;
  }
}

/* Simultaneous sinh and cosh for positive arguments, used above.
   Needs the C99 expm1 function, but I'm assuming most systems
   have this. */

static void pdsincosh (double x, double *s, double *c) {
  double exm, ex, rex;

  /* exp(x)-1 */
  exm = expm1(x);

  /* exp(x) */
  ex = exm + 1;

  /* exp(-x) */
  rex = 1.0 / ex;

  /* sinh(x) = (exp(x) - 1 + (exp(x) - 1) / exp(x)) / 2 
     Prevents <number close to 1> - <number close to 1> */
  *s = 0.5*(exm + exm * rex);

  if(x > 0.5) {
    /* cosh(x) = (exp(x) + exp(-x)) / 2 */
    *c = 0.5*(ex + rex);
  }
  else {
    /* cosh(x) = 1.0 + (exp(x)-1)**2 / (2 * exp(x))
       Prevents <number close to 1> - <number close to 1> */
    *c = 1.0 + 0.5*exm*exm*rex;
  }
}
