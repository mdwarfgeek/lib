#include <math.h>

/* For FORTRAN, compile with -DFORTRAN and call as if this was defined:
   SUBROUTINE DSINCOS (a, s, c)
   DOUBLE PRECISION a, s, c
 */

#ifdef FORTRAN
void dsincos_ (double *ain, double *s, double *c) {
#else
void dsincos (double a, double *s, double *c) {
#endif
  int k, r;
  double sr, cr;

#ifdef FORTRAN
  double a = *ain;  /* copy so our mods don't affect original */
#endif

  /* Check for NaN and infinity, and act accordingly */
  if(isnan(a) || isinf(a)) {
    *s = a-a;
    *c = a-a;
  }
  else {
    /* Reduce argument to a - k*pi/2 in [-pi/4, pi/4] */
    if(a >= 0) {
      k = (int) (a / M_PI_2 + 0.5);
      r = k % 4;
    }
    else {
      k = (int) (a / M_PI_2 - 0.5);
      r = k % 4;
      if(r < 0)
	r += 4;
    }

    a -= k*M_PI_2;

    /* Compute answer */
    __asm__ volatile("fsincos"
		     : "=t" (cr), "=u" (sr)
		     : "0" (a));
    
    /* Now figure out what to do with it based on k mod 4 */
    switch(r) {
    case 0:
      *s =  sr;
      *c =  cr;
      break;
    case 1:
      *s =  cr;
      *c = -sr;
      break;
    case 2:
      *s = -sr;
      *c = -cr;
      break;
    case 3:
      *s = -cr;
      *c =  sr;
      break;
    }
  }
}

