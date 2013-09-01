#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lfa.h"
#include "util.h"

#define MAX_SECZ 11.4737  /* 85 degrees */

/* Compute air mass using Hardie (1962) formula */

double airmass (double secz) {
  double szm;

  /* Hold constant above 85 degrees zenith distance */
  if(secz < 0 || secz > MAX_SECZ)
    secz = MAX_SECZ;

  szm = secz - 1.0;

  return(secz - (0.0018167 +
		(0.002875 + 0.0008083*szm) * szm) * szm);
}
