#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lfa.h"
#include "util.h"

/* Compute air mass using Hardie (1962) formula */

double v_airmass (double v[3]) {
  double z, secz, szm;

  /* Hold constant below about 3 degrees */
  if(v[2] < 0.05)
    z = 0.05;
  else
    z = v[2];

  secz = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) / z;
  szm = secz - 1.0;

  return(secz - ((0.0008083*szm + 0.002875) * szm + 0.0018167) * szm);
}

