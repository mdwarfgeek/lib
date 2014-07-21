#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lfa.h"

/* Geodetic coordinates to geocentric, WGS84.
   The result of this routine can be converted to a 3-vector as follows:
   x = u * cos(lambda)
   y = u * sin(lambda)
*/

void geoc (double sinphi,   /* sin, cos geodetic latitude */
	   double cosphi,
	   double height,   /* above geoid, m */
	   double *u,       /* returned: distance from Earth spin axis, m */
	   double *z) {     /* returned: distance N from equatorial plane, m */
  double cfsq, n;

  /* Square of compression factor: b^2/a^2 = (1 - f)^2 */
  cfsq = (1.0 - FEARTH);
  cfsq *= cfsq;

  /* Distance from the centre * b/a */
  n = AEARTH / sqrt(cosphi*cosphi + cfsq * sinphi*sinphi);

  /* Distance from Earth axis */
  (*u) = (n       + height) * cosphi;
  (*z) = (n*cfsq  + height) * sinphi;
}

