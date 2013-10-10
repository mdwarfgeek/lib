#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lfa.h"
#include "cvtunit.h"
#include "util.h"

/* Vector must be in (-h, delta) system, and not adjusted for refraction. */

double parallactic (double sinphi, double cosphi, double v[3]) {
  return atan2(v[1] * sinphi, (1.0 - v[2]*v[2]) * sinphi - v[0]*v[2] * cosphi);
}
