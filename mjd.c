#include <stdlib.h>
#include <math.h>

#include "lfa.h"

int date2mjd (int yr, int mn, int dy) {
  int a, m, rv;

  a = yr - (12 - mn) / 10;
  m = (mn + 9) % 12;

  rv  = (1461 * (a + 4712)) / 4;
  rv += (306 * m + 5) / 10;
  rv += dy - 2399904;
  rv -= (3 * ((a + 4900) / 100)) / 4;

  return(rv);
}

void mjd2date (int n, int *yr, int *mn, int *dy) {
  int a, b, dtmp;

  /* Convert to Julian day number */
  n += 2400001;

  /* Convert to 4*Gregorian */
  a = (4*n - 17918) / 146097;
  b = (3*a + 2) / 4;

  n = 4 * (n + b - 37);

  dtmp = 10 * (((n - 237) % 1461) / 4) + 5;

  *yr = n / 1461 - 4712;
  *mn = (dtmp / 306 + 2) % 12 + 1;
  *dy = (dtmp % 306) / 10 + 1;
}
