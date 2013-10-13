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

