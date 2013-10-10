#include <stdlib.h>
#include <math.h>

#include "lfa.h"

int date2mjd (int yr, int mn, int dy) {
  int a, m, rv;

  a = yr - (12 - mn) / 10;
  m = (mn - 3) % 12;

  rv  = 365.25 * (a + 4712);
  rv += 30.6 * m + 0.5;
  rv += dy - 2399904;
  rv += (a / 100 + 49) * 0.75;

  return(rv);
}

