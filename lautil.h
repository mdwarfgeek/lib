#ifndef LAUTIL_H
#define LAUTIL_H

#include "lfa.h"

int svdsolve (double *a, double *b, int m, double scale);
int svdsolcov (double *a, double *b, int m, double scale);
int mvgaussdev (struct rng_state *s,
                double *mean, double *cov,
                double *work,
                double *ans, int n);

#endif  /* LAUTIL_H */
