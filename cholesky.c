#include <stdlib.h>
#include <math.h>

#include "lfa.h"

/* In-place Cholesky decomposition for L in A = L L^T where L is a
   lower triangular matrix.  Uses the Cholesky-Crout algorithm,
   based on Stoer & Bulirsch "Introduction to numerical analysis"
   Sect. 4.3, but with the loops for i = j and i != j separated and
   changed to access only the lower triangular part of the matrix.
   Output is computed column by column so the accesses in the
   innermost loop (indexed by k) are sequential in memory. */

int cholesky (double *a, int n) {
  double x, r;
  int i, j, k;

  /* Loop over columns */
  for(j = 0; j < n; j++) {
    /* i = j */
    x = a[j*n+j];  /* A_jj */

    for(k = 0; k < j; k++)
      x -= a[j*n+k] * a[j*n+k];  /* L_jk L_jk */

    if(x < 0)
      return(-1);

    x = sqrt(x);

    a[j*n+j] = x;  /* L_jj */
    r = 1.0 / x;

    /* i != j */
    for(i = j+1; i < n; i++) {
      x = a[i*n+j];  /* A_ij */

      for(k = 0; k < j; k++)
        x -= a[i*n+k] * a[j*n+k];  /* L_ik L_ij */

      a[i*n+j] = x * r;  /* L_ij = x / L_jj */
    }
  }

  return(0);
}

