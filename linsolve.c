#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "lfa.h"
#include "util.h"

/* Solve linear system A x = b with A square (e.g. normal equations).
   Uses QR decomposition of A with column pivoting, Householder
   transformations.  During the calculation, the upper triangular
   matrix R is stored in the upper triangular part of the input.  Q^T
   b is accumulated as we go so there is no need to store Q, and the
   singular values are also calculated and stored so we can eliminate
   ill-conditioned rows during back substitution.  The procedure is
   equivalent to (but more economical than) using SVD with truncation
   of small singular values.

   This simple implementation is intended for small, light duty
   applications only.  For heavy duty applications, use LAPACK.  It
   is intended as a hopefully slightly more robust replacement for
   the original dsolve routine (which used Gaussian elimination). */

int linsolve (double *a, double *b, int n, double rcond) {
  double s[n], svsq[n];
  int perm[n];

  int i, j, c, cpiv;
  double spiv, tmp;

  double x, u, beta, sigma, mu;
  double sum;
  double diag, thresh;

  int rank;

  /* Supply default conditioning ratio if needed */
  if(rcond < 0)
    rcond = DBL_EPSILON;

  /* Compute sum of squares in each column */
  for(j = 0; j < n; j++) {
    sum = 0;

    for(i = 0; i < n; i++)
      sum += a[i*n+j] * a[i*n+j];

    s[j] = sum;
  }

  /* Loop over start column */
  thresh = 0;

  for(c = 0; c < n; c++) {
    /* Select pivot: largest s in remaining part */
    cpiv = c;
    spiv = s[c];
    for(j = c+1; j < n; j++)
      if(s[j] > spiv) {
        cpiv = j;
        spiv = s[j];
      }

    perm[c] = cpiv;

    if(cpiv != c) {
      /* Swap pivot into current column */
      for(i = 0; i < n; i++) {
        SWAP(a[i*n+cpiv], a[i*n+c], tmp);
      }

      SWAP(s[cpiv], s[c], tmp);
    }

    /* Householder transformation
       A -> A - beta v v^T A
       based on Golub & van Loan "Matrix Computations" 5.4.1. */

    /* First element */
    x = a[c*n+c];

    /* Sum of squares down column excluding x */
    sigma = 0;

    for(i = c+1; i < n; i++)
      sigma += a[i*n+c] * a[i*n+c];

    /* If norm is zero, we can skip the transformation */
    if(sigma > 0) {
      /* Norm */
      mu = sqrt(x*x + sigma);

      /* First element of Householder vector */
      if(x > 0)
        u = -sigma / (x + mu);
      else
        u = x - mu;
      
      beta = 2.0 / (u*u + sigma);
      
      /* A_cc = A_cc - beta v_c v_k A_kc
              = A_cc - beta v_c (v_c A_cc + A_c+1,c^2 + ...) */
      a[c*n+c] -= beta * u * (u*x + sigma);

      /* A_ij = A_ij - beta v_i v_k A_kj */
      for(j = c+1; j < n; j++) {
        sum = u * a[c*n+j];  /* v_k A_kj, k = c */
      
        for(i = c+1; i < n; i++)
          sum += a[i*n+c] * a[i*n+j];  /* v_k A_kj */
      
        sum *= beta;

        a[c*n+j] -= u * sum;  /* v_i sum, i = c */

        for(i = c+1; i < n; i++)
          a[i*n+j] -= a[i*n+c] * sum;  /* v_i sum */

        /* Update s */
        s[j] -= a[c*n+j] * a[c*n+j];
      }

      /* b_i = b_i - beta v_i v_k b_k */
      sum = u * b[c];  /* v_k b_k, k = c */

      for(i = c+1; i < n; i++)
        sum += a[i*n+c] * b[i];  /* v_k b_k */

      sum *= beta;

      b[c] -= u * sum;  /* v_i * sum, i = c */

      for(i = c+1; i < n; i++)
        b[i] -= a[i*n+c] * sum;  /* v_i * sum */
    }

    /* Compute square of singular value, this is just the norm^2 of
       the output row. */
    sum = 0;

    for(j = c; j < n; j++)
      sum += a[c*n+j] * a[c*n+j];

    svsq[c] = sum;

    /* Maintain maximum */
    if(sum > thresh)
      thresh = sum;
  }

  /* Set threshold in square of singular values */
  thresh *= rcond;
  thresh *= rcond;

  /* Compute R^-1 (Q^t b) by back substitution */
  rank = 0;

  for(i = n-1; i >= 0; i--) {
    diag = a[i*n+i];

    if(svsq[i] > thresh) {
      sum = b[i];
      for(j = i+1; j < n; j++)
        sum -= a[i*n+j] * b[j];
      b[i] = sum / diag;

      rank++;
    }
    else
      b[i] = 0;
  }

  /* Unwrap the permutations applied during pivoting */
  for(i = n-1; i >= 0; i--)
    if(perm[i] != i) {
      SWAP(b[perm[i]], b[i], tmp);
    }

  return(rank);
}

