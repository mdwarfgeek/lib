#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "lfa.h"
#include "util.h"

/* Default rank conditioning criterion */
#define DEFAULT_RCOND (FLT_RADIX * DBL_EPSILON)

/* QR decomposition of A with column pivoting, Householder
   transformations.  During the calculation, the upper triangular
   matrix R is stored in the upper triangular part of the input.
   Singular values are also calculated and stored.

   This simple implementation is intended for small, light duty
   applications only.  For heavy duty applications, use LAPACK. */

void qr (double *a, double *s, double *betaarr, int *perm, int n) {
  int i, j, c, cpiv;
  double spiv, tmp;

  double x, u, beta, sigma, mu;
  double sum;

  /* Compute sum of squares in each column */
  for(j = 0; j < n; j++) {
    sum = 0;

    for(i = 0; i < n; i++)
      sum += a[i*n+j] * a[i*n+j];

    s[j] = sum;
  }

  /* Loop over start column */
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

      /* Value of beta for first element = unity */
      tmp = u*u;
      beta = 2.0 * tmp / (tmp + sigma);

      /* Convert Householder vector to have first element = unity */
      tmp = 1.0 / u;
      for(i = c+1; i < n; i++)
        a[i*n+c] *= tmp;
      
      /* A_cc = A_cc - beta v_c v_k A_kc
              = A_cc - beta v_c (v_c A_cc + A_c+1,c^2 + ...) */
      a[c*n+c] -= beta * tmp * (u*x + sigma);

      /* A_ij = A_ij - beta v_i v_k A_kj */
      for(j = c+1; j < n; j++) {
        sum = a[c*n+j];  /* v_k A_kj, k = c */
      
        for(i = c+1; i < n; i++)
          sum += a[i*n+c] * a[i*n+j];  /* v_k A_kj */
      
        sum *= beta;

        a[c*n+j] -= sum;  /* v_i sum, i = c */

        for(i = c+1; i < n; i++)
          a[i*n+j] -= a[i*n+c] * sum;  /* v_i sum */

        /* Update s */
        s[j] -= a[c*n+j] * a[c*n+j];
      }

      betaarr[c] = beta;
    }
    else
      betaarr[c] = 0;

    /* Compute square of singular value, this is just the norm^2 of
       the output row.  We can store it in s because we're done
       with the 'c'th element. */
    sum = 0;

    for(j = c; j < n; j++)
      sum += a[c*n+j] * a[c*n+j];

    s[c] = sum;
  }
}

/* Solve linear system A x = b with A square (e.g. normal equations).
   x is returned in the array for b.  Uses QR routine above.  Singular
   values are used to eliminate ill-conditioned rows during back
   substitution.  This procedure is equivalent to (but more economical
   than) using SVD with truncation of small singular values. */

int qrsolve (double *a, double *s, double *betaarr, int *perm,
             double *b, int n, double rcond) {
  int i, j, c;
  double tmp;

  double beta;
  double sum;
  double diag, thresh;

  int rank;

  /* Supply default conditioning ratio if needed */
  if(rcond < 0)
    rcond = DEFAULT_RCOND;

  /* Apply Householder matrices to b */
  thresh = 0;

  for(c = 0; c < n; c++) {
    beta = betaarr[c];

    if(beta) {
      /* b_i = b_i - beta v_i v_k b_k */
      sum = b[c];  /* v_k b_k, k = c */

      for(i = c+1; i < n; i++)
        sum += a[i*n+c] * b[i];  /* v_k b_k */

      sum *= beta;

      b[c] -= sum;  /* v_i * sum, i = c */

      for(i = c+1; i < n; i++)
        b[i] -= a[i*n+c] * sum;  /* v_i * sum */
    }

    /* Maintain maximum singular value */
    if(s[c] > thresh)
      thresh = s[c];
  }

  /* Set threshold in square of singular values */
  thresh *= rcond;
  thresh *= rcond;

  if(thresh < DBL_MIN)  /* smallest positive number */
    thresh = DBL_MIN;

  /* Compute R^-1 (Q^t b) by back substitution */
  rank = 0;

  for(i = n-1; i >= 0; i--) {
    diag = a[i*n+i];

    if(s[i] > thresh) {
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

/* Pseudo-inverse using same methods as "qrsolve". */

int qrinvert (double *a, double *s, double *betaarr, int *perm,
              double *ainv, int n, double rcond) {
  int i, j, k, c;
  double *pin, *pout, tmp;
  double *qj;

  double beta;
  double sum;
  double diag, thresh;

  int rank;

  /* Supply default conditioning ratio if needed */
  if(rcond < 0)
    rcond = DEFAULT_RCOND;

  /* Form Q^T by backward accumulation */
  thresh = 0;

  for(c = n-1; c >= 0; c--) {
    beta = betaarr[c];

    if(beta) {
      /* Q_ij = Q_ij - beta v_i v_k Q_kj 
         so
         Q^T_ji = Q^T_ji - beta v_i v_k Q^T_jk 

         The first time we "see" an element, it is the element of
         the identity, so handle j=c separately. */
      
      /* j=c, Q_j is just 1, 0.... */
      qj = &(ainv[c*n]);  /* row Q^T_c */

      qj[c] = 1.0 - beta;  /* v_i * sum, i = c; qj_c = 1 initially */
      
      for(i = c+1; i < n; i++)
        qj[i] = -a[i*n+c] * beta;  /* v_i * sum, qj_i = 0 initially */

      for(j = c+1; j < n; j++) {
        qj = &(ainv[j*n]);  /* row Q^T_j */

        sum = 0;  /* v_k qj_k, k = c; qj_c = 0 */

        for(i = c+1; i < n; i++)
          sum += a[i*n+c] * qj[i];  /* v_k qj_k */
        
        sum *= beta;
        
        qj[c] = -sum;  /* v_i * sum, i = c; qj_c = 0 initially */
        
        for(i = c+1; i < n; i++)
          qj[i] -= a[i*n+c] * sum;  /* v_i * sum */
      }
    }
    else {
      /* Initialize elements of Q^T we're seeing for the first time */
      qj = &(ainv[c*n]);  /* row Q^T_c */

      qj[c] = 1.0;
      
      for(i = c+1; i < n; i++)
        qj[i] = 0.0;

      for(j = c+1; j < n; j++)
        ainv[j*n+c] = 0.0;
    }

    /* Maintain maximum singular value */
    if(s[c] > thresh)
      thresh = s[c];
  }

  /* Set threshold in square of singular values */
  thresh *= rcond;
  thresh *= rcond;

  if(thresh < DBL_MIN)  /* smallest positive number */
    thresh = DBL_MIN;

  /* Compute R^-1 Q^t by back substitution */
  rank = 0;

  for(i = n-1; i >= 0; i--) {
    diag = a[i*n+i];

    if(s[i] > thresh) {
      diag = 1.0 / diag;

      for(k = n-1; k >= 0; k--) {
        sum = ainv[i*n+k];
        for(j = i+1; j < n; j++)
          sum -= a[i*n+j] * ainv[j*n+k];
        ainv[i*n+k] = sum * diag;
      }

      rank++;
    }
    else {
      for(k = n-1; k >= 0; k--)
        ainv[i*n+k] = 0;
    }
  }

  /* Unwrap the permutations applied during pivoting */
  for(i = n-1; i >= 0; i--)
    if(perm[i] != i) {
      /* Swap rows */
      pin = &(ainv[perm[i]*n]);
      pout = &(ainv[i*n]);

      for(j = 0; j < n; j++) {
        SWAP(pin[j], pout[j], tmp);
      }
    }

  return(rank);
}

/* Simple driver for users needing only a single solve for small
   values of n that will fit on the stack. */

int linsolve (double *a, double *b, int n, double rcond) {
  VLAONSTACK(double, s, n);
  VLAONSTACK(double, betaarr, n);
  VLAONSTACK(int, perm, n);

  /* QR */
  qr(a, s, betaarr, perm, n);

  /* Solve */
  return qrsolve(a, s, betaarr, perm, b, n, rcond);
}
