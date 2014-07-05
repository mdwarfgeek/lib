#include <stdio.h>  /* pulls in __STDC_IEC_559__ in glibc */
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include "lfa.h"

#define MASKL UINT64_C(0x000ffafffffffb3f)
#define MASKH UINT64_C(0x000ffdfffc90fffd)
#define INIT_ITER(v, i) INT32_C(1812433253) * ((v) ^ ((v) >> 30)) + (i)
#define INIT_MSKH(v) ((v) & UINT32_C(0x000fffff)) | UINT32_C(0x3ff00000)

/* Random number generator using Mersenne Twister.  Based on dSFMT by
   Saito and Matsumoto, but implemented in plain C (for my applications,
   the performance is adequate as-is, without the complications of the
   original SIMD implementation).  The generator produces IEEE 754
   binary64 directly.  On non-IEEE machines we use a normal integer to
   floating point conversion on the 52-bit integers instead. */

void rng_init (struct rng_state *s, uint32_t vl) {
  uint64_t *p;
  int i, il;

#ifndef __STDC_IEC_559__
  double d = 1.0;
#endif

  p = &(s->a[0]);

  *p = vl;
  vl = INIT_ITER(vl, 1);
  *p++ |= ((uint64_t) INIT_MSKH(vl)) << 32;

  for(i = 1; i < RNG_RQ; i++) {
    il = i+i;
    
    vl = INIT_ITER(vl, il);
    *p = vl;
    vl = INIT_ITER(vl, il+1);
    *p++ |= ((uint64_t) INIT_MSKH(vl)) << 32;
  }

  for(; i < RNG_NQ; i++) {
    il = i+i;

    vl = INIT_ITER(vl, il);
    *p = vl;
    vl = INIT_ITER(vl, il+1);
    *p++ |= ((uint64_t) vl) << 32;
  }

  s->p = (uint64_t *) NULL;
  s->r = 0;

  s->gauss = 0;
  s->havegauss = 0;

  /* Try to detect IEEE 754 floating point.  The C99 standard specifies
     the preprocessor macro __STDC_IEC_559__ for this purpose, but GCC
     does not implement it.  The GNU C library defines it in the non 
     standard header file <features.h>, which is pulled in by <stdio.h>
     above.  On other GCC platforms, we have to try to detect it at
     runtime instead. */
#ifndef __STDC_IEC_559__
  /* Tests unity.  I'm sure this is not completely general, but it
     will spot VAX, probably the most likely of the alternative
     formats we could encounter. */
  if(*((uint64_t *) &d) == UINT64_C(0x3ff0000000000000))
    s->haveieee = 1;
  else
    s->haveieee = 0;
#endif
}

static void rng_gen (struct rng_state *s) {
  uint64_t *pq, aq, bq, ma, mb;
  int i, j, iq, jq;

  pq = s->a;

  aq = pq[RNG_RQ];
  bq = pq[RNG_RQ+1];

  for(i = 0; i < RNG_RR; i++) {
    j = (i + 117) % RNG_RR;

    iq = i+i;
    jq = j+j;

    ma = (bq << 32) | (bq >> 32);
    mb = (aq << 32) | (aq >> 32);

    aq = (pq[iq] << 19) ^ pq[jq] ^ ma;
    pq[iq] = (aq >> 12) ^ pq[iq] ^ (aq & MASKL);

    bq = (pq[iq+1] << 19) ^ pq[jq+1] ^ mb;
    pq[iq+1] = (bq >> 12) ^ pq[iq+1] ^ (bq & MASKH);
  }

  pq[RNG_RQ] = aq;
  pq[RNG_RQ+1] = bq;
}

void rng_fetch_uniform (struct rng_state *s, double *a, int n) {
  for(;;) {
    /* The generator produces IEEE 754 format numbers in [1,2)
       directly.  We can use them as-is (after subtracting unity)
       on an IEEE machine.  On non-IEEE machines, we convert from
       integer to floating point and divide by 2^52 instead.  The
       runtime haveieee flag is detected in rng_init() above to
       work around the problem of platforms that do use IEEE
       format but do not define the preprocessor macro. */
#ifdef __STDC_IEC_559__
    for(; n > 0 && s->r > 0; a++, n--, s->p++, s->r--)
      *a = *((double *) s->p) - 1.0;
#else
    if(s->haveieee)
      for(; n > 0 && s->r > 0; a++, n--, s->p++, s->r--)
        *a = *((double *) s->p) - 1.0;
    else
      /* Inefficient, given all the effort we went to above to
         make numbers already in IEEE format, but I don't expect
         to encounter many non-IEEE machines these days.  In
         other words, this version is not intended to be used. :) */
      for(; n > 0 && s->r > 0; a++, n--, s->p++, s->r--)
        *a = ldexp(*(s->p) & UINT64_C(0x000fffffffffffff), -52);
#endif

    if(n == 0)
      return;  /* we are done */

    rng_gen(s);
    s->p = &(s->a[0]);
    s->r = RNG_RQ;
  }
}

double rng_fetch_one_uniform (struct rng_state *s) {
  double rv;

  rng_fetch_uniform(s, &rv, 1);

  return(rv);
}

void rng_fetch_gauss (struct rng_state *s, double *a, int n) {
  double *p, v[2], rsq, fac;
  int i, x, g;

  /* Check if we already have one. */
  if(s->havegauss) {
    *a++ = s->gauss;
    n--;

    s->havegauss = 0;
  }

  /* This loop generates pairs of Gaussian deviates using the supplied
     array as a buffer.  The implementation takes a simple approach,
     generating only as many uniform deviates as it knows it will need.
     There is some overhead due to the repeated calls to the uniform
     deviate generator and the extra -1.0 done inside it. */
  x = n % 2;
  n -= x;
  
  while(n > 0) {
    rng_fetch_uniform(s, a, n);

    for(i = 0, g = 0, p = a; i < n; i += 2, p += 2) {
      v[0] = p[0] + p[0] - 1.0;
      v[1] = p[1] + p[1] - 1.0;
      rsq = v[0]*v[0] + v[1]*v[1];

      if(rsq > 0.0 && rsq < 1.0) {
        fac = sqrt(-2.0 * log(rsq) / rsq);
        *a++ = v[0] * fac;
        *a++ = v[1] * fac;
        g += 2;
      }
    }

    n -= g;
  }

  /* This loop produces the last Gaussian deviate in the odd n
     case.  Since we generate two doing this, the other is
     stashed for later use. */
  while(x > 0) {
    rng_fetch_uniform(s, v, 2);

    v[0] += v[0] - 1.0;
    v[1] += v[1] - 1.0;
    rsq = v[0]*v[0] + v[1]*v[1];

    if(rsq > 0.0 && rsq < 1.0) {
      fac = sqrt(-2.0 * log(rsq) / rsq);
      *a++ = v[0] * fac;
      x--;

      s->gauss = v[1] * fac;
      s->havegauss = 1;
    }
  }
}

double rng_fetch_one_gauss (struct rng_state *s) {
  double rv;

  rng_fetch_gauss(s, &rv, 1);

  return(rv);
}

/* Multivariate Gaussian deviates.  User-passed workspace can be used to
   store Cholesky decomposition of covariance matrix if multiple sets of
   deviates need to be generated for the same covariance. */

int rng_fetch_mvgauss (struct rng_state *s,
                       double *mean, double *cov,
                       double *work,
                       double *ans, int n) {
  int i, j, rv = 0;
  double z[n], x;

  if(cov) {
    /* Copy into workspace */
    for(i = 0; i < n; i++)
      for(j = 0; j <= i; j++)
        work[i*n+j] = cov[i*n+j];
    
    /* Cholesky decomposition of the covariance matrix */
    rv = cholesky(work, n);
  }
  /* otherwise, assume Cholesky is already done */

  if(ans) {
    /* ans = mean + work * Z where Z is a vector of Gaussian deviates */
    
    /* Generate Gaussian deviates */
    rng_fetch_gauss(s, z, n);
    
    for(i = 0; i < n; i++) {
      x = 0;

      for(j = 0; j <= i; j++)
        x += work[i*n+j] * z[j];

      ans[i] = mean[i] + x;
    }
  }

  return(rv);
}

