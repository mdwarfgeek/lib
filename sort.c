#include <stdlib.h>
#include <string.h>

#include "lfa.h"

#define PART(il, ir, DATATYPE) {                                \
  /* Select median of first, middle and last for pivot */       \
  vl = V(il, DATATYPE);                                         \
  vr = V(ir, DATATYPE);                                         \
                                                                \
  ipiv = il + (ir-il)/2;                                        \
  vpiv = V(ipiv, DATATYPE);                                     \
                                                                \
  if(vl < vr) {                                                 \
    if(vr < vpiv) {                                             \
      ipiv = ir;                                                \
      vpiv = vr;                                                \
    }                                                           \
    else {                                                      \
      if(vpiv < vl) {                                           \
        ipiv = il;                                              \
        vpiv = vl;                                              \
      }                                                         \
      /* else ipiv, already set */                              \
                                                                \
      /* Move pivot to end */                                   \
      X(ipiv, ir);                                              \
    }                                                           \
  }                                                             \
  else {                                                        \
    if(vpiv < vr) {                                             \
      ipiv = ir;                                                \
      vpiv = vr;                                                \
    }                                                           \
    else {                                                      \
      if(vl < vpiv) {                                           \
        ipiv = il;                                              \
        vpiv = vl;                                              \
      }                                                         \
      /* else ipiv, already set */                              \
                                                                \
      /* Move pivot to end */                                   \
      X(ipiv, ir);                                              \
    }                                                           \
  }                                                             \
                                                                \
  nsr = 1;                                                      \
                                                                \
  /* Move everything less than pivot to start, and everything   \
     equal to pivot to end. */                                  \
  isl = il;                                                     \
  for(i = il; i+nsr <= ir; ) {                                  \
    vcur = V(i, DATATYPE);                                      \
                                                                \
    if(vcur < vpiv) {                                           \
      X(isl, i);                                                \
      isl++;                                                    \
      i++;                                                      \
    }                                                           \
    else if(vcur == vpiv) {                                     \
      X(ir-nsr, i);                                             \
      nsr++;                                                    \
    }                                                           \
    else                                                        \
      i++;                                                      \
  }                                                             \
                                                                \
  /* Move pivot and any equal elements to final location */     \
  for(i = 0; i < nsr; i++) {                                    \
    X(isl+i, ir-i);                                             \
  }                                                             \
                                                                \
  /* Elements isl..isr (inclusive) equal to pivot */            \
  isr = isl + nsr - 1;                                          \
}

/* QuickSelect, using the pivot and partition macro above */

#define MAKE_quickselect(DATATYPE) {                                    \
  size_t il, ir, isl, isr, ipiv, nsr, i;                                \
  DATATYPE vl, vr, vpiv, vcur;                                          \
  TEMPDECL(DATATYPE);                                                   \
                                                                        \
  if(n <= 1)                                                            \
    return R(0);                                                        \
                                                                        \
  il = 0;                                                               \
  ir = n-1;                                                             \
                                                                        \
  for(;;) {                                                             \
    PART(il, ir, DATATYPE);                                             \
                                                                        \
    if(k < isl)                                                         \
      ir = isl-1;                                                       \
    else if(k > isr)                                                    \
      il = isr+1;                                                       \
    else  /* it's one of the elements equal to the pivot */             \
      return R(isl);  /* doesn't really matter which one */             \
  }                                                                     \
}

/* QuickSort */

#define MAKE_quicksort(FUNC, DATATYPE) {                        \
  size_t il, ir, isl, isr, nl, nr, ipiv, nsr, i;                \
  DATATYPE vl, vr, vpiv, vcur;                                  \
  TEMPDECL(DATATYPE);                                           \
                                                                \
  il = 0;                                                       \
  ir = n-1;                                                     \
                                                                \
  while(il < ir) {                                              \
    PART(il, ir, DATATYPE);                                     \
                                                                \
    if(il < isl) {  /* work to do on the left? */               \
      /* work to do on the right? */                            \
      if(ir > isr) {                                            \
        nl = isl-il;                                            \
        nr = ir-isr;                                            \
                                                                \
        /* Recurse into the smaller side first */               \
        if(nl > nr) {                                           \
          RECURSE(FUNC, B(isr+1), nr);                          \
          ir = isl-1;                                           \
        }                                                       \
        else {                                                  \
          RECURSE(FUNC, B(il), nl);                             \
          il = isr+1;                                           \
        }                                                       \
      }                                                         \
      else                                                      \
        ir = isl-1;                                             \
    }                                                           \
    else if(ir > isr) {  /* work to do on the right? */         \
      /* Do the right */                                        \
      il = isr+1;                                               \
    }                                                           \
    else  /* we are done */                                     \
      return;                                                   \
  }                                                             \
}

/* "Simple" versions to sort a straightforward array of the given datatype */
#define B(i) (&(list[i]))
#define V(i, t) (list[i])
#define R(i) (list[i])
#define X(i, j)                                 \
  tmp = list[i];                                \
  list[i] = list[j];                            \
  list[j] = tmp;
#define TEMPDECL(t) t tmp
#define RECURSE(f, p, n) f(p, n)

double dquickselect (double *list, size_t k, size_t n)
  MAKE_quickselect(double)

float fquickselect (float *list, size_t k, size_t n)
  MAKE_quickselect(float)

int iquickselect (int *list, size_t k, size_t n)
  MAKE_quickselect(int)

void dquicksort (double *list, size_t n)
  MAKE_quicksort(dquicksort, double);

void fquicksort (float *list, size_t n)
  MAKE_quicksort(fquicksort, float);

void iquicksort (int *list, size_t n)
  MAKE_quicksort(iquicksort, int);

#undef B
#undef V
#undef R
#undef X
#undef TEMPDECL
#undef RECURSE

/* "Generic" versions to sort arrays of larger items using keys contained
   at arbitrary (provided they are aligned) byte offsets within. */
#define B(i) (((char *) list) + (i)*s)
#define V(i, t) *((t *) (B(i)+o))
#define R(i) B(i)
#define X(i, j)                                 \
  memcpy(tmp, B(i), s);                         \
  memcpy(B(i), B(j), s);                        \
  memcpy(B(j), tmp, s)
#define TEMPDECL(t) 
#define RECURSE(f, p, n) f(p, tmp, n, s, o)

void *dquickselect_gen (void *list, void *tmp,
                        size_t k, size_t n, size_t s, size_t o)
  MAKE_quickselect(double)

void *fquickselect_gen (void *list, void *tmp,
                          size_t k, size_t n, size_t s, size_t o)
  MAKE_quickselect(float)

void *iquickselect_gen (void *list, void *tmp,
                        size_t k, size_t n, size_t s, size_t o)
  MAKE_quickselect(int)

void dquicksort_gen (void *list, void *tmp,
                     size_t n, size_t s, size_t o)
  MAKE_quicksort(dquicksort_gen, double);

void fquicksort_gen (void *list, void *tmp,
                     size_t n, size_t s, size_t o)
  MAKE_quicksort(fquicksort_gen, float);

void iquicksort_gen (void *list, void *tmp,
                     size_t n, size_t s, size_t o)
  MAKE_quicksort(iquicksort_gen, int);

/* My usual "medsig" routine using median and MAD scaled to Gaussian
   RMS equivalent.  Lives here because I couldn't think of a better
   place for it. */

#define MAKE_medsig(DATATYPE, TQS, TABS) {                      \
  size_t i, m;                                                  \
  DATATYPE median;                                              \
                                                                \
  if(n < 1) {                                                   \
    /* Special case - can't do anything */                      \
    if(median_r)                                                \
      *median_r = 0;                                            \
                                                                \
    if(sigma_r)                                                 \
      *sigma_r = 0;                                             \
                                                                \
    return;                                                     \
  }                                                             \
                                                                \
  m = n/2;                                                      \
                                                                \
  median = TQS(list, m, n);                                     \
                                                                \
  if(median_r)                                                  \
    *median_r = median;                                         \
                                                                \
  if(sigma_r) {                                                 \
    for(i = 0; i < n; i++)                                      \
      list[i] = TABS(list[i] - median);                         \
                                                                \
    *sigma_r = 1.482602218505601 * TQS(list, m, n);             \
  }                                                             \
}

void dmedsig (double *list, size_t n, double *median_r, double *sigma_r)
  MAKE_medsig(double, dquickselect, fabs)

void fmedsig (float *list, size_t n, float *median_r, float *sigma_r)
  MAKE_medsig(float, fquickselect, fabsf)
