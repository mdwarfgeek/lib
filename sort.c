#include <stdlib.h>
#include <string.h>

#define B(i) (((char *) list) + (i)*s)
#define V(i, t) *((t *) (B(i)+o))
#define X(i, j)                                 \
  memcpy(tmp, B(i), s);                         \
  memcpy(B(i), B(j), s);                        \
  memcpy(B(j), tmp, s)
#define PART(il, ir, DATATYPE) {                                \
  /* Use middle element as pivot */                             \
  ipiv = il + (ir-il)/2;                                        \
  vpiv = V(ipiv, DATATYPE);                                     \
                                                                \
  /* Move pivot to end */                                       \
  X(ipiv, ir);                                                  \
  isr = ir-1;                                                   \
                                                                \
  /* Move everything less than pivot to start, and everything   \
     equal to pivot to end. */                                  \
  isl = il;                                                     \
  for(i = il; i <= isr; i++) {                                  \
    vcur = V(i, DATATYPE);                                      \
                                                                \
    if(vcur < vpiv) {                                           \
      X(isl, i);                                                \
      isl++;                                                    \
    }                                                           \
    else if(vcur == vpiv) {                                     \
      X(isr, i);                                                \
      isr--;                                                    \
    }                                                           \
  }                                                             \
                                                                \
  /* Move pivot and any equal elements to final location */     \
  nsr = ir-isr;                                                 \
                                                                \
  for(i = 0; i < nsr; i++) {                                    \
    X(isl+i, ir-i);                                             \
  }                                                             \
                                                                \
  /* Elements isl..isr (inclusive) equal to pivot */            \
  isr = isl + nsr - 1;                                          \
}

/* QuickSelect, using the pivot and partition macro above */

#define MAKE_quickselect(FUNC, DATATYPE)                                \
void *FUNC (void *list, size_t k, size_t n, size_t s, size_t o) {       \
  size_t il, ir, isl, isr, ipiv, nsr, i;                                \
  DATATYPE vpiv, vcur;                                                  \
  char tmp[s];                                                          \
                                                                        \
  if(n <= 1)                                                            \
    return list;                                                        \
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
      return B(isl);  /* doesn't really matter which one */             \
  }                                                                     \
}

MAKE_quickselect(dquickselect, double)
MAKE_quickselect(fquickselect, float)
MAKE_quickselect(iquickselect, int)

/* QuickSort, with internal recursion to avoid function call overhead */

#define MAKE_quicksort(FUNC, DATATYPE)                          \
void FUNC (void *list, size_t n, size_t s, size_t o) {          \
  size_t il, ir, isl, isr, ipiv, nsr, i;                        \
  DATATYPE vpiv, vcur;                                          \
  char tmp[s];                                                  \
                                                                \
  struct {                                                      \
    size_t il;                                                  \
    size_t ir;                                                  \
  } stack[n];                                                   \
  int nstack = 0;                                               \
                                                                \
  /* Trap for single or zero element lists */                   \
  if(n <= 1)                                                    \
    return;                                                     \
                                                                \
  il = 0;                                                       \
  ir = n-1;                                                     \
                                                                \
  for(;;) {                                                     \
    PART(il, ir, DATATYPE);                                     \
                                                                \
    isl--;                                                      \
    isr++;                                                      \
                                                                \
    if(il < isl) {  /* work to do on the left? */               \
      /* Stack any work on the right */                         \
      if(ir > isr) {                                            \
        stack[nstack].il = isr;                                 \
        stack[nstack].ir = ir;                                  \
        nstack++;                                               \
      }                                                         \
                                                                \
      /* Do the left */                                         \
      ir = isl;                                                 \
    }                                                           \
    else if(ir > isr) {  /* work to do on the right? */         \
      /* Do the right */                                        \
      il = isr;                                                 \
    }                                                           \
    else if(nstack > 0) {  /* work to do on the stack? */       \
      nstack--;                                                 \
      il = stack[nstack].il;                                    \
      ir = stack[nstack].ir;                                    \
    }                                                           \
    else  /* we are done */                                     \
      return;                                                   \
  }                                                             \
}

MAKE_quicksort(dquicksort, double);
MAKE_quicksort(fquicksort, float);
MAKE_quicksort(iquicksort, int);
