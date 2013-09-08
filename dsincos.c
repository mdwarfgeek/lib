#include <stdint.h>
#include <math.h>

/* For FORTRAN, compile with -DFORTRAN and call as if this was defined:
   SUBROUTINE DSINCOS (a, s, c)
   DOUBLE PRECISION a, s, c
 */

#ifdef FORTRAN
void dsincos_ (double *a, double *s, double *c) {
#else
void dsincos (double a, double *s, double *c) {
#endif

  /* Check machine architecture */
#if defined(__GNUC__) && (defined(__i386) || defined(__amd64))
  /* x87 assembler routine */
  __asm__ (
      "fldl    %2"        "\n\t"  /* load angle */
      "fsincos"           "\n\t"  /* compute sin, cos */
      "fnstsw  %%ax"      "\n\t"  /* get status word */
      "btw     $10, %%ax" "\n\t"  /* check bit 10 = C2, computed? */
      "jnc     1f"        "\n\t"  /* keep going if need to reduce arguments */
      "fldpi"             "\n\t"  /* pi */
      "fadd    %%st(0)"   "\n\t"  /* 2*pi */
      "fxch    %%st(1)"   "\n\t"  /* put args in correct order */
"2: " "fprem1"            "\n\t"  /* angle modulo 2*pi */
      "fnstsw  %%ax"      "\n\t"
      "btw     $10, %%ax" "\n\t"
      "jc      2b"        "\n\t"  /* loop until reduction is complete */
      "fstp    %%st(1)"   "\n\t"  /* clean up stack */
      "fsincos"           "\n\t"  /* compute final sin, cos */
"1: " "fstpl   %1"        "\n\t"  /* store outputs */
      "fstpl   %0"
      : "=m" (*s), "=m" (*c)
#ifdef FORTRAN
      : "m" (*a)
#else
      : "m" (a)
#endif
      : "ax", "st", "st(1)", "cc");

#else
  /* Generic C implementation using library functions */
  *s = sin(*a);
  *c = cos(*a);
#endif
}

