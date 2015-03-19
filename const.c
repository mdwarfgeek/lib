/* const.c: this file contains various pre-computed matrices */

#include "lfa.h"

/* Precomputed GCRS to ecliptic of J2000, using IAU 2006 precession
   formulation in prenut.c */

double gcrs2ecl[3][3] = {
  {  9.99999999999994116e-1,  3.28970040774196531e-8, -1.02070447254843567e-7 },
  { -7.07836896097155613e-8,  9.17482129914958477e-1, -3.97776999444043045e-1 },
  {  8.05621397761318608e-8,  3.97776999444047985e-1,  9.17482129914955591e-1 }
};

/* Precomputed FK5 to ecliptic of J2000, using IAU 2006 precession
   formulation in prenut.c */

double fk52ecl[3][3] = {
  {  1.00000000000000000e+0,  0.00000000000000000e+0,  0.00000000000000000e+0 },
  {  0.00000000000000000e+0,  9.17482143065241784e-1, -3.97776969112606016e-1 },
  {  0.00000000000000000e+0,  3.97776969112606016e-1,  9.17482143065241784e-1 }
};

/* Precomputed equatorial ICRS to Galactic, following Hipparcos Vol. 1,
   Sect. 1.5.3, and using the constants given there:

   alpha(NGP) = 192.85948 deg
   delta(NGP) =  27.12825 deg
   l_omega    =  32.93192 deg
   
   As stated in the Hipparcos documentation, the distinction between
   ICRS and FK5 can be neglected given the uncertainties. */

double eq2gal[3][3] = {
  { -5.48755604162155797e-2, -8.73437090234885138e-1, -4.83835015548713165e-1 },
  {  4.94109427875583596e-1, -4.44829629960011241e-1,  7.46982244497218950e-1 },
  { -8.67666149019004740e-1, -1.98076373431201297e-1,  4.55983776175066968e-1 }
};
