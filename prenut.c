#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lfa.h"
#include "cvtunit.h"
#include "util.h"

static double ang_sum_poly (double t, const double *poly, int npoly) {
  double ang = 0;
  int ipoly;

  for(ipoly = npoly-1; ipoly >= 0; ipoly--)
    ang = ang * t + poly[ipoly];

  return(fmod(ang, AS_PER_REV) * AS_TO_RAD);
}

void pfb06ang (double t, double ang[NPNANG]) {
  /* Precession and frame bias, IAU 2006, Fukushima-Williams angles.
     Coefficients from Hilton et al. 2006 */
  static const double prec_coef[NPNANG][6] = {
    {    -0.052928,   10.556378,  0.4932044, -0.00031238, -2.788e-6,   2.60e-8  }, /* gam */
    { 84381.412819,  -46.811016,  0.0511268,  0.00053289, -4.40e-7,   -1.76e-8  }, /* phi */
    {    -0.041775, 5038.481484,  1.5584175, -0.00018522, -2.6452e-5, -1.48e-8  }, /* psi */
    { 84381.406000,  -46.836769, -0.0001831,  0.00200340, -5.76e-7,   -4.34e-8  }  /* epsa */
  };

  int np, iang;

  np = sizeof(prec_coef)/sizeof(prec_coef[0]);

  for(iang = 0; iang < NPNANG; iang++)
    ang[iang] = ang_sum_poly(t, prec_coef[iang], np);
}

void nut00b (double t,
	     double *dpsi_r, double *deps_r,
	     double *sxy_r, double *eo_r) {
  /* Delaunay arguments ("1992 values" of constants) and moon mean
     longitude referred to mean ecliptic and equinox of date.
     Abridged from Simon et al. 1994. */
  static const double del_coef[5][2] = {
    {  485868.249036, 1717915923.2178 },  /* l */
    { 1287104.79305,   129596581.0481 },  /* l' */
    {  335779.526232, 1739527262.8478 },  /* F */
    { 1072260.70369,  1602961601.2090 },  /* D */
    {  450160.398036,   -6962890.5431 }   /* omega */
  };

  /* Luni-solar nutation, IAU 2000B.
     Coefficients from IERS 2003 conventions, Chapter 5, NU2000B.f
     For some reason, the series in there omits the last two terms in
     the McCarthy & Luzum (2003) paper (and IAU2000B.f, which is the
     same) and adds one new term.  Comparisons with the SOFA IAU 2000A
     routine indicate that this treatment is indeed closer, so it has
     been adopted for the present implementation.  The difference may
     be related to use of the rigorous method with the 2000 precession,
     whereas the McCarthy & Luzum implementation gave corrections to the
     classical IAU 1976 precession. */
  static const struct {
    int m[5];
    double a, ap, b, bp, app, bpp;  /* units 0.1 uas */
  } nut_coef[] = {
    { {  0,  0,  0,  0,  1 }, -172064161.0, -174666.0, 92052331.0,  9086.0,  33386.0, 15377.0 },
    { {  0,  0,  2, -2,  2 },  -13170906.0,   -1675.0,  5730336.0, -3015.0, -13696.0, -4587.0 },
    { {  0,  0,  2,  0,  2 },   -2276413.0,    -234.0,   978459.0,  -485.0,   2796.0,  1374.0 },
    { {  0,  0,  0,  0,  2 },    2074554.0,     207.0,  -897492.0,   470.0,   -698.0,  -291.0 },
    { {  0,  1,  0,  0,  0 },    1475877.0,   -3633.0,    73871.0,  -184.0,  11817.0, -1924.0 },
    { {  0,  1,  2, -2,  2 },    -516821.0,    1226.0,   224386.0,  -677.0,   -524.0,  -174.0 },
    { {  1,  0,  0,  0,  0 },     711159.0,      73.0,    -6750.0,     0.0,   -872.0,   358.0 },
    { {  0,  0,  2,  0,  1 },    -387298.0,    -367.0,   200728.0,    18.0,    380.0,   318.0 },
    { {  1,  0,  2,  0,  2 },    -301461.0,     -36.0,   129025.0,   -63.0,    816.0,   367.0 },
    { {  0, -1,  2, -2,  2 },     215829.0,    -494.0,   -95929.0,   299.0,    111.0,   132.0 },
    { {  0,  0,  2, -2,  1 },     128227.0,     137.0,   -68982.0,    -9.0,    181.0,    39.0 },
    { { -1,  0,  2,  0,  2 },     123457.0,      11.0,   -53311.0,    32.0,     19.0,    -4.0 },
    { { -1,  0,  0,  2,  0 },     156994.0,      10.0,    -1235.0,     0.0,   -168.0,    82.0 },
    { {  1,  0,  0,  0,  1 },      63110.0,      63.0,   -33228.0,     0.0,     27.0,    -9.0 },
    { { -1,  0,  0,  0,  1 },     -57976.0,     -63.0,    31429.0,     0.0,   -189.0,   -75.0 },
    { { -1,  0,  2,  2,  2 },     -59641.0,     -11.0,    25543.0,   -11.0,    149.0,    66.0 },
    { {  1,  0,  2,  0,  1 },     -51613.0,     -42.0,    26366.0,     0.0,    129.0,    78.0 },
    { { -2,  0,  2,  0,  1 },      45893.0,      50.0,   -24236.0,   -10.0,     31.0,    20.0 },
    { {  0,  0,  0,  2,  0 },      63384.0,      11.0,    -1220.0,     0.0,   -150.0,    29.0 },
    { {  0,  0,  2,  2,  2 },     -38571.0,      -1.0,    16452.0,   -11.0,    158.0,    68.0 },
    { {  0, -2,  2, -2,  2 },      32481.0,       0.0,   -13870.0,     0.0,      0.0,     0.0 },
    { { -2,  0,  0,  2,  0 },     -47722.0,       0.0,      477.0,     0.0,    -18.0,   -25.0 },
    { {  2,  0,  2,  0,  2 },     -31046.0,      -1.0,    13238.0,   -11.0,    131.0,    59.0 },
    { {  1,  0,  2, -2,  2 },      28593.0,       0.0,   -12338.0,    10.0,     -1.0,    -3.0 },
    { { -1,  0,  2,  0,  1 },      20441.0,      21.0,   -10758.0,     0.0,     10.0,    -3.0 },
    { {  2,  0,  0,  0,  0 },      29243.0,       0.0,     -609.0,     0.0,    -74.0,    13.0 },
    { {  0,  0,  2,  0,  0 },      25887.0,       0.0,     -550.0,     0.0,    -66.0,    11.0 },
    { {  0,  1,  0,  0,  1 },     -14053.0,     -25.0,     8551.0,    -2.0,     79.0,   -45.0 },
    { { -1,  0,  0,  2,  1 },      15164.0,      10.0,    -8001.0,     0.0,     11.0,    -1.0 },
    { {  0,  2,  2, -2,  2 },     -15794.0,      72.0,     6850.0,   -42.0,    -16.0,    -5.0 },
    { {  0,  0, -2,  2,  0 },      21783.0,       0.0,     -167.0,     0.0,     13.0,    13.0 },
    { {  1,  0,  0, -2,  1 },     -12873.0,     -10.0,     6953.0,     0.0,    -37.0,   -14.0 },
    { {  0, -1,  0,  0,  1 },     -12654.0,      11.0,     6415.0,     0.0,     63.0,    26.0 },
    { { -1,  0,  2,  2,  1 },     -10204.0,       0.0,     5222.0,     0.0,     25.0,    15.0 },
    { {  0,  2,  0,  0,  0 },      16707.0,     -85.0,      168.0,    -1.0,    -10.0,    10.0 },
    { {  1,  0,  2,  2,  2 },      -7691.0,       0.0,     3268.0,     0.0,     44.0,    19.0 },
    { { -2,  0,  2,  0,  0 },     -11024.0,       0.0,      104.0,     0.0,    -14.0,     2.0 },
    { {  0,  1,  2,  0,  2 },       7566.0,     -21.0,    -3250.0,     0.0,    -11.0,    -5.0 },
    { {  0,  0,  2,  2,  1 },      -6637.0,     -11.0,     3353.0,     0.0,     25.0,    14.0 },
    { {  0, -1,  2,  0,  2 },      -7141.0,      21.0,     3070.0,     0.0,      8.0,     4.0 },
    { {  0,  0,  0,  2,  1 },      -6302.0,     -11.0,     3272.0,     0.0,      2.0,     4.0 },
    { {  1,  0,  2, -2,  1 },       5800.0,      10.0,    -3045.0,     0.0,      2.0,    -1.0 },
    { {  2,  0,  2, -2,  2 },       6443.0,       0.0,    -2768.0,     0.0,     -7.0,    -4.0 },
    { { -2,  0,  0,  2,  1 },      -5774.0,     -11.0,     3041.0,     0.0,    -15.0,    -5.0 },
    { {  2,  0,  2,  0,  1 },      -5350.0,       0.0,     2695.0,     0.0,     21.0,    12.0 },
    { {  0, -1,  2, -2,  1 },      -4752.0,     -11.0,     2719.0,     0.0,     -3.0,    -3.0 },
    { {  0,  0,  0, -2,  1 },      -4940.0,     -11.0,     2720.0,     0.0,    -21.0,    -9.0 },
    { { -1, -1,  0,  2,  0 },       7350.0,       0.0,      -51.0,     0.0,     -8.0,     4.0 },
    { {  2,  0,  0, -2,  1 },       4065.0,       0.0,    -2206.0,     0.0,      6.0,     1.0 },
    { {  1,  0,  0,  2,  0 },       6579.0,       0.0,     -199.0,     0.0,    -24.0,     2.0 },
    { {  0,  1,  2, -2,  1 },       3579.0,       0.0,    -1900.0,     0.0,      5.0,     1.0 },
    { {  1, -1,  0,  0,  0 },       4725.0,       0.0,      -41.0,     0.0,     -6.0,     3.0 },
    { { -2,  0,  2,  0,  2 },      -3075.0,       0.0,     1313.0,     0.0,     -2.0,    -1.0 },
    { {  3,  0,  2,  0,  2 },      -2904.0,       0.0,     1233.0,     0.0,     15.0,     7.0 },
    { {  0, -1,  0,  2,  0 },       4348.0,       0.0,      -81.0,     0.0,    -10.0,     2.0 },
    { {  1, -1,  2,  0,  2 },      -2878.0,       0.0,     1232.0,     0.0,      8.0,     4.0 },
    { {  0,  0,  0,  1,  0 },      -4230.0,       0.0,      -20.0,     0.0,      5.0,    -2.0 },
    { { -1, -1,  2,  2,  2 },      -2819.0,       0.0,     1207.0,     0.0,      7.0,     3.0 },
    { { -1,  0,  2,  0,  0 },      -4056.0,       0.0,       40.0,     0.0,      5.0,    -2.0 },
    { {  0, -1,  2,  2,  2 },      -2647.0,       0.0,     1129.0,     0.0,     11.0,     5.0 },
    { { -2,  0,  0,  0,  1 },      -2294.0,       0.0,     1266.0,     0.0,    -10.0,    -4.0 },
    { {  1,  1,  2,  0,  2 },       2481.0,       0.0,    -1062.0,     0.0,     -7.0,    -3.0 },
    { {  2,  0,  0,  0,  1 },       2179.0,       0.0,    -1129.0,     0.0,     -2.0,    -2.0 },
    { { -1,  1,  0,  1,  0 },       3276.0,       0.0,       -9.0,     0.0,      1.0,     0.0 },
    { {  1,  1,  0,  0,  0 },      -3389.0,       0.0,       35.0,     0.0,      5.0,    -2.0 },
    { {  1,  0,  2,  0,  0 },       3339.0,       0.0,     -107.0,     0.0,    -13.0,     1.0 },
    { { -1,  0,  2, -2,  1 },      -1987.0,       0.0,     1073.0,     0.0,     -6.0,    -2.0 },
    { {  1,  0,  0,  0,  2 },      -1981.0,       0.0,      854.0,     0.0,      0.0,     0.0 },
    { { -1,  0,  0,  1,  0 },       4026.0,       0.0,     -553.0,     0.0,   -353.0,  -139.0 },
    { {  0,  0,  2,  1,  2 },       1660.0,       0.0,     -710.0,     0.0,     -5.0,    -2.0 },
    { { -1,  0,  2,  4,  2 },      -1521.0,       0.0,      647.0,     0.0,      9.0,     4.0 },
    { { -1,  1,  0,  1,  1 },       1314.0,       0.0,     -700.0,     0.0,      0.0,     0.0 },
    { {  0, -2,  2, -2,  1 },      -1283.0,       0.0,      672.0,     0.0,      0.0,     0.0 },
    { {  1,  0,  2,  2,  1 },      -1331.0,       0.0,      663.0,     0.0,      8.0,     4.0 },
    { { -2,  0,  2,  2,  2 },       1383.0,       0.0,     -594.0,     0.0,     -2.0,    -2.0 },
    { { -1,  0,  0,  0,  2 },       1405.0,       0.0,     -610.0,     0.0,      4.0,     2.0 },
    { {  1,  1,  2, -2,  2 },       1290.0,       0.0,     -556.0,     0.0,      0.0,     0.0 }
  };

  /* Secular component of CIO locator */
  static const double sxy_coef[6] =
    { 0.0000940, 0.00380865, -0.00012268, -0.07257411,  2.798e-5, 1.562e-5 };

  double arg, sinarg, cosarg;
  double dpsi, deps;
  double nut_arg[5];
  int iarg, narg;
  int ils, nls;

  double scs, tsq;

  /* Nutation, IAU 2000B */
  dpsi = 0.0;
  deps = 0.0;
  
  /* Compute arguments */
  narg = sizeof(del_coef) / sizeof(del_coef[0]);

  for(iarg = 0; iarg < narg; iarg++)
    nut_arg[iarg] = fmod(del_coef[iarg][0] +
			 del_coef[iarg][1] * t, AS_PER_REV) * AS_TO_RAD;

  /* Sum luni-solar series, smallest terms first */
  nls = sizeof(nut_coef) / sizeof(nut_coef[0]);
  
  //printf("%d nutation terms\n", nls);
  
  for(ils = nls-1; ils >= 0; ils--) {
    /* Argument */
    arg = 0;
    
    for(iarg = 0; iarg < 5; iarg++)
      arg += nut_coef[ils].m[iarg] * nut_arg[iarg];
    
    dsincos(arg, &sinarg, &cosarg);

    /* Accumulate sums */
    dpsi += (nut_coef[ils].a + nut_coef[ils].ap * t) * sinarg + nut_coef[ils].app * cosarg;
    deps += (nut_coef[ils].b + nut_coef[ils].bp * t) * cosarg + nut_coef[ils].bpp * sinarg;
  }
  
  /* Offsets to compensate for planetary terms.  The paper gave corrections
     to the IAU 1976 precession, but we don't need this.  These values are
     appropriate for the modern precession formulae, and came from the NOVAS
     implementation. */
  dpsi -= 1.35e3;  /* units 0.1 uas */
  deps += 3.88e3;  /* units 0.1 uas */
  
  /* Change units to rad */
  dpsi *= 1.0e-7 * AS_TO_RAD;
  deps *= 1.0e-7 * AS_TO_RAD;

  /* Return */
  (*dpsi_r) = dpsi;
  (*deps_r) = deps;

  if(sxy_r) {
    /* Compute secular component of CIO locator */
    scs = ang_sum_poly(t, sxy_coef, sizeof(sxy_coef)/sizeof(sxy_coef[0]));

    /* Compute periodic components of CIO locator.  This series
       is based on the one provided in IERS TN 36, but has been
       further truncated, removing most of the terms smaller than
       20 uas in 100 yrs.  Two of the smaller quadratic terms have
       been retained to keep a lid on the behaviour at epochs
       very far from 2000.0, in case the routine is inadvertently
       used there.  The results agree with the SOFA routine
       to better than 100 uas from about 1900 - 2100, so the
       accuracy of the final matrix should be limited by the
       nutation series. */
    tsq = t*t;
    
    (*sxy_r) = (+  9.84*tsq * sin(2*nut_arg[2]+2*nut_arg[4])
		+ 56.91*tsq * sin(2*nut_arg[2]-2*nut_arg[3]+2*nut_arg[4])
		- (    -8.85*tsq +   63.53) * sin(2*nut_arg[4])
		+ (   743.52*tsq - 2640.73) * sin(nut_arg[4])
		) * 1.0e-6 * AS_TO_RAD + scs;
  }

  if(eo_r) {
    /* Equation of the origins, to be implemented */
  }
}

void makepnm (double ang[NPNANG],
	      double m[3][3]) {

  m_identity(m);
  euler_rotate(m, 3,  ang[PNANG_GAM]);
  euler_rotate(m, 1,  ang[PNANG_PHI]);
  euler_rotate(m, 3, -ang[PNANG_PSI]);
  euler_rotate(m, 1, -ang[PNANG_EPSA]);
}

void makecim (double ang[NPNANG], double sxy,
	      double dxnut, double dynut,
	      double m[3][3]) {
  double sinang[NPNANG], cosang[NPNANG];
  int iang;

  double acip, bcip, xcip, ycip, xcipsq, ycipsq, xycip;
  double scio, a;

  /* Precompute sine and cosine of angles */
  for(iang = 0; iang < NPNANG; iang++)
    dsincos(ang[iang], sinang+iang, cosang+iang);
  
  /* Compute X and Y coordinates of CIP */
  acip = sinang[PNANG_EPSA]*sinang[PNANG_PSI];
  bcip = sinang[PNANG_EPSA]*cosang[PNANG_PSI]*cosang[PNANG_PHI]
    - cosang[PNANG_EPSA]*sinang[PNANG_PHI];
  
  xcip = acip*cosang[PNANG_GAM] - bcip*sinang[PNANG_GAM] + dxnut;
  ycip = acip*sinang[PNANG_GAM] + bcip*cosang[PNANG_GAM] + dynut;

  xycip = xcip*ycip;

  scio = sxy - 0.5*xycip;

  /* Approximation to "a" from IERS TN36 accurate to 1 microarcsec */
  xcipsq = xcip*xcip;
  ycipsq = ycip*ycip;
  
  a = 0.5 + 0.125*(xcipsq+ycipsq);
  
  /* GCRS-CIRS matrix, transpose of Eq. 5.10 of IERS TN36 */
  m[0][0] = 1.0 - a*xcipsq;
  m[0][1] = -a*xycip;
  m[0][2] = -xcip;
  m[1][0] = m[0][1];
  m[1][1] = 1.0 - a*ycipsq;
  m[1][2] = -ycip;
  m[2][0] = xcip;
  m[2][1] = ycip;
  m[2][2] = 1.0 - a*(xcipsq+ycipsq);

  euler_rotate(m, 3, -scio);
}
