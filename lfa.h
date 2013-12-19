/* lfa.h: positional astronomy subroutines.
          various source files, see below. */

/* TODO: fix up implementation in iers.c and other ephemeris type routines
         to ensure correct handling of running off the end.
         add local math header and some macros so we can use sincos()
         portably.  this should make it easier to use CORDIC efficiently
	 when porting to hardware later.  possibly also use hypot()
	 svr4/bsd/c99 routine where applicable, and investigate ways
	 to implement robust vector norms. */

#ifndef LFA_H
#define LFA_H

#include <stdio.h>  /* for FILE */
#include <stdint.h>

/* -- Physical, Mathematical and Astronomical constants -- */

/* NOTE: integer constants are expressed here without the .0, so be
   careful when using them if the result needs to be a floating point
   number. */

/* Fundamental (defining) constants, mostly IAU and IERS */

#define GMSUN  1.32712440041e20  /* m^3 / s^2, IAU 2009 system, TDB compatible */
#define AU     149597870700.0    /* m, IAU 2009 system */
#define LIGHT  2.99792458e8      /* m/s, definition */

#define GMEARTH 3.986004356e14     /* m^3 / s^2, IAU 2009 system, TDB compatible */
#define AEARTH  6378137.0          /* m, WGS84 */
#define FEARTH (1.0/298.257223563) /* WGS84 */

#define DAY    86400             /* s */
#define DTT    32.184            /* s, TT-TAI */

#define ZMJD   2400000.5         /* JD at MJD=0 */
#define JYR    365.25            /* days */
#define J2K    51544.5           /* J2000 as MJD */
#define JUNIX  40587             /* UNIX epoch as MJD */
#define JTCB   43144.0003725     /* TT, TCG, TCB = JTCB at 1977 Jan 1 0h TAI */

#define LB     1.550519768e-8    /* IAU 2006 Resol. B3 */
#define LG     6.969290134e-10   /* IAU 2000 Resol. B1.9 */
#define ZTDB  -6.55e-5           /* s, IAU 2006 Resol. B3 */

#define TWOPI      (2.0*M_PI)
#define AS_PER_REV 1296000

/* Earth rotation angle from Capitaine et al. 2000.
   ERA = TWOPI*(ERA2K + (1+ERADAY) * (UT1 - J2K)) */

#define ERA2K  0.7790572732640      /* ERA at J2000.0 UT1 */
#define ERADAY 0.00273781191135448  /* ERA increment relative to UT1 per Julian day */

#define EOMEGA  (TWOPI * (1.0 + ERADAY))  /* Earth angular velocity, rads/day */

/* Other astrophysical quantities */

#define RSUN   6.95508e8         /* m, from Brown & Christensen-Dalsgaard 1998
                                  * and as adopted by Cox 2000 (Allen's
                                  * Astrophysical Quantities, 4th Ed.) */

/* -- Other constants -- */

/* Buffer size for reading files */
#define BUFSIZE 1024

/* Number of refraction constants */
#define NREFCO  5

/* Precession angles */
#define PNANG_GAM  0
#define PNANG_PHI  1
#define PNANG_PSI  2
#define PNANG_EPSA 3
#define NPNANG     4

/* Constants specifying corrections to apply.  Many of the possible
   combinations don't make any sense, but there are no restrictions
   implemented in the present routine. */
#define TR_MOTION    0x01
#define TR_PLX       0x02
#define TR_DEFL      0x04
#define TR_ANNAB     0x08
#define TR_TOPO      0x10
#define TR_LAT       0x20
#define TR_REFRO     0x40

/* Combinations for some common operations */
#define TR_TO_AST     (TR_MOTION | TR_PLX)
#define TR_TO_GCRS    (TR_TO_AST | TR_DEFL | TR_ANNAB)
#define TR_TO_TOPO_HD (TR_TO_GCRS | TR_TOPO)
#define TR_TO_OBS_HD  (TR_TO_TOPO_HD | TR_REFRO)
#define TR_TO_TOPO_AZ (TR_TO_TOPO_HD | TR_LAT)
#define TR_TO_OBS_AZ  (TR_TO_TOPO_AZ | TR_REFRO)

/* -- Data structures -- */

struct jpleph_table {
  FILE *fp;
  int32_t denum;
  int32_t recsize;
  long taboff;

  double mjd_start;
  double mjd_end;
  double mjd_step;

  int32_t *ipt;
#define JPLEPH_MERCURY    0
#define JPLEPH_VENUS      1
#define JPLEPH_EMB        2
#define JPLEPH_MARS       3
#define JPLEPH_JUPITER    4
#define JPLEPH_SATURN     5
#define JPLEPH_URANUS     6
#define JPLEPH_NEPTUNE    7
#define JPLEPH_PLUTO      8
#define JPLEPH_MOON       9  /* geocenter to moon */
#define JPLEPH_SUN       10
#define JPLEPH_NUTATION  11
#define JPLEPH_LIBRATION 12

/* Time ephemerides, see http://timeephem.sourceforge.net/
   and Irwin & Fukushima (1999).  To use the integral,
   TDB-TT = ZTDB + (TEI(Teph) - TEI(Teph_0)) / (1.0 - LC).
   Using TT rather than Teph as the argument on the RHS
   introduces negligible error, as does dropping the
   TEI(Teph_0) term.  Similarly, assuming LC = 0 would
   cause a maximum error of 30 ps. */
#define TIMEEPH_TEI      14  /* time ephemeris integral */
#define TIMEEPH_TEV      15  /* time ephemeris vector */
  int32_t npt;

  double au;
  double emratio;
  double emfac;

  /* Constants for time ephemerides */
  double lc;

  double *buf;
  int32_t brec;
};

struct dtai_entry {
  int mjd;  /* leap seconds can only happen at midnight so always integer */
  float dtai;
  float mjdzero;
  float scale;
};

struct dtai_table {
  struct dtai_entry *table;
  int ilast;
  int ntab;
};

struct iers_entry {
  double mjd;

  /* Values */
  double xp;
  double yp;
  double dut1;
  double dxnut;
  double dynut;
};

struct iers_table {
  FILE *fp;
  int recsize;
  int lrec;     /* last record */
  double mjd_start;
  double mjd_step;

  struct dtai_table *dtai_tab;

  /* Records bracketing desired MJD for interpolation */
  struct iers_entry buf[2];
  double h;

  int bufrec;   /* record before target in buffer pos 1 */
  int filerec;  /* current file position */
};

struct observer {
  /* Static quantities */
  struct jpleph_table *jpltab;
  struct jpleph_table *tetab;
  struct iers_table *iertab;

  double latitude;
  double longitude;
  double height;

  double lmm[3][3];        /* Observer longitude matrix */
  double phm[3][3];        /* Observer latitude matrix */
  double itgop[3];         /* ITRS geocentric position vector, AU */

  double gpearth;          /* Earth gravitational potential at observer / c^2 */

  double refco[NREFCO];    /* refraction coefficients */

  /* Time */
  double tt;
  double tdb;
  double dtdb;

  /* Earth orientation parameters from IERS tabulations */
  double dut1;             /* TT-UT1, s */
  double xp;               /* X polar motion, rad */
  double yp;               /* Y polar motion, rad */

  double pmm[3][3];        /* Polar motion matrix */
  double lpm[3][3];        /* Longitude and polar motion matrix */

  double tigop[3];         /* TIRS geocentric position vector, AU */
  double tigov[3];         /* TIRS geocentric velocity vector, AU/d */
  double tigoa[3];         /* TIRS geocentric acceleration vector, AU/d/d */

  /* Frame bias, precession and nutation */
  double ang_pfb[NPNANG];  /* Frame bias and precession */
  double ang_pnfb[NPNANG]; /* Frame bias, precession and nutation */

  double dpsi;             /* Nutation in longitude, rad */
  double deps;             /* Nutation in obliquity, rad */
  double sxy;              /* CIO locator + XY/2 */

  double cim[3][3];        /* GCRS to CIRS matrix */

  /* Earth rotation and derived quantites */
  double era;              /* ERA (rad) */
  double erm[3][3];        /* Earth rotation matrix */
  double dermdt[3][3];     /* Time derivative (/d) */

  double ctm[3][3];        /* Full GCRS - Topocentric matrix */
  double dctmdt[3][3];     /* Time derivative (/d) */

  double gop[3];           /* GCRS geocentric position vector, AU */
  double gov[3];           /* GCRS geocentric velocity vector, AU/d */
  double goa[3];           /* GCRS geocentric acceleration vector, AU/d/d */

  /* Solar-system ephemerides */
  double bep[3];           /* SSB to Earth, AU */
  double bev[3];           /* Earth velocity relative to SSB, AU/d */
  double hep[3];           /* Heliocentre to Earth, AU */
  double hev[3];           /* Earth velocity relative to Heliocentre, AU/d */

  double bsp[3];           /* SSB to Sun, AU */
  double bsv[3];           /* Sun velocity relative to SSB, AU/d */
  double emp[3];           /* Earth to Moon, AU */
  double emv[3];           /* Moon velocity relative to Earth, AU/d */

  /* Quantities for parallax, light deflection, aberration */
  double bop[3];           /* Barycentre to observer, AU */
  double hop[3];           /* Heliocentre to observer, AU */

  double vab[3];           /* Observer velocity relative to SSB / c */
  double aab[3];           /* Observer acceleration relative to SSB / c */
  double gab;              /* Lorentz factor */

  double hdist;            /* Modulus of "hop" vector, AU */
  double gpsun;            /* Gravitational potential of Sun at Geocentre / c^2 */
  double minbdefl;         /* Minimum impact parameter for deflection, AU */
};

struct source {
  /* Name (not touched at all internally, use is optional) */
  char name[128];

  /* What type of source is it?  The JPLEPH constants defined above
     are used for sources covered by the JPL ephemerides, so these
     constants must not overlap. */
  unsigned char type;
#define SOURCE_STAR       16
#define SOURCE_ELEM       17

  /* For stars, unit vector and velocity at catalogue epoch.  For elements,
     heliocentric position and velocity at reference epoch. */
  double ref_n[3];
  double ref_dndt[3];

  /* Reference epoch.  For stars, this is the epoch of the catalogue.
     For elements, epoch of osculation or perihelion passage. */
  double ref_tdb;

  /* Parallax or 1/distance */
  double pr;        /* 1 / AU */

  /* Universal orbital elements */
  double mu;        /* G(M_1+M_2), AU^3 / day^2 */
  double nn;        /* mean daily motion */
  double alpha;     /* total energy, mu/a, AU^2 / day^2 */
  double sqrtalpha; /* sqrt(|alpha|) */
  double rref;      /* heliocentric distance at reference epoch (ref_tdb), AU */
  double muorref;   /* mu / rref */
  double rvref;     /* scalar product of position and velocity vectors at ref epoch */
  double psi;       /* universal eccentric anomaly, last value to speed up iterations */
};

/* -- Miscellaneous useful macros and inline functions -- */

/* Inline sine and cosine in one operation, for arguments already
   in or reduced to appropriate range (usually using fmod).  In the
   library, this is called with arguments in (-TWOPI,TWOPI), but on
   x87 the allowed range for the fsincos instruction is much larger
   than this.  The inline routine allows us to avoid extra floating
   point load and store operations, which are not particularly quick
   on some x87 implementations, and a branch for argument reduction.
   It is intended for use with the GNU C compiler set to -ffast-math,
   which should inline many math.h functions, but doesn't seem to
   know about fsincos (even if one uses the sincos GNU extension
   provided by glibc, at least on my system). */
#if defined(__GNUC__) && (defined(__i386) || defined(__amd64))
/* x86 assembler routine */
#define rdsincos(a, s, c) __asm__ ("fsincos" : "=t" (c), "=u" (s) : "0" (a))
#else
/* Generic C implementation using library functions */
#define rdsincos(a, s, c) {			\
  (s) = sin(a);					\
  (c) = cos(a);					\
}
#endif

/* Wrap angle to [0, TWOPI) */
#define ranorm(a) ((a) >= 0 ? fmod((a), TWOPI) : TWOPI+fmod((a), TWOPI))

/* Wrap angle to (-PI, PI] */
#define range(a) remainder((a), TWOPI)

/* Evaluate multiple polynomials simultaneously using Horner's method.  The
   array "p" is two dimensional, the most rapidly varying dimension is "nout"
   and the less rapid is "ncoeff" - i.e. the coefficients for the same
   degree in the outputs are next to each other.  Degree decreases through
   the array.  The weird arrangement of the coefficients and loops are to
   make sure it can be vectorized.  Needs to be inline so the compiler
   notices the alignment and padding requirements on the array "p" when it
   is created.  The optimizer should unroll both loops, or at least the
   inner one. */
#define sum_poly(x, p, nout, ncoeff, r) {	\
  int iout, icoeff, i;				\
						\
  for(i = 0; i < (nout); i++)			\
    (r)[i] = (p)[i];				\
						\
  for(icoeff = 1; icoeff < (ncoeff); icoeff++)	\
    for(iout = 0; iout < (nout); iout++, i++)	\
      (r)[iout] = (r)[iout] * (x) + (p)[i];	\
}

/* -- airmass.c: Airmass -- */

double airmass (double secz);

/* -- bary.c: Barycentering -- */

double bary_delay (struct observer *obs, double s[3], double pr);
double bary_doppler (struct observer *obs, double s[3]);

/* -- dtai.c: TAI-UTC -- */

int dtai_read (struct dtai_table *tab, char *filename);
void dtai_free (struct dtai_table *tab);
double dtai (struct dtai_table *tab, int mjd, double frac);

/* -- dtdb.c: TDB-TT -- */

double dtdb (double mjd,        /* TDB as MJD (but TT is adequate) */
	     double dut1,       /* TT-UT1, s */
	     double longitude,  /* rad, East + */
	     double u,          /* distance from Earth spin axis, m */
	     double z);         /* distance from Earth equatorial plane, m */

/* -- dsincos.c: sine and cosine in one operation -- */

void dsincos (double a, double *s, double *c);

/* -- fpcoord.c: focal plane coordinates -- */

int vec2tp (double s[3], double tp[3], double *x, double *y);
void tp2vec (double x, double y, double tp[3], double s[3]);

/* -- geoc.c: geodetic to geocentric -- */

void geoc (double sinphi,     /* sin, cos geodetic latitude */
	   double cosphi,
	   double height,     /* above geoid, m */
	   double *u,         /* returned: distance from Earth spin axis, m */
	   double *z);        /* returned: distance from equatorial plane, m */

/* -- jpleph.c: JPL ephemerides -- */

int jpleph_open (struct jpleph_table *p, int type, char *filename);
void jpleph_close (struct jpleph_table *p);
int jpleph_fetch (struct jpleph_table *p, double mjd, int body,
		  double *pos, double *vel);

/* -- iers.c: TT-UT1 and polar motion -- */

int iers_open (struct iers_table *tab,
	       struct dtai_table *dtai_tab,
	       char *filename);
void iers_close (struct iers_table *tab);
int iers_fetch (struct iers_table *tab, double mjd,
		double *dut1, double *xp, double *yp, double *dxnut, double *dynut);

/* -- matvec.c: matrix and vector utilities */

/* Makes Euler rotation matrices about "axis" */
void euler_rotate (double m[3][3], int axis, double angle);
void euler_rotate_sc (double m[3][3], int axis, double sa, double ca);

/* Identity matrix */
#define m_identity(m) {				\
  (m)[0][0] = (m)[1][1] = (m)[2][2] = 1.0;	\
  (m)[0][1] = (m)[1][0] = (m)[2][1] = 0.0;	\
  (m)[0][2] = (m)[1][2] = (m)[2][0] = 0.0;	\
}

/* Matrix transpose */
#define m_transpose(a, b) {			\
  (b)[0][0] = (a)[0][0];			\
  (b)[0][1] = (a)[1][0];			\
  (b)[0][2] = (a)[2][0];			\
  (b)[1][0] = (a)[0][1];			\
  (b)[1][1] = (a)[1][1];			\
  (b)[1][2] = (a)[2][1];			\
  (b)[2][0] = (a)[0][2];			\
  (b)[2][1] = (a)[1][2];			\
  (b)[2][2] = (a)[2][2];			\
}

/* v' = M v */
void m_x_v (double m[3][3], double vi[3], double vo[3]);

/* v' = M^T v */
void mt_x_v (double m[3][3], double vi[3], double vo[3]);

/* C = A B */
void m_x_m (double a[3][3], double b[3][3], double c[3][3]);

/* Basic vector operations */
#define v_copy(u, v) {				\
  (u)[0] = (v)[0];				\
  (u)[1] = (v)[1];				\
  (u)[2] = (v)[2];				\
}

#define v_p_sv(u, s, v) {		\
  (u)[0] += (s)*(v)[0];			\
  (u)[1] += (s)*(v)[1];			\
  (u)[2] += (s)*(v)[2];			\
}

/* Scalar and vector products */
#define v_d_v(u, v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])
#define v_x_v(u, v, p) {			\
  (p)[0] = (u)[1]*(v)[2] - (u)[2]*(v)[1];		\
  (p)[1] = (u)[2]*(v)[0] - (u)[0]*(v)[2];		\
  (p)[2] = (u)[0]*(v)[1] - (u)[1]*(v)[0];		\
}

/* Normalises a vector, returns 1 / norm. */
double v_renorm (double v[3]);

/* Vector from/to spherical (-ha, delta) or (az, el) */
void ad_to_v (double a, double d, double v[3]);
void v_to_ad (double v[3], unsigned char flip, double *a, double *d);

/* Vectors to time derivatives of spherical coordinates */
void v_to_ad_dt (double v[3], double dvdt[3], unsigned char flip,
		 double *dadt, double *dddt);

/* Vector to spherical (az, zd) */
void v_to_az (double v[3], unsigned char flip, double *a, double *z);

/* Angle between two vectors, robust method */
double v_angle_v (double u[3], double v[3]);

/* -- mjd.c: generic Julian date functions -- */

int date2mjd (int yr, int mn, int dy);
void mjd2date (int n, int *yr, int *mn, int *dy);

/* -- mount.c: vectors to/from general mount roll and pitch angles -- */

int mount_ab2rp (double *aim, double *daimdt,
		 double *bore,
		 double snp, double cnp,
		 unsigned char flip,
		 double pos[3][3], double *r, double *p,
		 double dposdt[3][3], double *drdt, double *dpdt);

void mount_pa (double *aimp, double *daimpdt,
	       double *bore,
	       double pos[3][3], double dposdt[3][3],
	       double *a, double *dadt);

void mount_rp2pos (double r, double p,
		   double snp, double cnp,
		   double pos[3][3]);

/* -- mpc.c: routines to read MPC 1-line format -- */

int mpc_read (char *filename, struct source **srclist, int *nsrc);
int mpc_convert (char *line, struct source *src);

/* -- observer.c: observer structure setup and update -- */

void observer_init (struct observer *obs,
		    struct jpleph_table *jpltab,
		    struct jpleph_table *tetab,
		    struct iers_table *iertab,
		    double longitude,  /* WGS84 longitude, E positive */
		    double latitude,   /* WGS84 latitude, N positive */
		    double height);    /* Height above geoid, m */

#define OBSERVER_UPDATE_IERS   0x01
#define OBSERVER_UPDATE_NUT    0x02
#define OBSERVER_UPDATE_PFB    0x04
#define OBSERVER_UPDATE_ERA    0x08
#define OBSERVER_UPDATE_TDB    0x10
#define OBSERVER_UPDATE_SOLSYS 0x20

#define OBSERVER_UPDATE_ALL    0xff

int observer_update (struct observer *obs,
		     double tt,
		     unsigned char mask);

void observer_ast2obs (struct observer *obs,
		       double *s,
		       double *dsdt,
		       double pr,
		       unsigned char mask);

void observer_obs2ast (struct observer *obs,
		       double *s,
		       double pr,
		       unsigned char mask);

/* -- parallactic.c: parallactic angle -- */

double parallactic (double sinphi, double cosphi, double v[3]);

/* -- prenut.c: precession and nutation -- */

/* Precession and frame bias, IAU 2006, Fukushima-Williams angles */
void pfb06ang (double jctk, double ang[NPNANG]);

/* Nutation, IAU 2000B, and optionally, the CIO locator s+XY/2 and
   equation of the origins */
void nut00b (double jctk,
	     double *dpsi, double *deps,
	     double *sxy, double *eo);

/* Form traditional frame bias - precession - nutation matrix */
void makepnm (double ang[NPNANG],
	      double m[3][3]);

/* Form GCRS to CIRS matrix */
void makecim (double ang[NPNANG], double sxy,
	      double dxnut, double dynut,  /* nutation correction, from IERS Bull. */
	      double m[3][3]);

/* -- refract.c: atmospheric refraction -- */

/* Computes refraction constants using Saastamoinen (1972) approximate
   treatment.  Accuracy to about 0.1-0.2 arcsec down to 80 degrees ZD
   for most optical purposes. */
void refract_const (double temperat,   /* K */
		    double humidity,   /* 0-1 */
		    double pressure,   /* mbar */
		    double wavelength, /* um */
		    double height,     /* m */
		    double *refco);

/* Computes correction for refraction observed -> vacuum ZD, and derivatives
   with respect to ZD (for Newton-Raphson, etc.) */
void refract_corr (double *refco, double tanz, double *refr, double *deriv);

/* Apply refraction corrections to cartesian vector "vi".  If unref=0,
   converts vacuum ZD to observed ZD; if undef=1 converts observed ZD to
   vacuum ZD (this is the conventional way round in the refraction
   formula).  Correction is held back within about 3 degrees of the
   horizon to prevent overflows.  The treatment should be better than
   the refraction model itself above this point.  It is fine to pass the
   same vector for both arguments, in which case the input will be
   overwritten by the output. */
void refract_vec (double *refco, unsigned char unref,
		  double *vi, double *vo,
		  double *dvidt, double *dvodt);

/* -- source.c -- */

void source_init (void);

void source_star (struct source *src,
		  double ra, double de,      /* radians */
		  double pmra, double pmde,  /* sky projected, arcsec/yr */
		  double plx, double vrad,   /* arcsec, km/s */
		  double epoch);

#define SOURCE_ELEM_MAJOR  1
#define SOURCE_ELEM_MINOR  2
#define SOURCE_ELEM_COMET  3

void source_elem (struct source *src,
		  unsigned char eltype,
		  double epoch,
		  double incl,
		  double anode,
		  double longperi,
		  double aq,
		  double ecc,
		  double lm,
		  double nn);

void source_place (struct observer *obs,
		   struct source *src,
		   unsigned char mask,   /* parts we want */
		   double *s,
		   double *dsdt,
		   double *pr);

/* -- strutil.c: string utilities -- */

char *extractstr (char *str, int len,
		  int f, int l);
int extractdouble (char *str, int len,
		   int f, int l,
		   double *val);
int extractint (char *str, int len,
		int f, int l,
		int *val);
int extractintfrac (char *str, int len,
		    int f, int l,
		    int *ival, double *fval);
char *sstrip (char *str);

/* -- stumpff.c: Stumpff functions */

void stumpff (double s, double alpha, double sqrtalpha, double *c);

#endif  /* LFA_H */
