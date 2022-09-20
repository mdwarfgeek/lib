/* lfa.h: positional astronomy subroutines and other bits.
          various source files, see below. */

#ifndef LFA_H
#define LFA_H

#include <stdio.h>  /* for FILE */
#include <stdint.h>
#include <float.h>
#include <math.h>

#if defined(__GLIBC__)
#include <features.h>
#endif

#if defined(_MSC_VER)
#include <malloc.h>
#endif

/* -- Physical, Mathematical and Astronomical constants -- */

/* NOTE: integer constants are expressed here without the .0, so be
   careful when using them if the result needs to be a floating point
   number. */

/* Fundamental (defining) constants, mostly IAU and IERS */

#define GMSUN  1.32712440041e20  /* m^3 / s^2, IAU 2009, TDB compatible */
#define AU     149597870700.0    /* m, IAU 2009 system */
#define LIGHT  2.99792458e8      /* m/s, definition */

#define GMEARTH 3.986004356e14     /* m^3 / s^2, IAU 2009, TDB compatible */
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

#define RSUN   6.957e8           /* m, IAU 2015 Resol B3 conversion constant */
#define RJUP   7.1492e7          /* m, equ., IAU 2015 Resol B3 */
#define REARTH 6.3781e6          /* m, equ., IAU 2015 Resol B3 */

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

/* -- Miscellaneous macros -- */

/* Make a variable length array on the stack.  Uses _alloca on MSVC
   to patch around lack of support for C99 variable length automatic
   arrays.  Assumes C99 support on all other compilers. */
#if defined(_MSC_VER)
#define VLAONSTACK(t, i, n) t *i = (t *) _alloca((n) * sizeof(t)) 
#else
#define VLAONSTACK(t, i, n) t i[n]
#endif

/* -- Data structures -- */

struct jpleph_table {
  FILE *fp;
  size_t (*read) (void *, size_t, size_t, FILE *);

  int32_t denum;
  int32_t ncoeff;
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
#define JPLEPH_EULER     13

/* JPL time ephemeris integral, giving the quantity
   TT-TDB in seconds.  Note the sign compared to
   TIMEEPH_TEI below, and that the constant ZTDB and
   the conversion to TDB (the factor 1/(1-LC)) are
   already included. */
#define JPLEPH_TEI       14

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
  uint8_t has_time;
  double lc;

  double *buf;
  int32_t brec;
};

struct dtai_entry {
  int mjd;  /* leap seconds can only happen at midnight so always integer */
  double dtai;
  double mjdzero;
  double scale;
};

struct dtai_table {
  struct dtai_entry *table;
  int ilast;
  int ntab;
};

struct fsgp_fac {
  /* Inputs */
  double *kern;
  double sumaj;
  double *t;
  double *yvar;
  
  /* Factorization */
  double *phi;
  double *u;  /* ubar */
  double *v;  /* vbar */
  double *w;  /* wbar */
  double *d;
  double *sqrtd;

  /* Sizes */
  int nkern;
  int nckern;
  int ndp;
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
  double latitude;
  double longitude;
  double height;

  double lmm[3][3];        /* Observer longitude matrix */
  double phm[3][3];        /* Observer latitude matrix */
  double itgop[3];         /* ITRS geocentric position vector, AU */

  double gpearth;          /* Earth gravitational potential at observer / c^2 */

  double refco[NREFCO];    /* refraction coefficients */

  /* Time */
  double utc;
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
  double jbep[3];          /* original values and timestamp from JPLEPH */
  double jemp[3];          /* used to update estimates below */
  double jbsp[3];
  double jtdb;

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
  uint8_t type;
#define SOURCE_STAR       16
#define SOURCE_ELEM       17

  /* Pad to 64 bits explicitly */
  uint8_t _pad[7];

  /* For stars, unit vector and velocity at catalogue epoch.  For elements,
     heliocentric position and velocity at reference epoch. */
  double ref_n[3];
  double ref_dndt[3];

  /* Reference epoch.  For stars, this is the epoch of the catalogue.
     For elements, epoch of osculation or perihelion passage. */
  double ref_tdb;

  /* Parallax or 1/distance */
  double pr;        /* 1 / AU */

  /* Doppler factor */
  double df;

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

struct rng_state {
#define RNG_RR 191
#define RNG_RQ (RNG_RR*2)
#define RNG_RL (RNG_RR*4)

#define RNG_NR (RNG_RR+1)
#define RNG_NQ (RNG_NR*2)
#define RNG_NL (RNG_NR*4)

  union {
    uint64_t i[RNG_NQ];
    double   d[RNG_NQ];
  } a;

  uint64_t *ip;
  double *dp;
  int r;

  double gauss;
  unsigned char havegauss;
  unsigned char haveieee;
};

struct wcs_info {
  double a;
  double b;
  double c;
  double d;
  double e;
  double f;

  double tpa;
  double tpd;
  double sina;
  double cosa;
  double sind;
  double cosd;

  enum {
    PROJ_TAN,
    PROJ_SIN,
    PROJ_ARC
  } proj;

#define PC_MAX 20       /* defined in Calabretta et al. (2002) */
  double pc[PC_MAX+1];
  int mpc;              /* last populated index in pc array */
};

/* -- Inline routines -- */

/* inline_sincos(a, s, c)
   inline_bare_sincos(a, s, c)

   Inline sine and cosine in one operation.  The second routine is for
   arguments already in the appropriate range.  In the library, this
   is called with arguments in (-TWOPI,TWOPI), but on x87 the allowed
   range for the fsincos instruction is much larger than this. */
   
#if defined(__GLIBC__) && defined(_GNU_SOURCE)

/* Try to use glibc provided routines in math.h */
#define inline_sincos(a, s, c) sincos((a), &(s), &(c))
#define inline_bare_sincos(a, s, c) sincos((a), &(s), &(c))

#else

/* Generic implementation using standard math.h functions */
#define inline_sincos(a, s, c) (s) = sin(a); (c) = cos(a)
#define inline_bare_sincos(a, s, c) (s) = sin(a); (c) = cos(a)

#endif

/* ilogtwo(i, l)
   Computes l = floor(log2(i)).  The strange name is to avoid
   conflicting with various other implementations (e.g. NetBSD). */

#if defined(__GNUC__) && (defined(__i386) || defined(__amd64))

/* Inline x86 assembler routine, GNU C syntax */
#define ilogtwo(i, l) __asm__ ("bsr %1, %0" : "=r" (l) : "mr" (i))

#elif defined(__GNUC__) && (__GNUC__ >= 3 && __GNUC_MINOR__ >= 2 && __GNUC_PATCHLEVEL__ >= 2)

/* Use GCC count leading zeros builtin */
#define ilogtwo(i, l) (l) = 8*sizeof(i) - 1 - __builtin_clz(i)

#elif FLT_RADIX == 2

/* Using ilogb is probably better than a loop provided we have
   hardware floating point. */

#define ilogtwo(i, l) (l) = ilogb((double) (i))

#else

/* Weird architecture.  Note that this depends on automatic truncation. */
#define ilogtwo(i, l) (l) = log2((double) (i))

#endif

/* Wrap angle to [0, TWOPI).  The second fmod for negative cases prevents
   very small negative numbers from wrapping to +TWOPI. */
#define ranorm(a) ((a) >= 0 ? fmod((a), TWOPI) : fmod(fmod((a), TWOPI) + TWOPI, TWOPI))

/* Wrap [-PI, PI] angle to [0, TWOPI) */
#define ranormp(a) fmod((a) + TWOPI, TWOPI)

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

double v_airmass (double v[3]);

/* -- background.c: image background removal for object detection -- */

int backremove (float *mapin, unsigned char *mask, float *mapout,
                int nx, int ny, int nbsize);

/* -- bary.c: Barycentering -- */

double bary_delay (struct observer *obs, double s[3], double pr);
double bary_doppler (struct observer *obs, double sref[3],
                     double s[3], double dsdt[3], double pr);

/* -- cholesky.c: Cholesky decomposition -- */

int cholesky (double *a, int n);

/* -- const.c: Precomputed matrices -- */

extern double gcrs2ecl[3][3];
extern double ecl2gcrs[3][3];
extern double fk52ecl[3][3];
extern double ecl2fk5[3][3];
extern double eq2gal[3][3];
extern double gal2eq[3][3];

/* -- dict.c: Simple dictionary (hash table) */

struct dict_entry {
  uint32_t hash;
  uint32_t keylen;
  void *key;
  void *value;
  struct dict_entry *next;
};

struct dict {
  struct dict_entry **list;
  uint32_t nalloc;
  uint32_t nused;
};

int dict_init (struct dict *d, uint32_t nalloc);
void dict_free (struct dict *d);
int dict_delete (struct dict *d,
                 void *key, uint32_t keylen);
int dict_fetch (struct dict *d,
                void *key, uint32_t keylen,
                void **value);
int dict_store (struct dict *d,
                void *key, uint32_t keylen,
                void *value);

/* -- dplate.c: Linear transformation between 2-D coordinate systems */

int dplate (void *comxptr, size_t comxoff, size_t comxsz,
            void *comyptr, size_t comyoff, size_t comysz,
            void *refxptr, size_t refxoff, size_t refxsz,
            void *refyptr, size_t refyoff, size_t refysz,
            void *wtptr, size_t wtoff, size_t wtsz,
            int npt, int ncoeff, double *tr);

#define dplate_tr_x(x, y, tr) ((tr)[0]*(x) + (tr)[1]*(y) + (tr)[2])
#define dplate_tr_y(x, y, tr) ((tr)[3]*(y) + (tr)[4]*(x) + (tr)[5])

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

/* -- filt1d.c: standard 1-D filters -- */

int filt1d_nwork (int npt, int nkern);
void filt1d_boxcar (float *buf, unsigned char *mask, float *work,
                    int npt, int nstride, int nkern);
void filt1d_hanning (float *buf, unsigned char *mask, float *work,
                     int npt, int nstride);
void filt1d_kernel (float *buf, unsigned char *mask, float *work,
                    int npt, int nstride,
                    float *kern, int nkern);
void filt1d_median (float *buf, unsigned char *mask, float *work,
                    int npt, int nstride, int nkern);

/* -- fpcoord.c: focal plane coordinates -- */

int vec2tp (double s[3], double tp[3], double *x, double *y);
void tp2vec (double x, double y, double tp[3], double s[3]);

/* -- fsgp.c: "fast and scalable Gaussian processes" -- */

int fsgp_kern_sho (double s0, double w0, double q, double *kern);
int fsgp_kern_matern (double s0, double w0, double f, double *kern);
int fsgp_kern_valid (double *kern, int nkern);
int fsgp_compute (struct fsgp_fac *fac,
                  double *kern, int nkern,
                  double *t, double *yerr, int ndp);
int fsgp_predict (struct fsgp_fac *fac, double *y,
                  double *tpred, double *ypred, double *varpred, int npred);
int fsgp_residual (struct fsgp_fac *fac,
                   double *y, double *z, int nrhs);
int fsgp_sample (struct fsgp_fac *fac, double *q, double *y);
double fsgp_logdet (struct fsgp_fac *fac);
int fsgp_loglike (struct fsgp_fac *fac, double *y, double *loglike);
void fsgp_free (struct fsgp_fac *fac);

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

/* -- kepler.c: solve Kepler's equation -- */

double kepler (double ma, double ecc);

/* -- lrmatch.c: likelihood ratio catalogue matching -- */

int lrmatch (double *comx, double *comy,
             double *comlogrank, double *comerr, int ncom,
             double *refx, double *refy,
             double *reflogrank, double *referr, int nref,
             double searchrad, int sorted_y,
             int *best_ref_for_com, int *best_com_for_ref);

/* -- matvec.c: matrix and vector utilities */

/* Makes Euler rotation matrices about "axis" */
void euler_rotate (double m[3][3], int axis, double angle);
void euler_rotate_sc (double m[3][3], int axis, double sa, double ca);

/* Makes derivative of Euler rotation matrix about "axis" */
void euler_drot_sc (double m[3][3], int axis, double dsa, double dca);

/* Identity matrix */
#define m_identity(m) {				\
  (m)[0][0] = (m)[1][1] = (m)[2][2] = 1.0;	\
  (m)[0][1] = (m)[1][0] = (m)[2][1] = 0.0;	\
  (m)[0][2] = (m)[1][2] = (m)[2][0] = 0.0;	\
}

/* Matrix copy */
#define m_copy(a, b) {				\
  (b)[0][0] = (a)[0][0];			\
  (b)[0][1] = (a)[0][1];			\
  (b)[0][2] = (a)[0][2];			\
  (b)[1][0] = (a)[1][0];			\
  (b)[1][1] = (a)[1][1];			\
  (b)[1][2] = (a)[1][2];			\
  (b)[2][0] = (a)[2][0];			\
  (b)[2][1] = (a)[2][1];			\
  (b)[2][2] = (a)[2][2];			\
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
		    double longitude,  /* WGS84 longitude, E positive */
		    double latitude,   /* WGS84 latitude, N positive */
		    double height);    /* Height above geoid, m */

void observer_geoc (struct observer *obs);

#define OBSERVER_UPDATE_IERS   0x01
#define OBSERVER_UPDATE_NUT    0x02
#define OBSERVER_UPDATE_PFB    0x04
#define OBSERVER_UPDATE_ERA    0x08
#define OBSERVER_UPDATE_TDB    0x10
#define OBSERVER_UPDATE_SOLSYS 0x20

#define OBSERVER_UPDATE_ALL    0xff

int observer_update (struct observer *obs,
		     struct jpleph_table *jpltab,
		     struct jpleph_table *tetab,
		     struct iers_table *iertab,
		     double utc,
		     double ttmutc,
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

double v_parallactic (double sinphi, double cosphi, double v[3]);

/* -- pixovcirc.c: fraction of pixel inside circle -- */

double pixovcirc (double x, double y, double r);

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

/* Form GCRS to mean ecliptic of date matrix */
void makeeclm (double ang[NPNANG],
               double m[3][3]);

/* Form GCRS to CIRS matrix */
void makecim (double ang[NPNANG], double sxy,
	      double dxnut, double dynut,  /* nutation correction, from IERS Bull. */
	      double m[3][3]);

/* -- qr.c: QR decomposition and solution of linear equations -- */

void qr (double *a, double *s, double *betaarr, int *perm, int n);
int qrsolve (double *a, double *s, double *betaarr, int *perm,
             double *b, int n, double rcond);
int qrinvert (double *a, double *s, double *betaarr, int *perm,
              double *ainv, int n, double rcond);

int linsolve (double *a, double *b, int n, double rcond);

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

/* -- rng.c: random number generation -- */

/* User level routines: call rng_init() to initialize and seed the
   generator, and then rng_fetch() to fetch n random numbers.  The
   distribution is uniform in [0,1).  Order of the sequence is
   guaranteed regardless of how many numbers are fetched in each
   call.  For efficiency, fetch as many numbers at a time as
   possible, preferably multiples of RNG_RQ.  Internal buffering
   is used where necessary so we can satisfy any size of request. */
void rng_init (struct rng_state *s, uint32_t seed);
void rng_fetch_uniform (struct rng_state *s, double *a, int n);
double rng_fetch_one_uniform (struct rng_state *s);
void rng_fetch_gauss (struct rng_state *s, double *a, int n);
double rng_fetch_one_gauss (struct rng_state *s);
int rng_fetch_mvgauss (struct rng_state *s,
                       double *mean, double *cov,
                       double *work,
                       double *ans, int n);

/* -- shuffle.c -- */

void shuffle (struct rng_state *rs,
              void *list, void *tmp,
              size_t n, size_t s);

/* -- skylevel.c -- */

#define SKYLEVEL_DEFAULT_OFFSET 1024
#define SKYLEVEL_DEFAULT_ULIM   66559
#define SKYLEVEL_DEFAULT_SIZE   66560

void skylevel (int *ihist, int ihmin, int ihmax, int mpix,
	       float clip_low, float clip_high,
	       float *skylev_r, float *sigma_r);

void skylevel_image (float *map, unsigned char *mask, int npix,
                     float clip_low, float clip_high,
                     float *skylev, float *skynoise);

/* -- sort.c: sorting and selection -- */

/* Select the k'th smallest of il..ir elements of size s, each 
   containing data of the given datatype at offset o.  The offset
   within a structure can be determined using offsetof() from
   stddef.h.  Returns a pointer to the element.  */

double dquickselect (double *list, size_t k, size_t n);
float fquickselect (float *list, size_t k, size_t n);
int iquickselect (int *list, size_t k, size_t n);

void dmultquickselect (double *list, size_t n,
                       size_t *ind, size_t nind,
                       double *result);
void fmultquickselect (float *list, size_t n,
                       size_t *ind, size_t nind,
                       float *result);
void imultquickselect (int *list, size_t n,
                       size_t *ind, size_t nind,
                       int *result);

/* Sort n elements of size s, using recursion on the call stack. */

void dquicksort (double *list, size_t n);
void fquicksort (float *list, size_t n);
void iquicksort (int *list, size_t n);

void dquicksort_gen (void *list, void *tmp,
                     size_t n, size_t s, size_t o);
void fquicksort_gen (void *list, void *tmp,
                     size_t n, size_t s, size_t o);
void iquicksort_gen (void *list, void *tmp,
                     size_t n, size_t s, size_t o);

/* My usual "medsig" routine using median and MAD scaled to Gaussian
   RMS equivalent.  Implemented using selection internally, and
   doesn't take the mean of the middle two elements if n is even,
   unlike the original routines. */

void dmedsig (double *list, size_t n, double *median_r, double *sigma_r);
void fmedsig (float *list, size_t n, float *median_r, float *sigma_r);

/* -- source.c -- */

void source_init (void);

void source_star (struct source *src,
		  double ra, double de,      /* radians */
		  double pmra, double pmde,  /* sky projected, arcsec/yr */
		  double plx, double vrad,   /* arcsec, km/s */
		  double epoch);

void source_star_vec (struct source *src,
		      double *n,
		      double *dndt,
		      double pr,
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
		   struct jpleph_table *jpltab,
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

/* -- stumpff.c: Stumpff functions -- */

void stumpff (double s, double alpha, double sqrtalpha, double *c);

/* -- sysinfo.c: System information -- */

/* Returns number of CPUs available for threading, or -1 on error. */
int get_num_cpus (void);

/* Returns total system RAM in bytes, or 0 on error.  Please use this
   routine with discretion and common sense engaged.  It's intended for
   scaling block sizes in image processing routines that support
   partitioning the workload, to prevent them running into swap. */
uint64_t get_total_mem (void);

/* -- wcs.c: Support for a limited subset of FITS-WCS -- */

void wcs_vec2xy (struct wcs_info *wcs, double *vec, double *x, double *y);
void wcs_xy2vec (struct wcs_info *wcs, double x, double y, double *vec);
void wcs_ad2xy (struct wcs_info *wcs, double a, double d, double *x, double *y);
void wcs_xy2ad (struct wcs_info *wcs, double x, double y, double *a, double *d);
void wcs_xy2xy (struct wcs_info *wcs1, struct wcs_info *wcs2,
		double x1, double y1, double *x2, double *y2);
void wcs_tp2xy (struct wcs_info *wcs, double fx, double fy,
		double *x, double *y);
void wcs_xy2tp (struct wcs_info *wcs, double x, double y,
		double *fx, double *fy);

#endif  /* LFA_H */
