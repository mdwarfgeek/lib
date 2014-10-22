#include <stdio.h>
#include <math.h>

#include "lfa.h"

#include "cvtunit.h"
#include "util.h"

int main (int argc, char *argv[]) {
  struct dtai_table dtab;
  struct iers_table itab;
  struct jpleph_table jtab, ttab, *tptr = NULL;
  struct source src;
  struct observer obs;

  double utc = 56470.3;
  int uday;
  double dtt;

  /* MEarth */
  double longitude = -399161.0*AS_TO_RAD;  /* -110 52 41.0 */
  double latitude  =  114063.1*AS_TO_RAD;  /*  +31 41 03.1 */
  double height = 2384.0;
  double temperat = 283.16;
  double humidity = 0.3;
  double pressure = 767.1;
  double wavelength = 0.85;

  /* Barnard's star, from Hipparcos 2007 reduction */
  double catra = 269.45402263 * DEG_TO_RAD;
  double catde =   4.66828781 * DEG_TO_RAD;
  double pmra  =  -798.58 / 1000;
  double pmde  = 10328.12 / 1000;
  double plx   =   548.31 / 1000;
  double catep = 1991.25;

  /* CfA velocity from Dave Latham */
  double vrad  =  -110.3;

  double s[3], dsdt[3], a, d;
  char astr[64], dstr[64];

  int rv, i, n;

  double pr;
  double tbc, zbc;

  /* The entries here correspond to lines in the output.
     type=0 are RA/DEC, and type=2 are HA/DEC */
  struct {
    char *name;
    unsigned char mask;
  } ops[] = {
    { "Catalogue",    0 },
    { "BCRS",         TR_MOTION },
    { "Astrometric",  TR_MOTION | TR_PLX },
    { "Deflected",    TR_MOTION | TR_PLX | TR_DEFL },
    { "Local GCRS",   TR_MOTION | TR_PLX | TR_DEFL | TR_ANNAB },
    { "Topocentric",  TR_MOTION | TR_PLX | TR_DEFL | TR_ANNAB | TR_TOPO },
    { "Observed",     TR_MOTION | TR_PLX | TR_DEFL | TR_ANNAB | TR_TOPO | TR_REFRO },
  };

  /* Setup Earth orientation data and JPL ephemerides */
  rv = dtai_read(&dtab, (char *) NULL);
  if(rv)
    fatal(1, "dtai_open: error %d", rv);

  rv = iers_open(&itab, &dtab, (char *) NULL);
  if(rv)
    fatal(1, "iers_open: %d", rv);

  rv = jpleph_open(&jtab, 0, (char *) NULL);
  if(rv)
    fatal(1, "jpleph_open: %d", rv);

  if(!jtab.has_time) {
    rv = jpleph_open(&ttab, 1, (char *) NULL);
    if(rv == -2) {
      printf("Time ephemeris problem, continuing without\n");
      tptr = NULL;
    }
    else if(rv)
      fatal(1, "jpleph_open: %d", rv);
    else
      tptr = &ttab;
  }

  /* TT-UTC in seconds */
  uday = floor(utc);
  dtt = dtai(&dtab, uday, utc-uday) + DTT;

  /* Setup source structure */
  source_star(&src, catra, catde, pmra, pmde, plx, vrad, catep);

  /* Setup observer structure */
  observer_init(&obs, longitude, latitude, height);

  refract_const(temperat, humidity, pressure, wavelength, height,
		obs.refco);

  rv = observer_update(&obs, &jtab, tptr, &itab, utc, dtt, OBSERVER_UPDATE_ALL);
  if(rv)
    fatal(1, "observer_update_slow: %d", rv);

  /* Compute a few things */
  n = sizeof(ops) / sizeof(ops[0]);

  for(i = 0; i < n; i++) {
    source_place(&obs, &src, &jtab, ops[i].mask, s, dsdt, (double *) NULL);
    v_to_ad(s, 0, &a, &d);

    if(i >= 5) {
      /* Hour angle, so need to flip sign and put in [-pi,pi] */
      a = range(-a);  /* this is lazy and somewhat inefficient */
    }
    else {
      /* RA in [0,2*pi] */
      if(a < 0)
        a += TWOPI;
    }

    base10_to_60(a, UNIT_RAD, astr, sizeof(astr), " ", " ", 6, UNIT_HR);
    base10_to_60(d, UNIT_RAD, dstr, sizeof(dstr), " ", "+", 5, UNIT_DEG);

    printf("%-12s %s %s\n",
	   ops[i].name, astr, dstr);
  }

  /* Test of barycentric correction */
  source_place(&obs, &src, &jtab, TR_MOTION, s, dsdt, &pr);
  tbc = bary_delay(&obs, s, pr);
  zbc = bary_doppler(&obs, src.ref_n, s, dsdt, pr);

  printf("Bary Delta(T)  = %.10f Clock %.10f Total %.10f\n"
	 "Bary Vel Corr  = %.10f (z = %.15e)\n",
	 tbc,
	 dtt + obs.dtdb,
	 tbc + dtt + obs.dtdb,
         LIGHT * zbc, zbc);

  return(0);
}
