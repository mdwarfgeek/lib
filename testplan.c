#include <sys/types.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/select.h>
#endif
#include <sys/time.h>
#include <unistd.h>
#include <errno.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lfa.h"
#include "cvtunit.h"
#include "tcutil.h"
#include "util.h"

int main (int argc, char *argv[]) {
  struct dtai_table dtab;
  struct iers_table itab;
  struct jpleph_table jtab, ttab, *tptr = NULL;
  struct source src;
  struct observer obs;

  struct timeval tv, tvslow, tvdiff;
#ifndef _WIN32
  struct timeval tsl;
#endif
  int utcdn;
  double utcdf;

  double dtt;

  /* MEarth */
  double longitude = -399161.0*AS_TO_RAD;  /* -110 52 41.0 */
  double latitude  =  114063.1*AS_TO_RAD;  /*  +31 41 03.1 */
  double height = 2384.0;
  double temperat = 283.16;
  double humidity = 0.3;
  double pressure = 767.1;
  double wavelength = 0.85;

  char *names[] = { "Mercury", "Venus", "EMB", "Mars",
		    "Jupiter", "Saturn", "Uranus", "Neptune",
		    "Pluto", "Moon", "Sun" };

  double s[3], dsdt[3], a, d, dadt, dddt;
  char astr[64], dstr[64];
  double az, el;

  int rv, i;

  double pr;

  int doclr = 0, ncols, nrows, running = 0;

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

  /* Setup observer structure */
  observer_init(&obs, longitude, latitude, height);

  refract_const(temperat, humidity, pressure, wavelength, height,
		obs.refco);

  doclr = tcutil_init();
  if(doclr) {
    /* Check there is enough room for cursor jiggery pokery */
    tcutil_winsize(&ncols, &nrows);
    if(ncols < 62 || nrows < JPLEPH_SUN+2)
      doclr = 0;
  }

  tvslow.tv_sec = 0;
  tvslow.tv_usec = 0;

  for(;;) {
    /* Get UNIX time (UTC) */
    rv = gettimeofday(&tv, (struct timezone *) NULL);
    if(rv)
      error(1, "gettimeofday");

    /* Split into days since epoch and seconds since midnight.
       We want floor(tv_sec / DAY), computed here with integer
       arithmetic. */
    utcdn = tv.tv_sec / DAY;  /* truncation */
    if(tv.tv_sec < 0 && tv.tv_sec % DAY)
      utcdn--;  /* floor rather than truncating toward zero */

    /* Fraction, note subtle ordering to preserve precision */
    utcdf = ((tv.tv_sec - utcdn*DAY) +
	     tv.tv_usec / 1000000.0) / ((double) DAY);

    /* Convert to MJD */
    utcdn += JUNIX;

    /* Look up delta(TT) */
    dtt = dtai(&dtab, utcdn, utcdf) + DTT;

    /* How long since last update? */
    tvdiff.tv_sec = tv.tv_sec - tvslow.tv_sec;
    tvdiff.tv_usec = tv.tv_usec - tvslow.tv_usec;
    if(tvdiff.tv_usec < 0) {
      tvdiff.tv_sec--;
      tvdiff.tv_usec += 1000000;
    }

    /* Update */
    if(!running || tvdiff.tv_sec >= 10) {
      rv = observer_update(&obs, &jtab, tptr, &itab, utcdn+utcdf, dtt,
			   OBSERVER_UPDATE_ALL);
      if(rv)
	fatal(1, "observer_update: %d (%s)", rv, strerror(errno));

      tvslow = tv;
    }
    else {
      rv = observer_update(&obs, &jtab, tptr, &itab, utcdn+utcdf, dtt,
			   OBSERVER_UPDATE_PFB |
			   OBSERVER_UPDATE_ERA |
			   OBSERVER_UPDATE_SOLSYS);
    }

    a = obs.era+obs.longitude;
    a = ranorm(a);

    base10_to_60(utcdf*TWOPI, UNIT_RAD, astr, sizeof(astr), " ", "", 1, UNIT_HR);
    base10_to_60(a, UNIT_RAD, dstr, sizeof(dstr), " ", "", 1, UNIT_HR);

    if(doclr && running)
      for(i = -1; i <= JPLEPH_SUN; i++)
	tcutil_up();

    printf("Planet report for UTC %s LSA %s\n", astr, dstr);

#if 0
    printf("TT=%.10f TT-UT1=%.6f TDB-TT=%.6f XP=%.6f YP=%.6f\n",
	   tt, obs.dut1, obs.dtdb,
	   obs.xp*RAD_TO_AS, obs.yp*RAD_TO_AS);
#endif

    tcutil_attr(ATTR_BOLD);
    printf("Name    RA (Astrometric ICRS) DEC  d(RA)/dt d(DEC)/dt dist/AU el/d ");
    tcutil_attr(ATTR_NORM);
    printf("\n");

    for(i = 0; i <= JPLEPH_SUN; i++) {
      if(i == JPLEPH_EMB)
	continue;
      
      src.type = i;
      
      /* Astrometric place of source */
      source_place(&obs, &src, &jtab, TR_MOTION | TR_PLX, s, dsdt, &pr);

      if(i != JPLEPH_SUN) {
	/* Take off actual gravitational deflection, and put on
	   false deflection appropriate for star coordinates. */
	observer_ast2obs(&obs, s, dsdt, pr, TR_DEFL);
	observer_obs2ast(&obs, s, 0, TR_DEFL);
      }

      v_to_ad(s, 0, &a, &d);
      v_to_ad_dt(s, dsdt, 0, &dadt, &dddt);

      if(a < 0)
	a += TWOPI;

      base10_to_60(a, UNIT_RAD, astr, sizeof(astr), " ", "", 3, UNIT_HR);
      base10_to_60(d, UNIT_RAD, dstr, sizeof(dstr), " ", "+", 2, UNIT_DEG);

      /* Observed place */
      observer_ast2obs(&obs, s, dsdt, 0, TR_TO_OBS_AZ);
      v_to_ad(s, 0, &az, &el);
      
      printf("%-7s %s %s %9.6f %9.6f %7.4f %5.1f\n",
	     names[i], astr, dstr,
	     dadt*RAD_TO_AS/DAY, dddt*RAD_TO_AS/DAY,
	     1.0/pr, el*RAD_TO_DEG);
    }

#ifdef _WIN32
    Sleep(100);
#else
    tsl.tv_sec = 0;
    tsl.tv_usec = 100000;

    rv = select(0, (fd_set *) NULL, (fd_set *) NULL, (fd_set *) NULL, &tsl);
    if(rv < 0)
      error(1, "select");
#endif

    running = 1;
  }

  return(0);
}
