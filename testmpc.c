#include <sys/types.h>
#include <sys/time.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include "cvtunit.h"
#include "lfa.h"
#include "tcutil.h"
#include "util.h"

int main (int argc, char *argv[]) {
  struct dtai_table dtab;
  struct iers_table itab;
  struct jpleph_table jtab, ttab, *tptr = NULL;
  struct observer obs;

  struct timeval tv;
  int utcdn;
  double utcdf;

  double dtt, tt;

  /* MEarth */
  double longitude = -399161.0*AS_TO_RAD;  /* -110 52 41.0 */
  double latitude  =  114063.1*AS_TO_RAD;  /*  +31 41 03.1 */
  double height = 2384.0;
  double temperat = 283.16;
  double humidity = 0.3;
  double pressure = 767.1;
  double wavelength = 0.85;

  struct source *srclist = (struct source *) NULL, *src;
  int f, rv, isrc, nsrc = 0;

  double s[3], dsdt[3], a, d, dadt, dddt;
  char astr[64], dstr[64];
  double az, el, pr;

  if(argc < 2) {
    fprintf(stderr, "Usage:\t%s file [...]\n", argv[0]);
    exit(1);
  }

  tcutil_init();

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
  else
    tptr = NULL;

  /* Setup observer structure */
  observer_init(&obs, longitude, latitude, height);

  refract_const(temperat, humidity, pressure, wavelength, height,
		obs.refco);

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

  /* Easy to use but less precise TT as MJD */
  tt = utcdn + utcdf + dtt/DAY;

  printf("utc=%lf tt=%lf\n", utcdn+utcdf, tt);
  
  rv = observer_update(&obs, &jtab, tptr, &itab, utcdn+utcdf, dtt,
                       OBSERVER_UPDATE_ALL);
  if(rv)
    fatal(1, "observer_update: %d (%s)", rv, strerror(errno));

  for(f = 1; f < argc; f++) {
    rv = mpc_read(argv[f], &srclist, &nsrc);
    if(rv < 0)
      fatal(1, "mpc_read: %d", rv);
  }

  tcutil_attr(ATTR_BOLD);
  printf("Name             RA (Astrometric ICRS) DEC  d(RA)/dt d(DEC)/dt dist/AU el/d");
  tcutil_attr(ATTR_NORM);
  printf("\n");

  for(isrc = 0; isrc < nsrc; isrc++) {
    src = srclist + isrc;
    
#ifdef TEST
    fprintf(stderr,
            "%s: r_0=(%le, %le, %le) v_0=(%le, %le, %le) t_0=%lf alpha=%le r_0=%le\n",
            src->name,
            src->ref_n[0], src->ref_n[1], src->ref_n[2],
            src->ref_dndt[0], src->ref_dndt[1], src->ref_dndt[2],
            src->ref_tdb,
            src->alpha, src->rref);
#endif
    
    /* Astrometric place of source */
    source_place(&obs, src, &jtab, TR_MOTION | TR_PLX, s, dsdt, &pr);
    
    /* Take off actual gravitational deflection, and put on
       false deflection appropriate for star coordinates. */
    observer_ast2obs(&obs, s, dsdt, pr, TR_DEFL);
    observer_obs2ast(&obs, s, 0, TR_DEFL);
    
    v_to_ad(s, 0, &a, &d);
    v_to_ad_dt(s, dsdt, 0, &dadt, &dddt);
    
    a = ranormp(a);

    base10_to_60(a, UNIT_RAD, astr, sizeof(astr), " ", "", 3, UNIT_HR);
    base10_to_60(d, UNIT_RAD, dstr, sizeof(dstr), " ", "+", 2, UNIT_DEG);
    
    /* Observed place */
    observer_ast2obs(&obs, s, dsdt, 0, TR_TO_OBS_AZ);
    v_to_ad(s, 0, &az, &el);
    
    printf("%-16s %s %s %9.6f %9.6f %7.4f %5.1f\n",
           src->name, astr, dstr,
           dadt*RAD_TO_AS/DAY, dddt*RAD_TO_AS/DAY,
           1.0/pr, el*RAD_TO_DEG);
  }
  
  return(0);
}
