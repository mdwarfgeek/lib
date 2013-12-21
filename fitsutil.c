#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <math.h>

#include "lfa.h"
#include "fitsutil.h"

#include "cvtunit.h"
#include "util.h"

void fitsio_err (char *errstr, int status, const char *fmt, ...) {
  char errmsg[FLEN_STATUS];
  va_list ap;
  int rv;

  ffgerr(status, errmsg);
  
  va_start(ap, fmt);
  rv = vsnprintf(errstr, ERRSTR_LEN, fmt, ap);
  va_end(ap);
  
  if(rv != -1 && rv < (ERRSTR_LEN - 1))
    (void) snprintf(errstr + rv, ERRSTR_LEN - rv, ": %s", errmsg);
}

int read_wcs (fitsfile *fits, struct wcs_info *wcs, int verbose, char *errstr) {
  int status = 0, havecd, i;
  char ctype1[FLEN_VALUE], ctype2[FLEN_VALUE], key[FLEN_KEYWORD];
  double v, scl1, scl2, rho, sinrho, cosrho, atmp, btmp, dtmp, etmp;
 
  /* Figure out WCS type */
  wcs->mpc = -1;

  ffgkys(fits, "CTYPE1", ctype1, (char *) NULL, &status);
  ffgkys(fits, "CTYPE2", ctype2, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    if(verbose)
      printf("No WCS type keywords found, assuming TAN\n");
    wcs->proj = PROJ_TAN;
  }
  else if(status) {
    fitsio_err(errstr, status, "could not read WCS CTYPE keywords");
    goto error;
  }
  else {
    if(!strncmp(ctype1, "RA---TAN", 8) &&
       !strncmp(ctype2, "DEC--TAN", 8)) {
      if(verbose)
	printf("Found TAN WCS\n");
      wcs->proj = PROJ_TAN;
    }
    else if(!strncmp(ctype1, "RA---SIN", 8) &&
	    !strncmp(ctype2, "DEC--SIN", 8)) {
      if(verbose)
	printf("Found SIN WCS\n");
      wcs->proj = PROJ_SIN;
    }
    else if(!strncmp(ctype1, "RA---ARC", 8) &&
	    !strncmp(ctype2, "DEC--ARC", 8)) {
      if(verbose)
	printf("Found ARC WCS\n");
      wcs->proj = PROJ_ARC;
    }
    else if(!strncmp(ctype1, "RA---ZPN", 8) &&
	    !strncmp(ctype2, "DEC--ZPN", 8)) {
      if(verbose)
	printf("Found ZPN WCS\n");
      wcs->proj = PROJ_ARC;

      /* Read projection constants from header.  The default if there
	 are none is an ARC projection.  Technically, the FITS WCS
	 standard requires that all constants default to zero, but
	 the author believes the present default behaviour is more
	 sensible and useful (all constants zero would crash when
	 trying to invert the projection in xy2vec). */

      for(i = 0; i <= PC_MAX; i++) {
	/* Try WCS standard (2002) form */
	snprintf(key, sizeof(key), "PV2_%d", i);
	ffgkyd(fits, key, &v, (char *) NULL, &status);
	if(status == KEY_NO_EXIST) {
	  status = 0;

	  /* Try old form */
	  snprintf(key, sizeof(key), "PROJP%d", i);
	  ffgkyd(fits, key, &v, (char *) NULL, &status);
	  if(status == KEY_NO_EXIST) {
	    status = 0;
	    v = 0;
	  }
	  else if(status) {
	    fitsio_err(errstr, status, "ffgkyd: %s", key);
	    goto error;
	  }
	  else {
	    if(i == 0 && v != 0) {
	      report_err(errstr, "PROJP0 != 0 is not supported");
	      goto error;
	    }

	    wcs->mpc = i;
	  }
	}
	else if(status) {
	  fitsio_err(errstr, status, "ffgkyd: %s", key);
	  goto error;
	}
	else {
	  if(i == 0 && v != 0) {
	    report_err(errstr, "PV2_0 != 0 is not supported");
	    goto error;
	  }

	  wcs->mpc = i;
	}

	wcs->pc[i] = v;
      }
    }
    else {
      report_err(errstr, "unrecognised WCS type: %s %s", ctype1, ctype2);
      goto error;
    }
  }

  /* Tangent point */
  ffgkyd(fits, "CRPIX1", &(wcs->c), (char *) NULL, &status);
  ffgkyd(fits, "CRPIX2", &(wcs->f), (char *) NULL, &status);
  ffgkyd(fits, "CRVAL1", &(wcs->tpa), (char *) NULL, &status);
  ffgkyd(fits, "CRVAL2", &(wcs->tpd), (char *) NULL, &status);

  /* Try for sensible CD matrix form first */
  havecd = 4;

  ffgkyd(fits, "CD1_1", &(wcs->a), (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    wcs->a = 0.0;
    havecd--;
  }
  ffgkyd(fits, "CD1_2", &(wcs->b), (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    wcs->b = 0.0;
    havecd--;
  }
  ffgkyd(fits, "CD2_1", &(wcs->d), (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    wcs->d = 0.0;
    havecd--;
  }
  ffgkyd(fits, "CD2_2", &(wcs->e), (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    wcs->e = 0.0;
    havecd--;
  }

  if(havecd == 0) {
    /* No CD matrix.  More work then.  First try to read PC matrix.
     * Note default is the identity.
     */
    wcs->a = 1.0;
    wcs->b = 0.0;
    wcs->d = 0.0;
    wcs->e = 1.0;

    ffgkyd(fits, "PC1_1", &(wcs->a), (char *) NULL, &status);
    if(status == KEY_NO_EXIST)
      status = 0;
    ffgkyd(fits, "PC1_2", &(wcs->b), (char *) NULL, &status);
    if(status == KEY_NO_EXIST)
      status = 0;
    ffgkyd(fits, "PC2_1", &(wcs->d), (char *) NULL, &status);
    if(status == KEY_NO_EXIST)
      status = 0;
    ffgkyd(fits, "PC2_2", &(wcs->e), (char *) NULL, &status);
    if(status == KEY_NO_EXIST)
      status = 0;
  
    /* Now read and multiply in CDELTs */
    ffgkyd(fits, "CDELT1", &scl1, (char *) NULL, &status);
    ffgkyd(fits, "CDELT2", &scl2, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      scl1 = 1.0;
      scl2 = 1.0;
    }
    else if(status) {
      report_err(errstr, "ffgkyd: CDELT1,2");
      goto error;
    }
    
    wcs->a *= scl1;
    wcs->b *= scl1;
    wcs->d *= scl2;
    wcs->e *= scl2;

    /* Now read and apply the dreaded (and ill-defined) CROTAs */
    ffgkyd(fits, "CROTA1", &rho, (char *) NULL, &status);
    if(status == KEY_NO_EXIST)
      status = 0;
    else if(status) {
      fitsio_err(errstr, status, "ffgkyd: CROTA1");
      goto error;
    }
    else {
      /* Apply rotation matrix with angle rho (in degrees!) */
      dsincos(rho*DEG_TO_RAD, &sinrho, &cosrho);

      atmp = wcs->a * cosrho - wcs->d * sinrho;
      btmp = wcs->b * cosrho - wcs->e * sinrho;
      dtmp = wcs->a * sinrho + wcs->d * cosrho;
      etmp = wcs->b * sinrho + wcs->e * cosrho;

      wcs->a = atmp;
      wcs->b = btmp;
      wcs->d = dtmp;
      wcs->e = etmp;
    }

    ffgkyd(fits, "CROTA2", &rho, (char *) NULL, &status);
    if(status == KEY_NO_EXIST)
      status = 0;
    else if(status) {
      fitsio_err(errstr, status, "ffgkyd: CROTA2");
      goto error;
    }
    else {
      /* Apply rotation matrix with angle rho (in degrees!) */
      dsincos(rho*DEG_TO_RAD, &sinrho, &cosrho);

      atmp = wcs->a * cosrho + wcs->d * sinrho;
      btmp = wcs->b * cosrho + wcs->e * sinrho;
      dtmp = -wcs->a * sinrho + wcs->d * cosrho;
      etmp = -wcs->b * sinrho + wcs->e * cosrho;

      wcs->a = atmp;
      wcs->b = btmp;
      wcs->d = dtmp;
      wcs->e = etmp;
    }
  }
  
  if(status) {
    fitsio_err(errstr, status, "could not read WCS keywords");
    goto error;
  }

  wcs->tpa *= DEG_TO_RAD;
  wcs->tpd *= DEG_TO_RAD;
  dsincos(wcs->tpa, &(wcs->sina), &(wcs->cosa));
  dsincos(wcs->tpd, &(wcs->sind), &(wcs->cosd));

  wcs->a *= DEG_TO_RAD;
  wcs->b *= DEG_TO_RAD;
  wcs->d *= DEG_TO_RAD;
  wcs->e *= DEG_TO_RAD;
    
/*    if(verbose)
      printf("\nFrame transform constants:\n"
	     "%.5e %.5e %7.2f\n%.5e %.5e %7.2f\n"
	     "Tangent points: %.5f %.5f\n",
	     wcs->a, wcs->b, wcs->c, wcs->d, wcs->e, wcs->f,
	     wcs->tpa, wcs->tpd); */

  return(0);

 error:

  return(1);
}

