#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>

#include "lfa.h"

int iers_open (struct iers_table *tab,
	       struct dtai_table *dtai_tab,
	       char *filename) {
  char fnbuf[BUFSIZE], buf[BUFSIZE];
  double mjd;
  int rv;
  long off;

  /* Get default file from environment if none specified */
  if(!filename) {
    filename = getenv("IERS_DATA");
    if(!filename)
      return(-2);

    snprintf(fnbuf, sizeof(fnbuf), "%s/finals2000A.data", filename);
    filename = fnbuf;
  }

  /* Open file.  Note that on Windows we have to do this in binary
     mode to disable CRLF translation.  Otherwise, the CR is stripped
     off and defeats our efforts to determine (and check) the record
     length.  The extra CRs are not a problem. */
  tab->fp = fopen(filename, "rb");
  if(!(tab->fp)) {
    return(-1);
  }  

  /* Read first record to figure out record length and start date */
  if(!fgets(buf, sizeof(buf), tab->fp)) {
    fclose(tab->fp);
    return(-1);
  }

  tab->recsize = strlen(buf);

  if(extractdouble(buf, tab->recsize, 8, 15, &mjd)) {
    fclose(tab->fp);
    return(-2);
  }

  tab->mjd_start = mjd;

  /* Read second record to figure out step */
  if(!fgets(buf, sizeof(buf), tab->fp)) {
    fclose(tab->fp);
    return(-1);
  }

  /* Sanity check that record lengths are the same */
  if(strlen(buf) != tab->recsize) {
    fclose(tab->fp);
    return(-2);
  }

  if(extractdouble(buf, tab->recsize, 8, 15, &mjd)) {
    fclose(tab->fp);
    return(-2);
  }

  tab->mjd_step = mjd - tab->mjd_start;

  if(tab->mjd_step <= 0) {
    fclose(tab->fp);
    return(-2);
  }

  /* Seek to last record */
  rv = fseek(tab->fp, -tab->recsize, SEEK_END);
  if(rv < 0) {
    fclose(tab->fp);
    return(-1);
  }

  /* Where are we? */
  off = ftell(tab->fp);
  if(off < 0) {
    fclose(tab->fp);
    return(-1);
  }

  tab->lrec = off / tab->recsize;

  /* Adjust file pointer as necessary in case there is truncation
     during the last record.  This is a non-critical problem. */
  off -= ((long) tab->lrec) * tab->recsize;

  if(off > 0) {
    rv = fseek(tab->fp, -off, SEEK_CUR);
    if(rv < 0) {
      fclose(tab->fp);
      return(-1);
    }
  }

  /* Try to discover where there are useful measurements. */
  for(;;) {
    /* Read last record to check range */
    if(!fgets(buf, sizeof(buf), tab->fp)) {
      fclose(tab->fp);
      return(-1);
    }
    
    /* Sanity check that record lengths are the same */
    if(strlen(buf) != tab->recsize) {
      fclose(tab->fp);
      return(-2);
    }

    /* Column 17 blank = no prediction, so wait until it's not */
    if(!isspace(buf[16]))
      break;

    /* Seek back, if we can */
    tab->lrec--;
    if(tab->lrec <= 0)
      break;

    rv = fseek(tab->fp, -2*tab->recsize, SEEK_CUR);
    if(rv < 0) {
      fclose(tab->fp);
      return(-1);
    }
  }

  /* We should now be pointing at the last record with data. */
  if(extractdouble(buf, tab->recsize, 8, 15, &mjd)) {
    fclose(tab->fp);
    return(-2);
  }

  /* Sanity check MJD is what we think it should be.  Dates in
     files are given to 2 dp, so we require agreement to that
     level. */
  if(fabs(tab->mjd_start+tab->mjd_step*tab->lrec - mjd) > 0.01) {
    fclose(tab->fp);
    return(-2);
  }

  /* Init a few things */
  tab->bufrec = -1;
  tab->filerec = -1;
  tab->dtai_tab = dtai_tab;

  return(0);
}

void iers_close (struct iers_table *tab) {
  fclose(tab->fp);
}

static int iers_read (struct iers_table *tab,
		      struct iers_entry *entbuf,
		      int frec, int nrec) {
  struct iers_entry *pent;

  char buf[BUFSIZE];
  int irec, rv;

  double ttmutc;

  int imjd;
  double fmjd;

  /* Seek to requested record, if needed */
  if(tab->filerec != frec) {
    rv = fseek(tab->fp, frec * tab->recsize, SEEK_SET);
    if(rv)
      return(-1);

    tab->filerec = frec;
  }

  /* Read records */
  for(irec = 0, pent = entbuf; irec < nrec; irec++, pent++) {
    /* Read raw record */
    rv = fread(buf, 1, tab->recsize, tab->fp);
    if(rv != tab->recsize)
      return(-1);

    /* Make sure properly null-terminated */
    buf[tab->recsize] = '\0';

    /* Increment record counter */
    tab->filerec++;

    /* MJD: file is broken if this fails */
    if(extractintfrac(buf, tab->recsize, 8, 15, &imjd, &fmjd))
      return(-2);

    /* These three available for full predicted period */
    if(extractdouble(buf, tab->recsize, 19, 27, &(pent->xp)))
      return(-2);
    if(extractdouble(buf, tab->recsize, 38, 46, &(pent->yp)))
      return(-2);
    if(extractdouble(buf, tab->recsize, 59, 68, &(pent->dut1)))
      return(-2);

    /* These two not always available for predicted */
    if(extractdouble(buf, tab->recsize, 98, 106, &(pent->dxnut)))
      pent->dxnut = 0;  /* should really handle these better */
    if(extractdouble(buf, tab->recsize, 117, 125, &(pent->dynut)))
      pent->dynut = 0;

    /* The file gives UT1-UTC, whereas we want TT-UT1 (the standard
       "Delta T" quantity) for interpolation purposes, because UT1-UTC
       has discontinuities at leap seconds.  Put in the appropriate
       offset. */
    ttmutc = dtai(tab->dtai_tab, imjd, fmjd)+DTT;
    pent->dut1 = ttmutc - pent->dut1;

    /* Combine separate parts of MJD, now they have served their purpose */
    pent->mjd = imjd + fmjd;
  }

  return(0);
}

/* Routine to fetch and interpolate TT-UT1 for a given UTC (as MJD).
   TT is used in the output to avoid discontinuities due to leap
   seconds.  While technically the input argument should be UTC, it
   doesn't really matter what is used for most applications.
   Running off the ends of the table is allowed, and will cause
   the routine to return 1 rather than the usual value of 0. */

int iers_fetch (struct iers_table *tab, double mjd,
		double *dut1, double *xp, double *yp, double *dxnut, double *dynut) {
  struct iers_entry *p;
  int trec, clamp = 0, f, n, rv;
  double a = 0, b = 0;

  /* Record before target.  Clamped to be within file range. */
  trec = floor((mjd - tab->mjd_start) / tab->mjd_step);
  if(trec < 0) {
    trec = 0;
    a = 1;
    clamp = 1;
  }
  else if(trec > (tab->lrec-1)) {
    trec = tab->lrec-1;
    b = 1;
    clamp = 1;
  }

  /* Do we currently have it in memory? */
  if(trec != tab->bufrec) {  /* nope */
    p = &(tab->buf[0]);
    f = trec;
    n = 2;

    /* Figure out what we can reuse, if anything */
    if(tab->bufrec >= 0) {
      if(trec == tab->bufrec + 1) {
	memmove(&(tab->buf[0]), &(tab->buf[1]), sizeof(tab->buf[0]));
	p = &(tab->buf[1]);
	f = trec+1;
	n = 1;
      }
      else if(trec == tab->bufrec - 1) {
	memmove(&(tab->buf[1]), &(tab->buf[0]), sizeof(tab->buf[0]));
	n = 1;
      }
    }

    /* Read */
    rv = iers_read(tab, p, f, n);
    if(rv != 0)
      return(rv);

    /* Precompute local date step (in case different) */
    tab->h = tab->buf[1].mjd - tab->buf[0].mjd;
    if(tab->h == 0.0)
      return(-2);  /* bad table */
  }

  /* Interpolate, linear only for now */
  if(!clamp) {
    a = (tab->buf[1].mjd - mjd) / tab->h;
    b = (mjd - tab->buf[0].mjd) / tab->h;
  }

  if(dut1)
    *dut1  = (a * tab->buf[0].dut1  +  b * tab->buf[1].dut1);
  if(xp)
    *xp    = (a * tab->buf[0].xp    +  b * tab->buf[1].xp);
  if(yp)
    *yp    = (a * tab->buf[0].yp    +  b * tab->buf[1].yp);
  if(dxnut)
    *dxnut = (a * tab->buf[0].dxnut +  b * tab->buf[1].dxnut);
  if(dxnut)
    *dynut = (a * tab->buf[0].dynut +  b * tab->buf[1].dynut);

  return(clamp);  /* let user know if we ran off end */
}

