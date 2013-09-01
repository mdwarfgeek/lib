#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lfa.h"
#include "util.h"

static void cheby (double tc, double pfac, double vfac,
		   double *coef, int ncoef, int ncpt,
		   double *pos, double *vel);

int jpleph_open (struct jpleph_table *p, int type, char *filename) {
  char title[252];
  char *namebuf = (char *) NULL, *np;
  int32_t binvers, ncoeff, ncname, icon, ncon;
  int32_t ilc = -1;

#ifdef DEBUG
  int i;
#endif

  /* Get default file from environment if none specified */
  if(!filename) {
    if(type)
      filename = getenv("TIMEEPH_DATA");
    else
      filename = getenv("JPLEPH_DATA");
  }

  if(!filename)
    return(-2);

  /* Open ephemeris file */
  p->fp = fopen(filename, "rb");
  if(!(p->fp))
    return(-1);

  /* Read header */
  if(type == 1) {  /* Ephcom / TE405 format */
    fread(&binvers, sizeof(binvers), 1, p->fp);
    fread(&(p->denum), sizeof(p->denum), 1, p->fp);
    fread(&ncoeff, sizeof(ncoeff), 1, p->fp);
    fread(&(p->mjd_start), sizeof(p->mjd_start), 1, p->fp);
    fread(&(p->mjd_end), sizeof(p->mjd_end), 1, p->fp);
    fread(&(p->mjd_step), sizeof(p->mjd_step), 1, p->fp);
    fread(&(p->au), sizeof(p->au), 1, p->fp);
    fread(&(p->emratio), sizeof(p->emratio), 1, p->fp);
    fread(&(p->npt), sizeof(p->npt), 1, p->fp);

    if(ferror(p->fp)) {
      fclose(p->fp);
      return(-1);
    }

    if(feof(p->fp)) {
      fclose(p->fp);
      return(-3);  /* file damaged */
    }

    /* Make sure we know this version */
    if(binvers != 1) {
      fclose(p->fp);
      return(-3);  /* file damaged */
    }

    /* Allocate array for object pointers */
    p->ipt = (int32_t *) malloc(3*p->npt * sizeof(int32_t));
    if(!p->ipt) {
      fclose(p->fp);
      return(-1);
    }

    fread(p->ipt, sizeof(int32_t), 3*p->npt, p->fp);
    fread(&ncon, sizeof(ncon), 1, p->fp);
    fread(&ncname, sizeof(ncname), 1, p->fp);

    /* Allocate space for names */
    namebuf = (char *) malloc(ncname+1);
    if(!namebuf) {
      free((void *) p->ipt);
      fclose(p->fp);
      return(-1);
    }

    /* Read constant names and figure out useful ones */
    for(icon = 0; icon < ncon; icon++) {
      /* Read and nul-terminate */
      fread(namebuf, 1, ncname, p->fp);
      namebuf[ncname] = '\0';

      /* Strip off the junk */
      np = sstrip(namebuf);

      if(!strcmp(np, "LC"))
	ilc = icon;
    }
    
    free((void *) namebuf);

    /* Read title */
    fread(title, 1, sizeof(title), p->fp);

    /* Compute record size */
    p->recsize = ncoeff * sizeof(double);

    if(ilc >= 0) {
      /* Read constants */
      fseek(p->fp, p->recsize + ilc*sizeof(double), SEEK_SET);
      fread(&(p->lc), sizeof(double), 1, p->fp);
    }

    if(ferror(p->fp)) {
      free((void *) p->ipt);
      fclose(p->fp);
      return(-1);
    }
    
    if(feof(p->fp)) {
      free((void *) p->ipt);
      fclose(p->fp);
      return(-3);  /* file damaged */
    }
  }
  else {  /* original JPL format */
    p->npt = 13;

    p->ipt = (int32_t *) malloc(3*p->npt * sizeof(int32_t));
    if(!p->ipt) {
      fclose(p->fp);
      return(-1);
    }

    fread(title, 1, sizeof(title), p->fp);

    /* Skip constant names */
    fseek(p->fp, 2400, SEEK_CUR);

    fread(&(p->mjd_start), sizeof(p->mjd_start), 1, p->fp);
    fread(&(p->mjd_end), sizeof(p->mjd_end), 1, p->fp);
    fread(&(p->mjd_step), sizeof(p->mjd_step), 1, p->fp);
    fread(&ncon, sizeof(ncon), 1, p->fp);
    fread(&(p->au), sizeof(p->au), 1, p->fp);
    fread(&(p->emratio), sizeof(p->emratio), 1, p->fp);
    fread(p->ipt, sizeof(int32_t), 3*(p->npt-1), p->fp);
    fread(&(p->denum), sizeof(p->denum), 1, p->fp);
    fread(p->ipt + 3*(p->npt-1), sizeof(int32_t), 3, p->fp);

    if(ferror(p->fp)) {
      free((void *) p->ipt);
      fclose(p->fp);
      return(-1);
    }
    
    if(feof(p->fp)) {
      free((void *) p->ipt);
      fclose(p->fp);
      return(-3);  /* file damaged */
    }

    /* Check eph. type and set record length */
    switch(p->denum) {
    case 200:
      p->recsize = 6608;
      break;
    case 403:
    case 405:
    case 421:
      p->recsize = 8144;
      break;
    case 404:
    case 406:
      p->recsize = 5824;
      break;
    default:
      return(-2);
    }
  }    

  /* Convert dates to MJD */
  p->mjd_start -= ZMJD;
  p->mjd_end   -= ZMJD;

  /* Compute 1/(1+emratio), a commonly needed quantity */
  p->emfac = 1.0 / (1.0 + p->emratio);

  /* Allocate line buffer */
  p->buf = (double *) malloc(p->recsize);
  if(!p->buf) {
    free((void *) p->ipt);
    fclose(p->fp);
    return(-1);
  }

  p->brec = -1;

#ifdef DEBUG
  fprintf(stderr,
	  "JPL DE%d ephemeris opened, %d byte records\n"
	  "Title: %s\n"
	  "JD %.1lf-%.1lf @ %.1lf\n"
	  "%d constants: AU=%.3lf km\n"
	  "E/M=%.5lf\n",
	  p->denum, p->recsize,
	  title,
	  p->mjd_start, p->mjd_end, p->mjd_step,
	  ncon, p->au, p->emratio);

  for(i = 0; i < p->npt; i++)
    fprintf(stderr,
	    "ipt[%d] = { %d, %d, %d }\n",
	    i, p->ipt[i*3], p->ipt[i*3+1], p->ipt[i*3+2]);
#endif

  return(0);
}

void jpleph_close (struct jpleph_table *p) {
  free((void *) p->buf);
  p->buf = (double *) NULL;
  free((void *) p->ipt);
  p->ipt = (int32_t *) NULL;

  fclose(p->fp);
}

int jpleph_fetch (struct jpleph_table *p, double mjd, int body,
		  double *pos, double *vel) {
  int32_t iptr, ncoef, nsub, ncpt;

  size_t rec;
  int rv;

  double trec, tsub, tc, pfac, vfac;
  int irec, isub;

  double *tcoef;

  /* Range checks */
  if(mjd < p->mjd_start || mjd >= p->mjd_end)
    return(-2);

  if(body < 0 || body >= p->npt) {
    return(-2);
  }

  /* Fetch pointer information for this body */
  iptr = p->ipt[3*body] - 1;  /* convert to zero-indexed */
  ncoef = p->ipt[3*body+1];
  nsub = p->ipt[3*body+2];

  switch(body) {
  case JPLEPH_NUTATION:
    ncpt = 2;
    break;
  case TIMEEPH_TEI:
    ncpt = 1;
    break;
  default:
    ncpt = 3;
  }

  if(iptr < 0 || (iptr+ncpt*ncoef*nsub)*sizeof(double) > p->recsize)
    return(-3);  /* file damaged */

  /* Compute record number */
  trec = (mjd - p->mjd_start) / p->mjd_step;
  irec = floor(trec);

  /* Fetch record */
  if(p->brec < 0 || irec != p->brec) {
    /* Seek if needed */
    if(p->brec < 0 || irec != p->brec+1) {
      /* Compute offset in file */
      rec = (irec+2) * p->recsize;

      /* Seek to record */
      rv = fseek(p->fp, rec, SEEK_SET);
      if(rv)
	return(-1);
    }

    /* Read record */
    rv = fread(p->buf, p->recsize, 1, p->fp);
    if(rv != 1) {
      if(ferror(p->fp))
	return(-1);
      else
	return(-3);  /* file damaged */
    }

    p->brec = irec;
  }

  /* Compute subinterval */
  tsub = (trec - irec) * nsub;
  isub = floor(tsub);

  /* Chebyshev argument and scale factors */
  tc = 2.0 * (tsub - isub) - 1.0;

  if(body > JPLEPH_SUN &&
     body < TIMEEPH_TEV) {  /* angles or times, leave alone */
    pfac = 1.0;
    vfac = nsub / p->mjd_step;
  }
  else {  /* planets, convert to AU and AU/day */
    pfac = 1.0 / p->au;
    vfac = nsub / (p->mjd_step * p->au);
  }

  /* Evaluate */
  tcoef = p->buf + iptr + ncpt*ncoef*isub;
  cheby(tc, pfac, vfac, tcoef, ncoef, ncpt, pos, vel);

  return(0);
}

static void cheby (double tc, double pfac, double vfac,
		   double *coef, int ncoef, int ncpt,
		   double *pos, double *vel) {
  double pc[ncoef], vc[ncoef];
  double twot = 2*tc;

  int i, p, op;

  /* Evaluate polynomial and derivative */
  pc[0] = 1.0;
  pc[1] = tc;

  vc[0] = 0.0;
  vc[1] = 1.0;

  for(p = 2; p < ncoef; p++) {
    pc[p] = twot * pc[p-1] - pc[p-2];
    vc[p] = twot * vc[p-1] + 2*pc[p-1] - vc[p-2];
  }

  for(i = 0; i < ncpt; i++) {
    op = i*ncoef;

    /* Sum backwards */
    pos[i] = 0.0;
    vel[i] = 0.0;
    
    for(p = ncoef-1; p >= 0; p--) {
      pos[i] += pc[p]*coef[op+p];
      vel[i] += vc[p]*coef[op+p];
    }
    
    pos[i] *= pfac;
    vel[i] *= 2.0*vfac;
  }
}
