#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>

#include "cvtunit.h"
#include "lfa.h"
#include "util.h"

int mpc_read (char *filename, struct source **srclist, int *nsrc) {
  FILE *fp;
  char buf[1024], *p;
  int blank, nalloc;

  /* Open file */
  fp = fopen(filename, "r");
  if(!fp)
    return(-1);

  /* Read entries */
  while(fgets(buf, sizeof(buf), fp)) {
    /* Skip blank lines */
    blank = 1;
    for(p = &(buf[0]); *p; p++)
      if(!isspace((unsigned char) *p)) {
	blank = 0;
	break;
      }

    if(blank)
      continue;

    /* Allocate space for output */
    nalloc = (*nsrc)+1;

    *srclist = (struct source *) realloc(*srclist,
					 nalloc * sizeof(struct source));
    if(!(*srclist)) {
      fclose(fp);
      return(-1);
    }

    /* Convert to internal format (universal elements) */
    if(mpc_convert(buf, (*srclist) + (*nsrc)))
      goto error;

    (*nsrc)++;
  }

  if(ferror(fp)) {
    fclose(fp);
    return(-1);
  }

  fclose(fp);

  return(0);

 error:
  fclose(fp);

  return(-3);
}

static int packed_cent (char c) {
  if(islower((unsigned char) c))
    return(c - 'i' + 18);
  if(isupper((unsigned char) c))
    return(c - 'I' + 18);
  else
    return(-1);
}

static int packed_md (char c) {
  if(isdigit((unsigned char) c))
    return(c - '0');
  if(islower((unsigned char) c))
    return(c - 'a' + 10);
  if(isupper((unsigned char) c))
    return(c - 'A' + 10);
  else
    return(-1);
}

int mpc_convert (char *line, struct source *src) {
  int len;
  char *p;

  int yr, mn, di;
  double df, epoch;

  double ma, aq, ecc, perih, anode, incl, nn;

  len = strlen(line);

  /* We need to figure out if this is a minor planet or comet entry.
     Column 22 should be part of the (packed) epoch, which is always
     present, for a minor planet, but blank for a comet.  Numbering is
     from 1, so we test array index 21, after making sure there are
     enough characters. */
  if(len < 22)
    goto error;
  
  if(isspace((unsigned char) line[21])) {  /* comet */
    /* Epoch */
    if(extractint(line, len, 15, 18, &yr))
      goto error;
    if(extractint(line, len, 20, 21, &mn))
      goto error;
    if(extractintfrac(line, len, 23, 29, &di, &df))
      goto error;
    
    epoch = date2mjd(yr, mn, di) + df;  /* convert to MJD(TT) */
    
    /* Elements */
    if(extractdouble(line, len, 31, 39, &aq))  /* perihelion distance, AU */
      goto error;
    if(extractdouble(line, len, 42, 49, &ecc))  /* eccentricity */
      goto error;
    
    if(extractdouble(line, len, 52, 59, &perih))  /* argument of perihelion, deg */
      goto error;
    if(extractdouble(line, len, 62, 69, &anode))  /* longitude of ascending node, deg */
      goto error;
    if(extractdouble(line, len, 72, 79, &incl))  /* inclination, deg */
      goto error;
    
    /* Form source structure */
    source_elem(src, SOURCE_ELEM_COMET,
		epoch,
		incl * DEG_TO_RAD,
		anode * DEG_TO_RAD,
		perih * DEG_TO_RAD,
		aq,
		ecc,
		0.0,
		0.0);
    
    /* Designation */
    p = extractstr(line, len, 103, 158);
    if(p) {
      strncpy(src->name, p, sizeof(src->name)-1);
      src->name[sizeof(src->name)-1] = '\0';
    }
    else
      src->name[0] = '\0';

#ifdef TEST
    fprintf(stderr, "Got comet '%s': epoch=%lf i=%lf anode=%lf perih=%lf q=%lf e=%lf\n",
	    src->name, epoch, incl, anode, perih, aq, ecc);
#endif
  }
  else {  /* minor planet */
    /* Packed epoch */
    p = extractstr(line, len, 21, 25);
    
    yr = packed_cent(p[0]);
    if(yr <= 0)
      goto error;
    
    yr = (yr * 10 + (p[1] - '0')) * 10 + (p[2] - '0');
    
    mn = packed_md(p[3]);
    if(mn < 1 || mn > 12)
      goto error;
    
    di = packed_md(p[4]);
    if(di < 1 || di > 31)
      goto error;

    epoch = date2mjd(yr, mn, di) - 0.5;  /* convert to MJD(TT) */
    
    /* Elements */
    if(extractdouble(line, len, 27, 35, &ma))  /* mean anomaly at epoch, deg */
      goto error;
    if(extractdouble(line, len, 38, 46, &perih))  /* argument of perihelion, deg */
      goto error;
    if(extractdouble(line, len, 49, 57, &anode))  /* longitude of ascending node, deg */
      goto error;
    if(extractdouble(line, len, 60, 68, &incl))  /* inclination, deg */
      goto error;
    if(extractdouble(line, len, 71, 79, &ecc))  /* eccentricity */
      goto error;
    if(extractdouble(line, len, 81, 91, &nn))  /* mean motion, deg/day */
      goto error;
    if(extractdouble(line, len, 93, 103, &aq))  /* semimajor axis, AU */
      goto error;
    
    /* Form source structure */
    source_elem(src, SOURCE_ELEM_MINOR,
		epoch,
		incl * DEG_TO_RAD,
		anode * DEG_TO_RAD,
		perih * DEG_TO_RAD,
		aq,
		ecc,
		ma * DEG_TO_RAD,
		nn * DEG_TO_RAD);

    /* Designation */
    p = extractstr(line, len, 167, 194);
    if(p) {
      strncpy(src->name, p, sizeof(src->name)-1);
      src->name[sizeof(src->name)-1] = '\0';
    }
    else
      src->name[0] = '\0';

#ifdef TEST
    fprintf(stderr, "Got MP '%s': epoch=%lf i=%lf anode=%lf perih=%lf a=%lf e=%lf ma=%lf nn=%lf\n",
	    src->name, epoch, incl, anode, perih, aq, ecc, ma, nn);
#endif
  }

  return(0);

 error:

  return(-3);
}
