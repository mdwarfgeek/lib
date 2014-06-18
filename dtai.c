#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>

#include "lfa.h"
#include "util.h"

int dtai_read (struct dtai_table *tab, char *filename) {
  FILE *fp;
  char fnbuf[1024], buf[1024], *p;
  int rv;

  double jd, val, mjdzero, scale;

  tab->table = (struct dtai_entry *) NULL;
  tab->ilast = 0;
  tab->ntab = 0;

  /* Get default file from environment if none specified */
  if(!filename) {
    filename = getenv("IERS_DATA");
    if(!filename)
      return(-2);

    snprintf(fnbuf, sizeof(fnbuf), "%s/tai-utc.dat", filename);
    filename = fnbuf;
  }

  /* Open file */
  fp = fopen(filename, "r");
  if(!fp)
    return(-1);

  /* Read sorted entries */
  while(fgets(buf, sizeof(buf), fp)) {
    p = sstrip(buf);

    /* This format is not particularly machine friendly.  When parsing,
       we should really try to guard against anything they might change. */

    /* First skip the Gregorian date and try to find "JD" */
    p = strstr(p, "JD");
    if(!p)
      goto error;

    /* Simple approach for now */
    rv = sscanf(p, "JD %lf TAI-UTC= %lf S + (MJD - %lf) X %lf",
		&jd, &val, &mjdzero, &scale);
    if(rv != 4)
      goto error;

    /* Allocate new entry */
    tab->table = (struct dtai_entry *) realloc(tab->table,
					       (tab->ntab+1) * sizeof(struct dtai_entry));
    if(!tab->table) {
      tab->ntab = 0;
      fclose(fp);
      return(-1);
    }
    
    tab->table[tab->ntab].mjd = rint(jd-ZMJD);
    tab->table[tab->ntab].dtai = val;
    tab->table[tab->ntab].mjdzero = mjdzero;
    tab->table[tab->ntab].scale = scale;
    tab->ntab++;
  }

  if(ferror(fp)) {
    fclose(fp);
    return(-1);
  }

  /* Start off at the end, this will usually be right */
  tab->ilast = tab->ntab;

  fclose(fp);

  return(0);

 error:
  fclose(fp);

  return(-3);
}

void dtai_free (struct dtai_table *tab) {
  /* Clear table, if there is one */
  if(tab->ntab > 0) {
    free((void *) tab->table);
    tab->table = (struct dtai_entry *) NULL;
    tab->ilast = 0;
    tab->ntab = 0;
  }
}

/* This routine returns TAI-UTC (in seconds) for a given UTC.
   UTC is supplied as integer MJD number and fraction of day in
   order to avoid problems during leap seconds.  During a
   positive leap second, the MJD number should be held back
   until the first second of the next day. */

double dtai (struct dtai_table *tab, int mjd, double frac) {
  int i;

  /* Check if we're already there */
  i = tab->ilast;

  /* Step backwards */
  while(i > 1 && mjd < tab->table[i-1].mjd)
    i--;

  /* Step forwards */
  while(i < tab->ntab && mjd >= tab->table[i].mjd)
    i++;
  
  /* Store pointer */
  tab->ilast = i;

  /* Compute answer */
  return(tab->table[i-1].dtai +
	 (frac + (mjd - tab->table[i-1].mjdzero))*tab->table[i-1].scale);
}


