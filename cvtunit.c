#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>

#include "cvtunit.h"

/* Array of conversion factors.  These convert from degrees, to each unit. */
static double cfactor[] = {
  M_PI / 180.0,  /* Rad */
  1.0 / 15.0,    /* Hours */
  4.0,           /* Minutes */
  240.0,         /* Seconds */
  14400.0,       /* Milliseconds */
  1.0,           /* Degrees */
  60.0,          /* Arcmin */
  3600.0,        /* Arcsec */
  3600000.0,     /* MAS */
};

int base60_to_10 (char *inangle, char **endptr, char *sep, short inunit,
		  double *outangle, short outunit) {
  double td, inangle_d = 0.0, sign = +1.0;
  long tl;
  char *str = (char *) NULL, *tokens[3], *p, *ep;
  size_t seplen;
  int i = 0, ilim;
  short sepspace = 0;

  /* Check the unit parameters */
  if(inunit < 0 || inunit > UNIT_MAS || outunit < 0 || outunit > UNIT_MAS)
    goto error;

  if(!sep)
    sep = ":";

  /* Make a copy of the string to work on. */
  str = strdup(inangle);
  if(!str) {
    fprintf(stderr, "malloc: %s", strerror(errno));
    exit(1);
  }

  /* First attempt to convert to base 10, still in 'inunit'.  We try three
   * possible formats, with the given separator 'sep'.  Note: be sensible
   * with the 'sep' parameter, things will generally go wrong if it is
   * numeric (decimals, '.', etc)!
   */

  /* Skip any leading whitespace */
  p = str;

  while(*p && isspace((unsigned char) *p))
    p++;

  if(*p == '\0')
    goto error;

  /* Remove any trailing whitespace */
  /* XXX - todo */

  /* Check for a sign at the start of the string ('+' or '-'). */
  if(*p == '+')
    p++;
  else if(*p == '-') {
    sign = -1.0;
    p++;
  }

  /* Split the string into upto 3 tokens, on 'sep'.  If 'sep' is a single
   * space character, accept any number of space characters.
   */
  seplen = strlen(sep);

  if(seplen == 1 && *sep == ' ')
    sepspace = 1;

  while(i < 3) {
    tokens[i] = p;

    ep = strstr(p, sep);
    if(!ep)
      break;

    *ep = '\0';
    
    i++;

    p = ep + seplen;

    if(sepspace)
      while(*p == ' ')
	p++;

    if(!*p)  /* end of string */
      break;
  }

  if(i == 0)
    goto error;
  if(i == 3)
    i--;

  /* Check and parse the tokens as numbers.  The last one should be
   * parsed as a double, the rest as longs.  Signs are not allowed
   * within tokens.
   */
  ilim = i;

  for(i = 0; i < ilim; i++) {
    p = tokens[i];

    if(!isdigit((unsigned char) *p))
      goto error;

    tl = strtol(p, &ep, 10);
    if(*ep != '\0')
      goto error;

    inangle_d += (double) tl;
    inangle_d *= 60;
  }

  p = tokens[ilim];

  if(!isdigit((unsigned char) *p) && *p != '.')
    goto error;

  td = strtod(p, &ep);
  if(endptr)
    *endptr = inangle + (int) (ep - str);
  else if(*ep != '\0')
    goto error;

  inangle_d += td;

  /* Apply the sign and convert to units of 'inunit' */
  inangle_d *= sign / pow(60.0, (double) ilim);

  /* The variable inangle_d now contains the value in units of 'inunit'.
   * convert to the relevant output unit 'outunit'. */
  *outangle = inangle_d * (cfactor[outunit] / cfactor[inunit]);

  /* Done, free the temporary string */
  free((void *) str);
  str = (char *) NULL;

  return(0);

 error:
  if(str) {
    free((void *) str);
    str = (char *) NULL;
  }

  return(1);
}

int base10_to_60 (double inangle, short inunit, char *outangle, size_t outlen,
		  char *sep, char *sign, short dp, short outunit) {
  double outangle_d, round, min, sec;
  int rv = 0;

  /* Check the unit parameters */
  if(inunit < 0 || inunit > UNIT_MAS || outunit < 0 || outunit > UNIT_MAS)
    goto error;

  if(outlen < 1)
    goto error;

  if(!sep)
    sep = ":";
  if(!sign)
    sign = "";
  if(dp == -1)
    dp = 2;

  /* Extract the sign and convert to units of 'outunit'. */
  if(inangle < 0) {
    sign     = "-";
    inangle *= -1;
  }

  outangle_d = inangle * (cfactor[outunit] / cfactor[inunit]);

  round = 0.5 / pow(10, (double) dp);
  outangle_d += round / 3600;

  /* Split into hours, minutes, seconds */
  min = (outangle_d - ((int) outangle_d)) * 60.0;
  sec = (min        - ((int) min))        * 60.0;

  /* Write out the angle. */
  rv = snprintf(outangle, outlen, "%s%02d%s%02d%s%0*.*f", sign,
		(int) outangle_d, sep, (int) min, sep,
		(dp > 0) ? dp + 3 : 2, dp, fabs(sec - round));
  if(rv == -1 || rv >= outlen)
    goto error;

  return(rv);

 error:
  return(-1);
}

double get_cfactor (short inunit, short outunit) {
  /* Check the unit parameters */
  if(inunit < 0 || inunit > UNIT_MAS || outunit < 0 || outunit > UNIT_MAS)
    return(0.0);  /* XXX - need a better way to handle this */

  /* Return conversion factor */
  return(cfactor[outunit] / cfactor[inunit]);
}

int parsedate (char *indate, char *sep, int *yr, int *mn, int *dy) {
  long tl[3];
  char *str = (char *) NULL, *tokens[3], *p, *ep;
  size_t seplen;
  int i = 0, ilim;

  /* Check the unit parameters */
  if(!sep)
    sep = "-";

  /* Make a copy of the string to work on. */
  str = strdup(indate);
  if(!str) {
    fprintf(stderr, "malloc: %s", strerror(errno));
    exit(1);
  }

  /* Skip any leading whitespace */
  p = str;

  while(*p && isspace((unsigned char) *p))
    p++;

  if(*p == '\0')
    goto error;

  /* Remove any trailing whitespace */
  /* XXX - todo */

  /* Split the string into upto 3 tokens, on 'sep'. */
  seplen = strlen(sep);

  while(i < 3) {
    tokens[i] = p;

    ep = strstr(p, sep);
    if(!ep)
      break;

    *ep = '\0';
    
    i++;

    p = ep + seplen;
    if(!*p)  /* end of string */
      break;
  }

  /* Make sure we have exactly 3 tokens */
  if(i != 2)
    goto error;

  /* Check and parse the tokens as unsigned numbers. */
  ilim = i;

  for(i = 0; i <= ilim; i++) {
    p = tokens[i];

    if(!isdigit((unsigned char) *p))
      goto error;

    tl[i] = strtol(p, &ep, 10);
    if(*ep != '\0')
      goto error;
  }

  /* Copy out to the user. */
  if(yr)
    *yr = (int) tl[0];
  if(mn)
    *mn = (int) tl[1];
  if(dy)
    *dy = (int) tl[2];

  /* Done, free the temporary string */
  free((void *) str);
  str = (char *) NULL;

  return(0);

 error:
  if(str) {
    free((void *) str);
    str = (char *) NULL;
  }

  return(1);
}
