#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>

#include "lfa.h"

char *extractstr (char *str, int len,
		  int f, int l) {
  static char tmp[BUFSIZE];
  int i, n;

  /* Check args */
  if(f < 1 || f >= l || l > len)
    return((char *) NULL);

  /* Compute length excluding terminating NUL */
  n = l - f + 1;

  /* Copy */
  for(i = 0; i < n; i++)
    tmp[i] = str[f+i-1];

  tmp[n] = '\0';

  /* Strip off white space when returning */
  return(sstrip(tmp));
}

int extractdouble (char *str, int len,
		   int f, int l,
		   double *val) {
  static char *p, *ep;

  p = extractstr(str, len, f, l);
  if(!p)
    return(-1);

  /* Now convert */
  *val = strtod(p, &ep);
  if(ep == p)
    return(-1);

  return(0);
}

/* A simple routine to read integer and fractional parts of a decimal
   number separately, such that they give the original value when added
   together (i.e. both parts will be negative if the input was negative).
   Used to read MJDs from the file.  Does not handle scientific notation. */

int extractintfrac (char *str, int len,
		    int f, int l,
		    int *ival, double *fval) {
  static char *p, *dp, *ep;
  int v, ndig;

  p = extractstr(str, len, f, l);
  if(!p)
    return(-1);

  /* Locate the decimal point (if there is one) */
  dp = strchr(p, '.');
  if(dp)
    *dp = '\0';

  /* Convert integer part */
  if(!dp || dp > p) {
    *ival = strtol(p, &ep, 0);
    if(ep == p)
      return(-1);
  }
  else  /* dp = p */
    *ival = 0;

  /* Convert fractional part, if needed */
  if(dp) {
    dp++;
    ndig = strlen(dp);

    if(ndig > 0) {
      v = strtol(dp, &ep, 0);
      if(ep == dp)
	return(-1);
      
      if(*p == '-')
	v *= -1;

      *fval = v * pow(10.0, -ndig);
    }
    else
      *fval = 0;
  }
  else
    *fval = 0;

  return(0);
}

char *sstrip (char *str) {
  char *p;

  /* Stop at the first comment character */
  p = strchr(str, '#');
  if(p)
    *p = '\0';

  /* First remove whitespace from start of string */
  while(*str != '\0' && isspace((unsigned char) *str))
    str++;

  if(*str == '\0')
    return(str);

  /* Remove whitespace from end of string */
  p = str + strlen(str) - 1;

  while(p > str && isspace((unsigned char) *p))
    p--;

  if(p == str && isspace((unsigned char) *p))
    *p = '\0';
  else
    *(p+1) = '\0';

  return(str);
}
