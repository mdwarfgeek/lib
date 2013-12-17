#include <sys/types.h>
#include <unistd.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <pwd.h>

#include "cnf.h"
#include "lfa.h"

#include "cvtunit.h"
#include "util.h"

int cnf_read (char *filename,
	      struct cnf_field *fields, int nfields,
	      char *errstr) {
  FILE *fp;
  char buf[16384], *p, *ep, *nep, *sp, *dp;
  int iline;

  char section[1024] = { '\0' };

  int f;

  int *iptr;
  float *fptr;
  double *dptr;
  char *sptr;

  int rv;
  struct passwd *pw;

  /* Open file */
  fp = fopen(filename, "r");
  if(!fp) {
    report_syserr(errstr, "open: %s", filename);
    goto error;
  }

  /* Read */
  iline = 0;
  while(fgets(buf, sizeof(buf), fp)) {
    /* Remove comments */
    p = strchr(buf, '#');
    if(p)
      *p = '\0';

    /* Remove leading and trailing whitespace */
    p = sstrip(buf);

    /* ini file format: either a section in [] */
    if(*p == '[') {
      p++;

      ep = strchr(p, ']');
      if(!ep) {
	report_err(errstr, "parse error on %s line %d: [ without ] in %s",
		   filename, iline+1, buf);
	goto error;
      }

      *ep = '\0';

      /* Interpret as section name (note length limit) */
      strncpy(section, p, sizeof(section)-1);
      section[sizeof(section)-1] = '\0';
    }
    else if(*p) {  /* or a key = value pair */
      ep = strchr(p, '=');
      if(!ep) {
	report_err(errstr, "parse error on %s line %d: no = in %s",
		   filename, iline+1, p);
	goto error;
      }

      *ep = '\0';
      ep++;

      p = sstrip(p);  /* keyword */
      ep = sstrip(ep);  /* value */

      /* Table lookup */
      for(f = 0; f < nfields; f++) {
	if(!strcasecmp(section, fields[f].section) &&
	   !strcasecmp(p, fields[f].key)) {
	  /* Match; repeats consume one item at a time for numeric fields */
	  switch(fields[f].type) {
	  case CNF_FIELD_INT:
	    iptr = (int *) fields[f].ptr;
	    *iptr = strtol(ep, &nep, 0);
	    if(ep == nep) {
	      report_err(errstr,
			 "parse error on %s line %d: "
			 "invalid numeric value %s for %s",
			 filename, iline+1, ep, p);
	      goto error;
	    }

	    ep = nep;  /* advance pointer */

	    break;
	  case CNF_FIELD_FLOAT:
	    fptr = (float *) fields[f].ptr;
	    *fptr = strtod(ep, &nep);
	    if(ep == nep) {
	      report_err(errstr,
			 "parse error on %s line %d: "
			 "invalid numeric value %s for %s",
			 filename, iline+1, ep, p);
	      goto error;
	    }

	    ep = nep;  /* advance pointer */

	    break;
	  case CNF_FIELD_DOUBLE:
	    dptr = (double *) fields[f].ptr;
	    *dptr = strtod(ep, &nep);
	    if(ep == nep) {
	      report_err(errstr,
			 "parse error on %s line %d: "
			 "invalid numeric value %s for %s",
			 filename, iline+1, ep, p);
	      goto error;
	    }

	    ep = nep;  /* advance pointer */

	    break;
	  case CNF_FIELD_DBL60:
	    dptr = (double *) fields[f].ptr;
	    rv = base60_to_10(ep, &nep, ":", UNIT_DEG, dptr, UNIT_DEG);
	    if(rv) {
	      rv = base60_to_10(ep, &nep, " ", UNIT_DEG, dptr, UNIT_DEG);
	      if(rv) {
		*dptr = strtod(ep, &nep);
		if(ep == nep) {
		  report_err(errstr,
			     "parse error on %s line %d: "
			     "invalid angle %s for %s",
			     filename, iline+1, ep, p);
		  goto error;
		}
	      }
	    }

	    ep = nep;  /* advance pointer */

	    break;
	  case CNF_FIELD_STRING:
	    sptr = (char *) fields[f].ptr;
	    strncpy(sptr, ep, fields[f].len-1);
	    sptr[fields[f].len-1] = '\0';
	    
	    break;
	  case CNF_FIELD_PATH:
	    sptr = (char *) fields[f].ptr;

	    /* Do we need to do tilde expansion? */
	    if(*ep == '~') {
	      ep++;

	      /* Find the next /, if there is one */
	      sp = strchr(ep, '/');
	      if(sp)
		*sp = '\0';

	      if(*ep == '\0') {
		/* Current user, first try $HOME */
		dp = getenv("HOME");
		if(!dp) {
		  /* Nope, try password file */
		  pw = getpwuid(getuid());
		  if(!pw)
		    goto error;
		  else
		    dp = pw->pw_dir;
		}
	      }
	      else {
		pw = getpwnam(ep);
		if(!pw)
		  goto error;
		else
		  dp = pw->pw_dir;
	      }

	      if(sp) {
		/* There was some more, concatenate (the lazy way) */
		sp++;
		snprintf(sptr, fields[f].len, "%s/%s", dp, sp);
	      }
	      else {
		strncpy(sptr, dp, fields[f].len-1);
		sptr[fields[f].len-1] = '\0';
	      }
	    }
	    else {
	      strncpy(sptr, ep, fields[f].len-1);
	      sptr[fields[f].len-1] = '\0';
	    }	    

	    break;
	  }

	  fields[f].found++;
	}
      }

      /* ignore the ones we don't know */
    }
    /* or nothing, in which case we do nothing */

    iline++;  /* increment running line count */
  }

  if(ferror(fp)) {
    report_syserr(errstr, "read");
    goto error;
  }

  fclose(fp);

  return(0);

 error:
  return(-1);
}

