#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <ctype.h>
#include <errno.h>

#include "util.h"

char progname[PATH_MAX+1] = { '\0' };

void setprogname (const char *name) {
  (void) strncpy(progname, name, PATH_MAX);
  progname[PATH_MAX] = '\0';
}

void fatal (int code, const char *fmt, ...) {
  va_list ap;

  if(progname[0])
    (void) fprintf(stderr, "%s: ", progname);
  va_start(ap, fmt);
  if(fmt)
    (void) vfprintf(stderr, fmt, ap);
  va_end(ap);
  (void) fprintf(stderr, "\n");

  exit(code);
}

void warning (const char *fmt, ...) {
  va_list ap;

  if(progname[0])
    (void) fprintf(stderr, "%s: ", progname);
  va_start(ap, fmt);
  if(fmt)
    (void) vfprintf(stderr, fmt, ap);
  va_end(ap);
  (void) fprintf(stderr, "\n");
}

void error (int code, const char *fmt, ...) {
  va_list ap;

  if(progname[0])
    (void) fprintf(stderr, "%s: ", progname);
  va_start(ap, fmt);
  if(fmt) {
    (void) vfprintf(stderr, fmt, ap);
    (void) fprintf(stderr, ": ");
  }
  va_end(ap);
  (void) fprintf(stderr, "%s\n", strerror(errno));

  exit(code);
}

void report_err (char *errstr, const char *fmt, ...) {
  va_list ap;

  va_start(ap, fmt);
  (void) vsnprintf(errstr, ERRSTR_LEN, fmt, ap);
  va_end(ap);
}

void report_syserr (char *errstr, const char *fmt, ...) {
  va_list ap;
  int rv;

  va_start(ap, fmt);
  rv = vsnprintf(errstr, ERRSTR_LEN, fmt, ap);
  va_end(ap);
  
  if(rv != -1 && rv < ERRSTR_LEN)
    (void) snprintf(errstr + rv, ERRSTR_LEN - rv, ": %s", strerror(errno));
}

