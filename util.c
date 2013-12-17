#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>

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

int daemonify (void) {
  int rv, pid, fd;
  struct sigaction sa, osa;

  /* Block SIGHUP temporarily so we don't blow up the child when
     the parent exits. */
  memset(&sa, 0, sizeof(sa));
  sa.sa_handler = SIG_IGN;
  sigemptyset(&(sa.sa_mask));
  sa.sa_flags = 0;

  rv = sigaction(SIGHUP, &sa, &osa);
  if(rv < 0)
    return(rv);

  /* Now fork */
  pid = fork();
  if(pid < 0)
    return(pid);
  else if(pid)
    /* Parent exits */
    _exit(0);
  else {
    /* In the child: become the process group leader */
    rv = setsid();
    if(rv < 0)
      return(rv);

    /* Change to root directory in case somebody unmounts */
    rv = chdir("/");
    if(rv < 0)
      return(rv);

    /* Open /dev/null */
    fd = open("/dev/null", O_RDWR);
    if(fd < 0)
      return(fd);

    /* Replace stdin, stdout and stderr with /dev/null to detach
       from controlling tty. */
    dup2(fd, 0);
    dup2(fd, 1);
    dup2(fd, 2);

    close(fd);

    /* Restore SIGHUP */
    rv = sigaction(SIGHUP, &osa, (struct sigaction *) NULL);
    if(rv < 0)
      return(rv);
  }

  return(0);
}
