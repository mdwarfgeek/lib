#ifndef __UTIL_H__
#define __UTIL_H__

/* --Constants-- */

/* Length of error message buffers */
#define ERRSTR_LEN 1024

/* --Macros-- */

#undef MIN
#define MIN(a, b) ((a) > (b) ? (b) : (a))

#undef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#undef ABS
#define ABS(x) ((x) >= 0 ? (x) : -(x))

#undef SWAP
#define SWAP(x, y, t) (t) = (x); (x) = (y); (y) = (t)

#define ARRAYLEN(a) (sizeof(a)/sizeof((a)[0]))

#define NINT(a) ((int) ((a)+((a) < 0 ? -0.5 : 0.5)))
#define NLNG(a) ((long) ((a)+((a) < 0 ? -0.5 : 0.5)))

/* --Variables-- */

extern char progname[];

/* --Functions-- */

/* Utility functions */

void setprogname (const char *name);

void fatal (int code, const char *fmt, ...);
void warning (const char *fmt, ...);
void error (int code, const char *fmt, ...);

void report_err (char *errstr, const char *fmt, ...);
void report_syserr (char *errstr, const char *fmt, ...);

/* In libc, but not the headers */

char *basename (const char *);
char *dirname (const char *);

#endif /* __UTIL_H__ */
