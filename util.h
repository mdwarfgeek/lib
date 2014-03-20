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

#undef SWAP
#define SWAP(x, y, t) (t) = (x); (x) = (y); (y) = (t)

#undef BSWAP16
#define BSWAP16(p, t) \
  SWAP(*(((uint8_t *) (p)) + 0), *(((uint8_t *) (p)) + 1), t)

#undef BSWAP32
#define BSWAP32(p, t) \
  SWAP(*(((uint8_t *) (p)) + 0), *(((uint8_t *) (p)) + 3), t); \
  SWAP(*(((uint8_t *) (p)) + 1), *(((uint8_t *) (p)) + 2), t);

#undef BSWAP64
#define BSWAP64(p, t) \
  SWAP(*(((uint8_t *) (p)) + 0), *(((uint8_t *) (p)) + 7), t); \
  SWAP(*(((uint8_t *) (p)) + 1), *(((uint8_t *) (p)) + 6), t); \
  SWAP(*(((uint8_t *) (p)) + 2), *(((uint8_t *) (p)) + 5), t); \
  SWAP(*(((uint8_t *) (p)) + 3), *(((uint8_t *) (p)) + 4), t);

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

int daemonify (void);

/* In libc, but not the headers */

char *basename (const char *);
char *dirname (const char *);

#endif /* __UTIL_H__ */
