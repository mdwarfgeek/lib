#ifndef UTIL_H
#define UTIL_H

#include <stdint.h>

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
#define BSWAP16(p) {                                            \
  uint8_t t;                                                    \
  SWAP(*(((uint8_t *) (p)) + 0), *(((uint8_t *) (p)) + 1), t);  \
}

#undef BSWAP32

#if defined(__GNUC__) && (defined(__i386) || defined(__amd64))
#define BSWAP32(p) \
  asm("bswapl %0" : "=r" (*((uint32_t *) (p))) : "0" (*((uint32_t *) (p))))
#else
#define BSWAP32(p) {                                            \
  uint8_t t;                                                    \
  SWAP(*(((uint8_t *) (p)) + 0), *(((uint8_t *) (p)) + 3), t);  \
  SWAP(*(((uint8_t *) (p)) + 1), *(((uint8_t *) (p)) + 2), t);  \
}
#endif

#undef BSWAP64

#if defined(__GNUC__) && defined(__amd64)
#define BSWAP64(p) \
  asm("bswapq %0" : "=r" (*((uint64_t *) (p))) : "0" (*((uint64_t *) (p))))
#else
#define BSWAP64(p) {                                            \
  uint8_t t;                                                    \
  SWAP(*(((uint8_t *) (p)) + 0), *(((uint8_t *) (p)) + 7), t);  \
  SWAP(*(((uint8_t *) (p)) + 1), *(((uint8_t *) (p)) + 6), t);  \
  SWAP(*(((uint8_t *) (p)) + 2), *(((uint8_t *) (p)) + 5), t);  \
  SWAP(*(((uint8_t *) (p)) + 3), *(((uint8_t *) (p)) + 4), t);  \
}
#endif

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

#endif /* UTIL_H */
