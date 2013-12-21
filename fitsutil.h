#ifndef FITSUTIL_H
#define FITSUTIL_H

#include <fitsio.h>

void fitsio_err (char *errstr, int status, const char *fmt, ...);
int read_wcs (fitsfile *fits, struct wcs_info *wcs, int verbose, char *errstr);

#endif  /* FITSUTIL_H */
