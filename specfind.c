#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "specfind.h"
#include "lfa.h"

/* Extend arrays by this amount at a time */
#define NEXTRA_TERM    16
#define NEXTRA_SOURCES 16

/* Some macros to save tedious typing and long lines */
#undef ALLOC
#undef EXTEND

#define ALLOC(p, n, s)     (p) = (s *) malloc((n) * sizeof(s))
#define EXTEND(p, n, e, s) (p) = (s *) realloc((p), ((n) + (e)) * sizeof(s))

/* Internal structures */
struct term {
  int xl;
  int xh;
};

struct termlist {
  float *line;
  unsigned char *mask;
  int nx;
  int minpix;
  float sky;
  float thresh;

  struct term *list;
  int nterm;   /* stack pointer */
  int nalloc;  /* stack length */

  struct specfind_line *sources;
  int nsources;
  int nalloc_sources;
};

static int do_overlp (float *filtline, int xl, int xh,
                      int minpix, float zthr,
                      struct termlist *termlist,
                      int depth);
static int do_params (float *filtline, int xl, int xh,
                      int minpix, float zthr,
                      struct termlist *termlist);

int specfind (float *line, unsigned char *mask, int nx,
              int minpix, float sky, float thresh,
              struct specfind_line **list, int *nlist) {
  float *filtline = (float *) NULL;
  double zfilt, zthr;
  int x, xnext, npix;

  struct termlist termlist;

  termlist.line = line;
  termlist.mask = mask;
  termlist.nx = nx;
  termlist.minpix = minpix;
  termlist.sky = sky;
  termlist.thresh = thresh;

  termlist.list = (struct term *) NULL;
  termlist.nterm = 0;
  termlist.nalloc = NEXTRA_TERM;

  termlist.sources = (struct specfind_line *) NULL;
  termlist.nsources = 0;
  termlist.nalloc_sources = NEXTRA_SOURCES;

  /* Line buffers */
  ALLOC(filtline, nx, float);
  if(!filtline)
    goto error;

  /* Termination list */
  ALLOC(termlist.list, termlist.nalloc, struct term);
  if(!termlist.list)
    goto error;

  /* Source list */
  ALLOC(termlist.sources, termlist.nalloc_sources, struct specfind_line);
  if(!termlist.sources)
    goto error;

  /* Scaled threshold */
  zthr = 4 * (sky + thresh);

  if(mask) {
    /* Compute filtered line */
    for(x = 1; x < nx-1; x++)
      filtline[x] = mask[x-1]*line[x-1]
                  + 2*mask[x]*line[x]
                  + mask[x+1]*line[x+1];
  }
  else {
    /* Compute filtered line */
    for(x = 1; x < nx-1; x++)
      filtline[x] = line[x-1] + 2*line[x] + line[x+1];
  }

  for(x = 1; x < nx-1; x = xnext) {
    xnext = x+1;

    /* Filtered value */
    zfilt = filtline[x];

    if(zfilt > zthr) {
      for(; xnext < nx-1; xnext++) {
        /* Above threshold? */
        zfilt = filtline[xnext];
        if(zfilt <= zthr)
          break;  /* nope, terminate strip */
      }

      /* Termination */
      npix = xnext - x;

      if(npix >= minpix) {
        do_overlp(filtline, x, xnext-1,
                  minpix, zthr,
                  &termlist,
                  0);
      }
    }
  }

  *list = termlist.sources;
  *nlist = termlist.nsources;

  free((void *) filtline);
  filtline = (float *) NULL;
  free((void *) termlist.list);
  termlist.list = (struct term *) NULL;

  return(0);

 error:
  if(filtline)
    free((void *) filtline);
  if(termlist.list)
    free((void *) termlist.list);
  if(termlist.sources)
    free((void *) termlist.sources);

  return(-1);
}

#define THRSTEP  M_SQRT2

static int do_overlp (float *filtline, int xl, int xh,
                      int minpix, float zthr,
                      struct termlist *termlist,
                      int depth) {
  float zfilt;
  int x, xnext, npix, nsrc, rv;

  int iterm, bterm;

  nsrc = 0;

  bterm = termlist->nterm;  /* base pointer in termination stack */

  for(x = xl; x <= xh; x = xnext) {
    xnext = x+1;
    
    /* Filtered value */
    zfilt = filtline[x];
    
    if(zfilt > zthr) {
      for(; xnext <= xh; xnext++) {
        /* Above threshold? */
        zfilt = filtline[xnext];
        if(zfilt <= zthr)
          break;  /* nope, terminate strip */
      }
      
      npix = xnext - x;

      if(npix >= minpix) {
        rv = do_overlp(filtline, x, xnext-1,
                       minpix, zthr * THRSTEP,
                       termlist,
                       depth+1);
        if(rv == 0)
          rv = 1;

        if(rv == 1) {
          /* Add to termination stack */
          if(termlist->nterm+1 >= termlist->nalloc) {
            EXTEND(termlist->list, termlist->nalloc, NEXTRA_TERM, struct term);
            if(!termlist->list)
              return(-1);

            termlist->nalloc += NEXTRA_TERM;
          }

          termlist->list[termlist->nterm].xl = x;
          termlist->list[termlist->nterm].xh = xnext-1;
          termlist->nterm++;
        }

        nsrc += rv;
      }
    }
  }

  if(depth == 0 || nsrc > 1) {
    /* Process any pending terminations */
    for(iterm = bterm; iterm < termlist->nterm; iterm++) {
      do_params(filtline, termlist->list[iterm].xl, termlist->list[iterm].xh,
                minpix, zthr,
                termlist);
    }
  }

  termlist->nterm = bterm;

  return(nsrc);
}

static int do_params (float *filtline, int xl, int xh,
                      int minpix, float zthr,
                      struct termlist *termlist) {
  int ix;
  double x, z;
  double sx, sz, zmax;
  unsigned char flags;

  struct specfind_line *obj;

  /* Loop through pixels */
  sx = 0;
  sz = 0;

  zmax = 0;

  flags = 0;

  for(ix = xl; ix <= xh; ix++) {
    x = ix + 1;

    /* Check for bad pixels */
    if(termlist->mask && !termlist->mask[ix]) {
      flags = 1;
      continue;  /* skip */
    }

    /* Uses original image pixels (without convolution) for computing
       source parameters. */
    z = termlist->line[ix] - termlist->sky;

    if(z > 0) {
      sx += x*z;
      sz += z;
      
      /* Accumulate peak height */
      if(z > zmax)
        zmax = z;
    }
  }

  /* Allocate source */
  if(termlist->nalloc_sources <= termlist->nsources) {
    EXTEND(termlist->sources,
           termlist->nalloc_sources, NEXTRA_SOURCES,
           struct specfind_line);
    if(!termlist->sources)
      goto error;
    
    termlist->nalloc_sources += NEXTRA_SOURCES;
  }
  
  obj = termlist->sources + termlist->nsources;
  
  /* 1st moments */
  obj->x = sx / sz;
  
  /* Isophotal flux, peak height */
  obj->ziso = sz;
  obj->zpk = zmax;
  
  /* Bounding box */
  obj->xl = xl + 1;
  obj->xh = xh + 1;

  /* Flags */
  obj->flags = flags;
  
  termlist->nsources++;

  return(0);

 error:

  return(-1);
}

