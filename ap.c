#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

#include "ap.h"

/* Extend arrays by this amount at a time */
#define NEXTRA_PARENTS  1024
#define NEXTRA_PIXELS   8192
#define NEXTRA_SOURCES  1024

/* Some macros to save tedious typing and long lines */
#undef ALLOC
#undef EXTEND

#define ALLOC(p, n, s)     (p) = (s *) malloc((n) * sizeof(s))
#define EXTEND(p, n, e, s) (p) = (s *) realloc((p), ((n) + (e)) * sizeof(s))

static int ap_params (struct ap *ap,
                      float *map, unsigned char *mask,
                      int minpix, double minflux,
                      float sky, float thresh,
                      int parent);

int ap_init (struct ap *ap, int nx, int ny) {
  /* Set size */
  ap->nx = nx;
  ap->ny = ny;

  ALLOC(ap->line, nx, float);

  /* Init stacks and arrays */
  ap->parents = (struct ap_parent *) NULL;
  ap->free_parents = (int *) NULL;
  ap->pixels = (struct ap_pixel *) NULL;
  ap->free_pixels = (int *) NULL;
  ap->line_parents = (int *) NULL;
  ap->sources = (struct ap_source *) NULL;

  ALLOC(ap->parents, NEXTRA_PARENTS, struct ap_parent);
  ALLOC(ap->free_parents, NEXTRA_PARENTS, int);
  if(!ap->parents || !ap->free_parents)
    goto error;

  ap->nalloc_parents = NEXTRA_PARENTS;
  ap->nfree_parents = NEXTRA_PARENTS;

  ALLOC(ap->pixels, NEXTRA_PIXELS, struct ap_pixel);
  ALLOC(ap->free_pixels, NEXTRA_PIXELS, int);
  if(!ap->pixels || !ap->free_pixels)
    goto error;

  ap->nalloc_pixels = NEXTRA_PIXELS;
  ap->nfree_pixels = NEXTRA_PIXELS;

  ALLOC(ap->line_parents, nx, int);
  if(!ap->line_parents)
    goto error;

  ALLOC(ap->sources, NEXTRA_SOURCES, struct ap_source);
  if(!ap->sources)
    goto error;

  ap->nsources = 0;
  ap->nalloc_sources = NEXTRA_SOURCES;

  ap->output = NULL;
  ap->output_user_data = NULL;

  return(0);

 error:
  ap_free(ap);

  return(-1);
}

void ap_free (struct ap *ap) {
  if(ap->line)
    free((void *) ap->line);
  if(ap->parents)
    free((void *) ap->parents);
  if(ap->free_parents)
    free((void *) ap->free_parents);
  if(ap->pixels)
    free((void *) ap->pixels);
  if(ap->free_pixels)
    free((void *) ap->free_pixels);
  if(ap->line_parents)
    free((void *) ap->line_parents);
  if(ap->sources)
    free((void *) ap->sources);
}

int ap_image (struct ap *ap,
              float *map, float *filtmap, unsigned char *mask,
              int minpix, float sky, float thresh) {
  float *prevline, *currline, *nextline, *filtline;
  unsigned char *prevlmsk, *currlmsk, *nextlmsk;
  double zfilt, zthr, minflux;
  int x, xnext, xt, y;
  int down, thisdown, left, parent, ipix, iext, ntocorr;

  int nmerge = 0;
  int npar = 0;
  int npix = 0;

  /* Put everything on the free list */
  for(parent = 0; parent < ap->nalloc_parents; parent++)
    ap->free_parents[parent] = parent;

  ap->nfree_parents = ap->nalloc_parents;

  for(ipix = 0; ipix < ap->nalloc_pixels; ipix++)
    ap->free_pixels[ipix] = ipix;

  ap->nfree_pixels = ap->nalloc_pixels;

  /* Reset last line */
  for(x = 0; x < ap->nx; x++)
    ap->line_parents[x] = -1;

  /* Clear source list */
  ap->nsources = 0;

  /* Scaled threshold */
  if(filtmap)
    zthr = sky + thresh;
  else
    zthr = 16 * (sky + thresh);

  /* Minimum isophotal flux for consideration as a source, following
     Mike's programs. */
  minflux = thresh * minpix * 1.5;

  /* Loop over image lines */
  for(y = 1; y < ap->ny-1; y++) {
    if(filtmap)
      filtline = filtmap + y * ap->nx;
    else {
      prevline = map + (y-1)*ap->nx;
      currline = map +     y*ap->nx;
      nextline = map + (y+1)*ap->nx;
      
      filtline = ap->line;
      
      if(mask) {
        prevlmsk = mask + (y-1)*ap->nx;
        currlmsk = mask +     y*ap->nx;
        nextlmsk = mask + (y+1)*ap->nx;
        
        /* Compute filtered line */
        for(x = 1; x < ap->nx-1; x++)
          filtline[x] =   prevlmsk[x-1]*prevline[x-1]
                      + 2*prevlmsk[x]*prevline[x]
                      +   prevlmsk[x+1]*prevline[x+1]
                      + 2*currlmsk[x-1]*currline[x-1]
                      + 4*currlmsk[x]*currline[x]
                      + 2*currlmsk[x+1]*currline[x+1]
                      +   nextlmsk[x-1]*nextline[x-1]
                      + 2*nextlmsk[x]*nextline[x]
                      +   nextlmsk[x+1]*nextline[x+1];
      }
      else {
        /* Compute filtered line */
        for(x = 1; x < ap->nx-1; x++)
          filtline[x] =   prevline[x-1] + 2*prevline[x] +   prevline[x+1]
                      + 2*currline[x-1] + 4*currline[x] + 2*currline[x+1]
                      +   nextline[x-1] + 2*nextline[x] +   nextline[x+1];
      }
    }

    for(x = 1; x < ap->nx-1; x = xnext) {
      xnext = x+1;

      /* Filtered value */
      zfilt = filtline[x];

      down = ap->line_parents[x];

      if(zfilt > zthr) {
        left = ap->line_parents[x-1];

        /* Collect pixels above threshold with same parents */
        for(; xnext < ap->nx-1; xnext++) {
          /* Above threshold? */
          zfilt = filtline[xnext];
          if(zfilt <= zthr)
            break;  /* nope, terminate strip */

          /* Has parent on previous row? */
          thisdown = ap->line_parents[xnext];

          if(thisdown >= 0) {
            /* Do we have one yet? */
            if(down >= 0) {
              /* Check it's the same */
              if(thisdown != down)
                break;  /* terminate strip */
            }
            else
              down = thisdown;  /* set parent for whole strip */
          }
        }

        if(down >= 0) {
          parent = down;
          
          if(left >= 0 && left != down) {
            /* Merge left into down by appending left to down */
            ap->pixels[ap->parents[down].last].next = ap->parents[left].first;
            ap->parents[down].last = ap->parents[left].last;
            
            /* Correct the line_parents array */
            ap->line_parents[x-1] = down;  /* we know this = left */
            ntocorr = ap->parents[left].refcount-1;

            /* Search for others, educated guess - work outward from
               the current position, left first. */
            for(xt = x-2; ntocorr > 0 && xt >= 1; xt--)
              if(ap->line_parents[xt] == left) {
                ap->line_parents[xt] = down;
                ntocorr--;
              }

            for(xt = x+1; ntocorr > 0 && xt < ap->nx-1; xt++)
              if(ap->line_parents[xt] == left) {
                ap->line_parents[xt] = down;
                ntocorr--;
              }

            /* Add up number of pixels and reference count */
            ap->parents[down].npixels += ap->parents[left].npixels;
            ap->parents[down].refcount += ap->parents[left].refcount;
            
            /* Put left on the free list */
            ap->free_parents[ap->nfree_parents] = left;
            ap->nfree_parents++;

            nmerge++;
          }
        }
        else {
          if(left >= 0)
            parent = left;
          else {
            /* New parent */
            if(ap->nfree_parents <= 0) {
              /* Allocate more */
              EXTEND(ap->parents,
                     ap->nalloc_parents, NEXTRA_PARENTS,
                     struct ap_parent);
              EXTEND(ap->free_parents,
                     ap->nalloc_parents, NEXTRA_PARENTS,
                     int);
              if(!ap->parents || !ap->free_parents)
                goto error;

              for(iext = 0; iext < NEXTRA_PARENTS; iext++) {
                ap->free_parents[ap->nfree_parents] = ap->nalloc_parents;
                ap->nalloc_parents++;
                ap->nfree_parents++;
              }
            }
            
            ap->nfree_parents--;
            parent = ap->free_parents[ap->nfree_parents];

            /* Init parent */
            ap->parents[parent].first = -1;
            ap->parents[parent].last = -1;
            ap->parents[parent].refcount = 0;
            ap->parents[parent].npixels = 0;

            npar++;
          }
        }
        
        /* Allocate pixel list element */
        if(ap->nfree_pixels <= 0) {
          /* Allocate more */
          EXTEND(ap->pixels,
                 ap->nalloc_pixels, NEXTRA_PIXELS,
                 struct ap_pixel);
          EXTEND(ap->free_pixels,
                 ap->nalloc_pixels, NEXTRA_PIXELS,
                 int);
          if(!ap->pixels || !ap->free_pixels)
            goto error;

          for(iext = 0; iext < NEXTRA_PIXELS; iext++) {
            ap->free_pixels[ap->nfree_pixels] = ap->nalloc_pixels;
            ap->nalloc_pixels++;
            ap->nfree_pixels++;
          }
        }
        
        ap->nfree_pixels--;
        ipix = ap->free_pixels[ap->nfree_pixels];
        
        /* Populate */
        ap->pixels[ipix].next = -1;
        ap->pixels[ipix].xa = x;
        ap->pixels[ipix].xb = xnext-1;
        ap->pixels[ipix].y = y;
        
        /* Add to tail of parent's list */
        if(ap->parents[parent].last >= 0)
          ap->pixels[ap->parents[parent].last].next = ipix;
        else  /* parent is empty, set head pointer */
          ap->parents[parent].first = ipix;

        /* Update tail pointer */
        ap->parents[parent].last = ipix;

        /* Increment counts */
        ap->parents[parent].npixels += xnext-x;
        
        for(xt = x; xt < xnext; xt++)
          if(ap->line_parents[xt] != parent) {
            /* Store reference in pointers for line */
            ap->line_parents[xt] = parent;
            ap->parents[parent].refcount++;
          }

        npix++;
      }
      else if(down >= 0) {
        /* Unref */
        ap->parents[down].refcount--;
        if(ap->parents[down].refcount <= 0) {
          /* Terminate */
          if(ap->parents[down].npixels >= minpix) {
            if(ap_params(ap, map, mask, minpix, minflux, sky, thresh, down))
              goto error;
          }
          
          /* Put all resources on free list */
          for(ipix = ap->parents[down].first;
              ipix >= 0;
              ipix = ap->pixels[ipix].next) {
            ap->free_pixels[ap->nfree_pixels] = ipix;
            ap->nfree_pixels++;
          }
          
          ap->free_parents[ap->nfree_parents] = down;
          ap->nfree_parents++;
        }
      
        ap->line_parents[x] = -1;
      }
    }
  }

  /* Terminate all remaining parents.  No need to free anything. */
  for(x = 1; x < ap->nx-1; x++) {
    down = ap->line_parents[x];

    if(down >= 0) {
      if(ap->parents[down].npixels >= minpix) {
        if(ap_params(ap, map, mask, minpix, minflux, sky, thresh, down))
          goto error;

        ap->parents[down].npixels = 0;  /* prevents repeat */
      }
    }
  }

  return(0);

 error:

  return(-1);
}

static int ap_params (struct ap *ap,
                      float *map, unsigned char *mask,
                      int minpix, double minflux,
                      float sky, float thresh,
                      int parent) {
  int ipix, ix, p;
  int xmin, xmax, ymin, ymax;
  double x, y, z;
  double sx, sy, sxx, syy, sxy, sz, zmax;
  int areal[AP_NAREAL], b, as;
  unsigned char flags;

  struct ap_source *obj, objbuf;

  /* Walk list of pixels, accumulating sums */
  sx = 0;
  sy = 0;
  sxx = 0;
  syy = 0;
  sxy = 0;
  sz = 0;

  zmax = 0;

  for(b = 0; b < AP_NAREAL; b++)
    areal[b] = 0;

  flags = 0;

  /* Init bounding box calculation */
  ipix = ap->parents[parent].first;

  xmin = ap->nx-1;
  xmax = 0;
  ymin = ap->ny-1;
  ymax = 0;

  for(ipix = ap->parents[parent].first;
      ipix >= 0;
      ipix = ap->pixels[ipix].next) {
    y = ap->pixels[ipix].y + 1;

    /* Loop through pixels */
    p = ap->pixels[ipix].y * ap->nx + ap->pixels[ipix].xa;
    for(ix = ap->pixels[ipix].xa; ix <= ap->pixels[ipix].xb; ix++, p++) {
      x = ix + 1;

      /* Check for bad pixels */
      if(mask && !mask[p]) {
        flags |= AP_SOURCE_FLAG_BADPIX;
        continue;  /* skip */
      }

      /* Uses original image pixels (without convolution) for computing
         source parameters. */
      z = map[p] - sky;

      if(z > 0) {
        sx  += x*z;
        sy  += y*z;
        sxx += (x*x + 1.0/12.0)*z;  /* 1/12 finite pixel size correction */
        syy += (y*y + 1.0/12.0)*z;  /* See Irwin 1985 */
        sxy += x*y*z;
        sz  += z;
        
        /* Accumulate peak height */
        if(z > zmax)
          zmax = z;
        
        /* Accumulate areal profiles.  We need only fill in the highest
           this object contributes to, and then sum down later. */
#if FLT_RADIX == 2
        b = ilogb(z / thresh);
#else
        b = floor(log2(z / thresh));
#endif

        if(b >= 0) {
          if(b >= AP_NAREAL)
            b = AP_NAREAL-1;
          
          areal[b]++;

          /* Update bounding box */
          if(ix < xmin)
            xmin = ix;
          
          if(ix > xmax)
            xmax = ix;

          if(ap->pixels[ipix].y < ymin)
            ymin = ap->pixels[ipix].y;
          
          if(ap->pixels[ipix].y > ymax)
            ymax = ap->pixels[ipix].y;
        }
      }
    }
  }

  /* Areal profiles, sum down from the top. */
  as = 0;
  
  for(b = AP_NAREAL-1; b >= 0; b--) {
    as += areal[b];
    areal[b] = as;
  }

  if(as >= minpix && sz >= minflux) {
    if(ap->output)
      obj = &objbuf;
    else {
      /* Allocate source */
      if(ap->nalloc_sources <= ap->nsources) {
        EXTEND(ap->sources,
               ap->nalloc_sources, NEXTRA_SOURCES,
               struct ap_source);
        if(!ap->sources)
          goto error;
        
        ap->nalloc_sources += NEXTRA_SOURCES;
      }
      
      obj = ap->sources + ap->nsources;
    }

    /* 1st moments */
    obj->x = sx / sz;
    obj->y = sy / sz;

    /* 2nd moments */
    obj->sxx = sxx / sz - obj->x*obj->x;
    obj->syy = syy / sz - obj->y*obj->y;
    obj->sxy = sxy / sz - obj->x*obj->y;

    /* Isophotal flux, peak height */
    obj->ziso = sz;
    obj->zpk = zmax;

    /* Areal profiles */
    for(b = 0; b < AP_NAREAL; b++)
      obj->areal[b] = areal[b];

    /* Bounding box */
    obj->bb[0] = xmin+1;
    obj->bb[1] = xmax+1;
    obj->bb[2] = ymin+1;
    obj->bb[3] = ymax+1;

    /* Flags */
    obj->flags = flags;

    ap->nsources++;

    if(ap->output) {
      if(ap->output(ap, obj, ap->output_user_data))
        goto error;
    }
  }

  return(0);

 error:

  return(-1);
}

void ap_ellipse (struct ap_source *obj,
                 double *gfwhm, double *ell, double *spa, double *cpa) {
  double t, d, tmp, tspa, tcpa, norm;

  /* Trace and determinant of covariance matrix */
  t = obj->sxx + obj->syy;
  d = obj->sxx*obj->syy - obj->sxy*obj->sxy;
  
  /* Gaussian FWHM */
  *gfwhm = sqrt(8*M_LN2*t);
  
  /* Ellipticity is 1 - sqrt(ratio of eigenvalues)
     Eigenvalues are 0.5*(t + tmp) and 0.5*(t - tmp) */
  tmp = sqrt(t*t - 4*d);
  *ell = 1.0 - sqrt((t - tmp) / (t + tmp));
  
  /* sin and cos of PA */
  tspa = obj->sxy + obj->sxy;
  tcpa = tmp + obj->sxx - obj->syy;
  
  norm = hypot(tspa, tcpa);
  
  if(norm > 0) {
    *spa = tspa / norm;
    *cpa = tcpa / norm;
  }
  else {
    *spa = 0;
    *cpa = 0;
  }
}

