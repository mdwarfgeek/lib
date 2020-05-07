#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "lfa.h"
#include "util.h"

int backremove (float *mapin, unsigned char *mask, float *mapout,
                int nx, int ny, int nbsize) {
  int nbx, nby, nbin, nbsizex, nbsizey, hbsizex, hbsizey, nxbord, nybord;
  int minpix, nwork;

  int medwin;

  float *bgmap = (float *) NULL;
  unsigned char *bgmask = (unsigned char *) NULL;
  float *tmpbuf = (float *) NULL;
  float *work = (float *) NULL;

  int hist[SKYLEVEL_DEFAULT_SIZE], nhist, hmin, hmax;
  int nmed;

  int bx, by, ox, oy, x, y, p, pb, v;
  float f, skylev, skynoise;

  float rsizex, rsizey, dx, dy;
  int bxp, byp;

  float bval;

  /* Decide number of bins */
  nbx = nx / nbsize;
  if(nbx < 1)
    nbx = 1;

  nby = ny / nbsize;
  if(nby < 1)
    nby = 1;

  /* Actual bin size */
  nbsizex = nx / nbx;
  nbsizey = ny / nby;

  hbsizex = nbsizex / 2;
  hbsizey = nbsizey / 2;

  /* Crop border size to ensure exact number of pixels */
  nxbord = (nx - nbx * nbsizex) / 2;
  nybord = (ny - nby * nbsizey) / 2;

  /* Minimum pixels need to be populated in each bin */
  minpix = nbsizex*nbsizey / 4;
  if(minpix < 1)
    minpix = 1;

  /* Window size for median filter */
  medwin = MIN(nbx, nby) / 2;
  if(medwin < 1)
    medwin = 1;
  else if(medwin > 5)
    medwin = 5;
  else if(medwin % 2 == 0)
    medwin++;

  /* Allocate workspace */
  nbin = nbx * nby;
  nwork = filt1d_nwork(MAX(nbx, nby), MAX(medwin, 3));

  bgmap = (float *) malloc(nbin * sizeof(float));
  bgmask = (unsigned char *) malloc(nbin * sizeof(unsigned char));
  tmpbuf = (float *) malloc(nbin * sizeof(float));
  work = (float *) malloc(nwork * sizeof(float));
  if(!bgmap || !bgmask || !tmpbuf || !nwork)
    goto error;

  hmin = 0;
  hmax = SKYLEVEL_DEFAULT_ULIM;

  /* Compute background in each bin */
  nmed = 0;

  for(by = 0; by < nby; by++) {
    oy = nybord + by * nbsizey;

    for(bx = 0; bx < nbx; bx++) {
      ox = nxbord + bx * nbsizex;

      pb = by*nbx + bx;

      /* Clear histogram */
      for(v = hmin; v <= hmax; v++)
        hist[v] = 0;

      hmin = SKYLEVEL_DEFAULT_ULIM;
      hmax = 0;

      nhist = 0;

      /* Accumulate pixels */
      for(y = 0; y < nbsizey; y++)
        for(x = 0; x < nbsizex; x++) {
          p = (oy+y)*nx + (ox+x);

          if(!mask || mask[p]) {
            f = rintf(mapin[p]) + SKYLEVEL_DEFAULT_OFFSET;

            if(f < 0)
              v = 0;
            else if(f > SKYLEVEL_DEFAULT_ULIM)
              v = SKYLEVEL_DEFAULT_ULIM;
            else
              v = f;

            hist[v]++;
            nhist++;

            if(v < hmin)
              hmin = v;
            if(v > hmax)
              hmax = v;
          }
        }

      if(nhist >= minpix) {
        skylevel(hist, hmin, hmax, nhist, -FLT_MAX, 3.0, &skylev, &skynoise);

        bgmap[pb] = skylev - SKYLEVEL_DEFAULT_OFFSET;
        bgmask[pb] = 1;

        tmpbuf[nmed] = skylev - SKYLEVEL_DEFAULT_OFFSET;
        nmed++;
      }
      else {
        bgmap[pb] = 0;
        bgmask[pb] = 0;
      }
    }
  }

  /* Compute median sky level */
  fmedsig(tmpbuf, nmed, &skylev, (float *) NULL);

  if(medwin > 1) {
    /* Median filter rows */
    for(by = 0; by < nby; by++)
      filt1d_median(bgmap+by*nbx, bgmask+by*nbx, work, nbx, 1, medwin);
    
    /* Median filter columns */
    for(bx = 0; bx < nbx; bx++)
      filt1d_median(bgmap+bx, bgmask+bx, work, nby, nbx, medwin);

    /* Hanning filter rows */
    for(by = 0; by < nby; by++)
      filt1d_hanning(bgmap+by*nbx, bgmask+by*nbx, work, nbx, 1);
    
    /* Hanning filter columns */
    for(bx = 0; bx < nbx; bx++)
      filt1d_hanning(bgmap+bx, bgmask+bx, work, nby, nbx);
  }

  /* Make map of sky background variations using bilinear interpolation.
     Edges are held constant, rather than extrapolating. */
  rsizex = 1.0 / nbsizex;
  rsizey = 1.0 / nbsizey;

  for(y = 0; y < ny; y++) {
    by = (y - hbsizey - nybord) / nbsizey;
    byp = (y + hbsizey - nybord) / nbsizey;

    if(by < 0)
      by = 0;

    if(byp > nby-1)
      byp = nby-1;

    dy = (y - nybord - nbsizey * by - hbsizey) * rsizey;

    for(x = 0; x < nx; x++) {
      bx = (x - hbsizex - nxbord) / nbsizex;
      bxp = (x + hbsizex - nxbord) / nbsizex;

      if(bx < 0)
        bx = 0;

      if(bxp > nbx-1)
        bxp = nbx-1;

      dx = (x - nxbord - nbsizex * bx - hbsizex) * rsizex;

      bval = (1.0 - dx) * ((1.0 - dy) * bgmap[by*nbx+bx] +
                                   dy * bgmap[byp*nbx+bx])
           +         dx * ((1.0 - dy) * bgmap[by*nbx+bxp] +
                                   dy * bgmap[byp*nbx+bxp])
           - skylev;

      mapout[y*nx+x] = mapin[y*nx+x] - bval;
    }
  }

  free((void *) bgmap);
  free((void *) bgmask);
  free((void *) tmpbuf);
  free((void *) work);

  return(0);

 error:
  if(bgmap)
    free((void *) bgmap);
  if(bgmask)
    free((void *) bgmask);
  if(tmpbuf)
    free((void *) tmpbuf);
  if(work)
    free((void *) work);

  return(-1);
}

