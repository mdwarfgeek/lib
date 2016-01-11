#include <stdlib.h>
#include <stddef.h>

#include "lfa.h"

struct filt1d_median_sort_tmp {
  float v;
  int p;
};

int filt1d_nwork (int npt, int nkern) {
  int nwork;

  nwork = npt + nkern;
  if(nwork < 3)
    nwork = 3;

  return(nwork);
}

static inline void filt1d_linear_prep (float *buf, unsigned char *mask,
                                       float *work,
                                       int npt, int nstride, int nkern) {
  float sum, xmns, xmnf;
  int i, il, ilow;

  /* Set first and last edges equal */
  il   = nkern/2;
  ilow = nkern/4;
  if(ilow < 3)
    ilow = 3;
  else
    ilow = (ilow/2)*2 + 1;

  sum = 0.0;
  for(i = 0; i < ilow; i++)
    sum += buf[i*nstride];

  xmns = sum / ilow;

  sum = 0.0;
  for(i = 0; i < ilow; i++)
    sum += buf[(npt-1-i)*nstride];

  xmnf = sum / ilow;

  /* Reflect edges before filtering */
  for(i = 0; i < il; i++) {
    work[i] = 2.0 * xmns - buf[(il+ilow-1-i)*nstride];
    work[npt+i+il] = 2.0 * xmnf - buf[(npt-i-ilow-1)*nstride];
  }

  for(i = 0; i < npt; i++)
    work[i+il] = buf[i*nstride];
}

static inline void filt1d_median_prep (float *buf, unsigned char *mask,
                                       float *work,
                                       int npt, int nstride, int nkern) {

  float xmns, xmnf;
  int i, imed, ilow;

  /* Set first and last edges equal */
  imed = nkern/2;
  ilow = nkern/4;
  if(ilow < 3)
    ilow = 3;
  else
    ilow = (ilow/2)*2 + 1;

  for(i = 0; i < ilow; i++)
    work[i] = buf[i*nstride];

  xmns = fquickselect(work, ilow/2, ilow);

  for(i = 0; i < ilow; i++)
    work[i] = buf[(npt-1-i)*nstride];

  xmnf = fquickselect(work, ilow/2, ilow);

  /* Reflect edges before filtering */
  for(i = 0; i < imed; i++) {
    work[i] = 2.0 * xmns - buf[(imed+ilow-1-i)*nstride];
    work[npt+i+imed] = 2.0 * xmnf - buf[(npt-i-ilow-1)*nstride];
  }

  for(i = 0; i < npt; i++)
    work[i+imed] = buf[i*nstride];
}

void filt1d_boxcar (float *buf, unsigned char *mask, float *work,
                    int npt, int nstride, int nkern) {
  float sum, norm;
  int i, po;

  if(npt < 4 || npt <= nkern)
    return;

  filt1d_linear_prep(buf, mask, work, npt, nstride, nkern);

  /* Boxcar filter on rest, inefficiently - can't do the obvious
     optimization of adding and subtracting in FP to avoid loss of
     precision. */
  norm = 1.0 / nkern;

  for(po = 0; po < npt; po++) {
    sum = 0;
    for(i = 0; i < nkern; i++)
      sum += work[i+po];

    buf[po*nstride] = sum * norm;
  }
}

void filt1d_hanning (float *buf, unsigned char *mask, float *work,
                     int npt, int nstride) {
  int i;

  if(npt < 4)
    return;

  filt1d_linear_prep(buf, mask, work, npt, nstride, 3);

  /* Hanning filter on rest */
  for(i = 0; i < npt; i++)
    buf[i*nstride] = 0.25*(work[i]+2*work[i+1]+work[i+2]);
}

void filt1d_median (float *buf, unsigned char *mask, float *work,
                    int npt, int nstride, int nkern) {
  struct filt1d_median_sort_tmp sortbuf[nkern], sorttmp;
  int revptr[nkern];

  float vmed, vnew, vmid;
  int i, il, ih, imed, imid, pi, po, r;

  if(npt < 4 || npt <= nkern)
    return;

  filt1d_median_prep(buf, mask, work, npt, nstride, nkern);

  /* Set first and last edges equal */
  imed = nkern/2;

  /* Load up first kernel size worth */
  for(i = 0; i < nkern; i++) {
    sortbuf[i].v = work[i];
    sortbuf[i].p = i;
  }

  /* Sort */
  fquicksort_gen(sortbuf, &sorttmp, nkern,
                 sizeof(struct filt1d_median_sort_tmp),
                 offsetof(struct filt1d_median_sort_tmp, v));

  /* Extract reverse lookup pointers */
  for(i = 0; i < nkern; i++)
    revptr[sortbuf[i].p] = i;

  /* First output pixel */
  vmed = sortbuf[imed].v;
  buf[0] = vmed;

  /* Loop over output pixels */
  for(po = 1; po < npt; po++) {
    pi = po - 1;
    r = pi % nkern;

    /* Fetch the next input pixel */
    vnew = work[nkern+pi];

    if(vnew == vmed)
      il = imed;
    else if(vnew > sortbuf[nkern-1].v)
      il = nkern;
    else {
      /* Binary chop, starting from the median we already computed,
         to locate insertion point for the new pixel value "vnew". */
      if(vnew < vmed) {
        il = 0;
        ih = imed;
      }
      else {
        il = imed+1;
        ih = nkern-1;
      }

      while(il < ih) {
        imid = (il+ih)/2;
        vmid = sortbuf[imid].v;
        if(vnew < vmid)
          ih = imid;
        else if(vnew > vmid)
          il = imid+1;
        else {
          il = imid;
          break;
        }
      }
    }

    /* Element we're removing */
    ih = revptr[r];

    if(ih < il)
      il--;

    /* il now points to insertion point, ih to element to remove */

    if(il < ih)
      /* Shift everything up to make room */
      for(i = ih; i > il; i--) {
        sortbuf[i] = sortbuf[i-1];
        revptr[sortbuf[i].p] = i;
      }
    else if(il > ih)
      /* Shift everything down to make room */
      for(i = ih; i < il; i++) {
        sortbuf[i] = sortbuf[i+1];
        revptr[sortbuf[i].p] = i;
      }
    else
      i = il;

    /* Do the insertion */
    sortbuf[i].v = vnew;
    sortbuf[i].p = r;
    revptr[r] = i;

    /* Finding the median is easy after all that work */
    vmed = sortbuf[imed].v;
    buf[po*nstride] = vmed;
  }
}
