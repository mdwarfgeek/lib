/* fftutil.c:  Various FFT utility routines that depend on FFTW. */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <fftw3.h>

int fftwf_convolve_2d (float *kern, int nkx, int nky,
                       float **mapsin, float **mapsout,
                       int nxmap, int nymap, int nmaps) {
  int ntx, nty, ntxr, ntxs, nallocc, nallocf;
  float norm;

  fftwf_complex *abuf = (fftwf_complex *) NULL;
  fftwf_complex *bbuf = (fftwf_complex *) NULL;
  fftwf_complex *cbuf = (fftwf_complex *) NULL;
  float *rabuf, *rbbuf, *rcbuf;
  fftwf_plan pa, pb, pc;

  int hkx, hky;

  float *mapin, *mapout;
  int imap;

  int x, y, xo, yo, p;

  /* Size of convolution */
  ntx = nxmap + nkx - 1;
  nty = nymap + nky - 1;
  
  /* Make sure even */
  if(ntx % 2)
    ntx++;

  if(nty % 2)
    nty++;

  ntxr = (ntx/2) + 1;
  ntxs = 2*ntxr;

  nallocc = nty*ntxr;
  nallocf = 2*nallocc;

  /* Normalisation */
  norm = 1.0 / (ntx*nty);

  /* Allocate workspace */
  abuf = (fftwf_complex *) fftwf_malloc(nallocc * sizeof(fftwf_complex));
  bbuf = (fftwf_complex *) fftwf_malloc(nallocc * sizeof(fftwf_complex));
  cbuf = (fftwf_complex *) fftwf_malloc(nallocc * sizeof(fftwf_complex));
  if(!abuf || !bbuf || !cbuf)
    goto error;

  rabuf = (float *) abuf;
  rbbuf = (float *) bbuf;
  rcbuf = (float *) cbuf;

  /* Make plans */
  pa = fftwf_plan_dft_r2c_2d(nty, ntx, rabuf, abuf, FFTW_ESTIMATE);
  pb = fftwf_plan_dft_r2c_2d(nty, ntx, rbbuf, bbuf, FFTW_ESTIMATE);
  pc = fftwf_plan_dft_c2r_2d(nty, ntx, cbuf, rcbuf, FFTW_ESTIMATE);

  if(!pa || !pb || !pc)
    goto error;

  /* Kernel, zero-padded to size */
  memset(bbuf, 0, nallocc*sizeof(fftwf_complex));

  hkx = (nkx-1) / 2;
  hky = (nky-1) / 2;

  for(y = 0; y < nky; y++) {
    yo = (y + nty - hky) % nty;
    for(x = 0; x < nkx; x++) {
      xo = (x + ntx - hkx) % ntx;
      rbbuf[yo*ntxs+xo] = kern[y*nkx+x];
    }
  }

  /* Forward transform */
  fftwf_execute(pb);

  for(imap = 0; imap < nmaps; imap++) {
    mapin = mapsin[imap];
    mapout = mapsout[imap];

    /* Image, zero-padded */
    for(y = 0; y < nymap; y++) {
      for(x = 0; x < nxmap; x++)
        rabuf[y*ntxs+x] = mapin[y*nxmap+x];
      
      for(; x < ntxs; x++)
        rabuf[y*ntxs+x] = 0;
    }
    for(; y < nty; y++)
      for(x = 0; x < ntxs; x++)
        rabuf[y*ntxs+x] = 0;
    
    /* Forward transform */
    fftwf_execute(pa);
    
    /* Convolution = A B */
    for(p = 0; p < nallocc; p++) {
      cbuf[p][0] = abuf[p][0] * bbuf[p][0] - abuf[p][1] * bbuf[p][1];
      cbuf[p][1] = abuf[p][1] * bbuf[p][0] + abuf[p][0] * bbuf[p][1];
    }
    
    /* Inverse transform */
    fftwf_execute(pc);
    
    /* Copy out result */
    for(y = 0; y < nymap; y++)
      for(x = 0; x < nxmap; x++)
        mapout[y*nxmap+x] = rcbuf[y*ntxs+x] * norm;
  }

  /* Free */
  fftwf_destroy_plan(pa);
  fftwf_destroy_plan(pb);
  fftwf_destroy_plan(pc);
  
  fftwf_free(abuf);
  abuf = (fftwf_complex *) NULL;
  
  fftwf_free(bbuf);
  bbuf = (fftwf_complex *) NULL;
  
  fftwf_free(cbuf);
  cbuf = (fftwf_complex *) NULL;

  return(0);

 error:
  if(abuf)
    fftwf_free((void *) abuf);
  if(bbuf)
    fftwf_free((void *) bbuf);
  if(cbuf)
    fftwf_free((void *) cbuf);

  return(-1);
}

