#include <stdlib.h>

#include "lfa.h"

/* Mike's sky level estimation method, for histogrammed 16-bit
   unsigned integer data.  Outputs returned in single precision FP.
   For efficiency, the ihmin and ihmax arguments pass in the lowest
   and highest occupied indices in the histogram, which the caller
   can keep track of when creating it (or simply give the full range
   if you don't care).  The use of "ihmin" prevents having to waste
   cycles looping past the bias / pedestal level.

   The original routine did 3 sigma upper clipping, to get this
   specify clip_low = -large number and clip_high = 3.

   This routine is optimised for sky levels close to ihmin and low
   noise levels relative to the total range.  Both are typically
   true for astronomical data. */

void skylevel (int *ihist, int ihmin, int ihmax, int mpix,
	       float clip_low, float clip_high,
	       float *skylev_r, float *sigma_r) {
  int i, s, n, iloop;
  int tclipl, tcliph, iclipl, icliph, irej, mcpix;
  float skymed, sigmed, sigma, skymedc, sigmedc, sigmac;
  float fclipl, fcliph;

  /* Find points for sigma, then median */
  n = (mpix+3)/4;

  i = ihmin;
  s = ihist[i];

  while(s <= n) {
    i++;
    s += ihist[i];
  }

  sigmed = i + (n - s) / ((float) ihist[i]) + 0.5;

  n = (mpix+1)/2;  /* guaranteed >= previous "n" */

  while(s <= n) {
    i++;
    s += ihist[i];
  }

  skymed = i + (n - s) / ((float) ihist[i]) + 0.5;

  /* the 1.48 converts MAD to rms equivalent */
  sigma = 1.48 * (skymed - sigmed);
  if(sigma < 0.5)
    sigma = 0.5;

  /* do an iterative clip to give a more robust estimator */
  iclipl = ihmin;
  icliph = ihmax;
  mcpix = mpix;
  skymedc = skymed;
  sigmac = sigma;

  for(iloop = 0; iloop < 5; iloop++) {
    /* Clip thresholds, range check and then convert to integer */
    fclipl = skymedc + clip_low * sigmac + 1;
    fcliph = skymedc + clip_high * sigmac - 1;

    if(fclipl < ihmin)
      tclipl = ihmin;
    else if(fclipl > ihmax)
      tclipl = ihmax;
    else
      tclipl = fclipl + 0.5;   /* nearest integer, positive quantity */

    if(fcliph < ihmin)
      tcliph = ihmin;
    else if(fcliph > ihmax)
      tcliph = ihmax;
    else
      tcliph = fcliph + 0.5;   /* nearest integer, positive quantity*/

    if(iloop == 0) {
      /* First time, just count between the limits.  This should
	 usually be quicker provided the clipping range is small
	 compared to [ihmin, ihmax]. */
      mcpix = 0;

      for(i = tclipl; i <= tcliph; i++)
        mcpix += ihist[i];

      iclipl = tclipl;
      icliph = tcliph;
    }
    else {
      /* How many points to reject? */
      irej = 0;

      for(i = iclipl; i < tclipl; i++)
	irej += ihist[i];

      for(i = tcliph+1; i <= icliph; i++)
	irej += ihist[i];
      
      /* Is there any work to do? */
      if(irej == 0)
	break;

      iclipl = tclipl;
      icliph = tcliph;

      /* New number of points */
      mcpix -= irej;
    }

    /* Find points for sigma, then median */
    n = (mcpix+3)/4;

    i = iclipl;
    s = ihist[i];
    
    while(s <= n) {
      i++;
      s += ihist[i];
    }
    
    sigmedc = i + (n - s) / ((float) ihist[i]) + 0.5;
    
    n = (mcpix+1)/2;  /* guaranteed >= previous "n" */
    
    while(s <= n) {
      i++;
      s += ihist[i];
    }
    
    skymedc = i + (n - s) / ((float) ihist[i]) + 0.5;

    sigmac = 1.48 * (skymedc - sigmedc);
    if(sigmac < 0.5)
      sigmac = 0.5;
  }

  *skylev_r = skymedc;
  *sigma_r = sigmac;
}

/* Simple driver routine to calculate sky level of contiguous image
   pixels, the most common case. */

void skylevel_image (float *map, unsigned char *mask, int npix,
                     float clip_low, float clip_high,
                     float *skylev, float *skynoise) {
  int hist[SKYLEVEL_DEFAULT_SIZE], nhist, hmin, hmax, v, p;
  float f;

  /* Clear histogram */
  for(v = 0; v <= SKYLEVEL_DEFAULT_ULIM; v++)
    hist[v] = 0;
  
  hmin = SKYLEVEL_DEFAULT_ULIM;
  hmax = 0;

  nhist = 0;
  
  /* Accumulate pixels */
  for(p = 0; p < npix; p++) {
    if(!mask || mask[p]) {
      f = rintf(map[p]) + SKYLEVEL_DEFAULT_OFFSET;
      
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
  
  skylevel(hist, hmin, hmax, nhist, clip_low, clip_high, skylev, skynoise);

  *skylev -= SKYLEVEL_DEFAULT_OFFSET;
}
