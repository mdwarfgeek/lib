#ifndef FFTUTIL_H
#define FFTUTIL_H

#include "lfa.h"

int fftwf_convolve_2d (float *kern, int nkx, int nky,
                       float **mapsin, float **mapsout,
                       int nxmap, int nymap, int nmaps);

#endif  /* FFTUTIL_H */
