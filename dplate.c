#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lfa.h"

/* Solve for standard 4 or 6-coefficient linear transformation between
   lists of 2-D x,y coordinates.  The forward transformation is from
   "ref" to "com".  In matrices, the transformation is:

   (x_com) = (tr[0] tr[1]) (x_ref) + (tr[2])
   (y_com) = (tr[4] tr[3]) (y_ref)   (tr[5])

   Please note the "reversed" ordering of indices 3 and 4.  This was
   kept for compatibility with the original (Fortran) cplate and has
   become somewhat ingrained.  If ncoeff=6 the full 6-coefficient
   solution is done; if ncoeff=4, tr[0] = tr[3] and tr[4] = -tr[1]
   but the result is still returned as 6 coefficients.

   Inputs are supplied as generic pointers, which are assumed to 
   be arrays of elements with size given by the "sz" parameters.
   Double precision values of the relevant quantities are assumed
   to be at offsets given by the "off" parameters (these can be
   determined using "offsetof" for C structure members) within
   each array element.  This form allows quite a lot of flexibility
   when choosing the storage layout for the inputs. */

int dplate (void *comxptr, size_t comxoff, size_t comxsz,
            void *comyptr, size_t comyoff, size_t comysz,
            void *refxptr, size_t refxoff, size_t refxsz,
            void *refyptr, size_t refyoff, size_t refysz,
            void *wtptr, size_t wtoff, size_t wtsz,
            int npt, int ncoeff, double *tr) {
  double xbar = 0.0, ybar = 0.0, xref = 0.0, yref = 0.0, sw = 0.0;
  double sx1sq = 0.0, sy1sq = 0.0;
  double sx1y1 = 0.0, sx1x2 = 0.0, sy1x2 = 0.0;
  double sy1y2 = 0.0, sx1y2 = 0.0;
  double wt, xx1, yy1, xx2, yy2;
  double denom;

  int i, rv = 0;

#define DDREF(p) *((double *) (p))
#define COM_X(i) DDREF(comxptr + comxsz*(i) + comxoff)
#define COM_Y(i) DDREF(comyptr + comysz*(i) + comyoff)
#define REF_X(i) DDREF(refxptr + refxsz*(i) + refxoff)
#define REF_Y(i) DDREF(refyptr + refysz*(i) + refyoff)
#define GET_W(i) DDREF(wtptr  + wtsz*(i)  + wtoff)

  /* Compute mean x and y */
  for(i = 0; i < npt; i++) {
    wt = wtptr ? GET_W(i) : 1.0;

    if(wt) {
      xbar += COM_X(i) * wt;
      ybar += COM_Y(i) * wt;
      xref += REF_X(i) * wt;
      yref += REF_Y(i) * wt;
      sw += wt;
    }
  }

  if(sw > 0) {
    xbar /= sw;
    ybar /= sw;
    xref /= sw;
    yref /= sw;
  }
  else {
    /* No points = cannot proceed */
    tr[0] = 1.0;
    tr[1] = 0.0;
    tr[2] = 0.0;
    tr[3] = 1.0;
    tr[4] = 0.0;
    tr[5] = 0.0;

    return(-1);
  }

  /* Accumulate sums in zero mean system */
  for(i = 0; i < npt; i++) {
    wt = wtptr ? GET_W(i) : 1.0;

    if(wt) {
      xx1 = COM_X(i) - xbar;
      xx2 = REF_X(i) - xref;
      yy1 = COM_Y(i) - ybar;
      yy2 = REF_Y(i) - yref;
      
      sx1sq += xx1*xx1 * wt;
      sy1sq += yy1*yy1 * wt;
      sx1y1 += xx1*yy1 * wt;
      sx1x2 += xx1*xx2 * wt;
      sy1x2 += yy1*xx2 * wt;
      sy1y2 += yy1*yy2 * wt;
      sx1y2 += xx1*yy2 * wt;
    }
  }

  /* What solution do they want? */
  if(ncoeff == 4) {  /* 4-coeff */
    denom = sx1sq + sy1sq;

    if(denom > 0) {
      denom = 1.0 / denom;

      tr[0] = (sx1x2+sy1y2) * denom;
      tr[1] = (sy1x2-sx1y2) * denom;
      tr[3] = tr[0];
      tr[4] = -tr[1];
    }
    else {
      tr[0] = 1.0;
      tr[1] = 0.0;
      tr[3] = 1.0;
      tr[4] = 0.0;

      rv = -2;
    }
  }
  else {  /* 6-coeff */
    denom = sx1sq*sy1sq - sx1y1*sx1y1;

    if(denom != 0) {
      denom = 1.0 / denom;

      tr[0] = (sx1x2*sy1sq - sy1x2*sx1y1) * denom;
      tr[1] = (sy1x2*sx1sq - sx1x2*sx1y1) * denom;
      tr[3] = (sy1y2*sx1sq - sx1y2*sx1y1) * denom;
      tr[4] = (sx1y2*sy1sq - sy1y2*sx1y1) * denom;
    }
    else {
      tr[0] = 1.0;
      tr[1] = 0.0;
      tr[3] = 1.0;
      tr[4] = 0.0;

      rv = -2;
    }
  }

  tr[2] = xref - xbar*tr[0] - ybar*tr[1];
  tr[5] = yref - xbar*tr[4] - ybar*tr[3];

  return(rv);
}

