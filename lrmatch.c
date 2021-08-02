#include <stdlib.h>
#include <math.h>

#include "lfa.h"

static int bserchmult (double coma, double comb,
                       double *refa, double *refb, int nref,
                       double errlim,
                       int *refrows, double *refsepsq) {
  int isp, ifp, index;
  double errsq, da, db, dsq;
  int nmatch;
  
  isp = 1;
  ifp = nref;
  errsq = errlim*errlim;
  index = (isp+ifp)/2;
  
  while(ifp-isp >= 2) {
    if(refa[index-1] < coma - errlim) {
      isp = index;
      index = (index+ifp)/2;
    }
    else if(refa[index-1] > coma - errlim) {
      ifp = index;
      index = (index+isp)/2;
    }
    else {
      isp = index;
      break;
    }
  }

  nmatch = 0;
  
  for(index = isp; index <= nref; index++) {
    if(refa[index-1] > coma + errlim)
      break;

    da = coma - refa[index-1];
    db = comb - refb[index-1];
    dsq = da*da + db*db;
    if(dsq < errsq) {
      refrows[nmatch] = index-1;
      refsepsq[nmatch] = dsq;
      nmatch++;
    }
  }

  return(nmatch);
}

int lrmatch (double *comx, double *comy,
             double *comlogrank, double *comerr, int ncom,
             double *refx, double *refy,
             double *reflogrank, double *referr, int nref,
             double searchrad, int sorted_y,
             int *best_ref_for_com, int *best_com_for_ref) {
  double *best_lr_for_ref = (double *) NULL;

  int *refrows = (int *) NULL;
  double *refsepsq = (double *) NULL;

  int icom, iref, imatch, nmatch;

  double sepsq, rsq, lr;
  
  int best_iref, prev_icom;
  double best_lr;
  
  /* Allocate workspace */
  best_lr_for_ref = (double *) malloc(nref * sizeof(double));
  refrows = (int *) malloc(nref * sizeof(int));
  refsepsq = (double *) malloc(nref * sizeof(double));
  if(!best_lr_for_ref || !refrows || !refsepsq)
    goto error;
  
  /* Initialize output */
  for(icom = 0; icom < ncom; icom++)
    best_ref_for_com[icom] = -1;
  
  for(iref = 0; iref < nref; iref++)
    best_com_for_ref[iref] = -1;

  /* For each comparison object, search for reference objects */
  for(icom = 0; icom < ncom; icom++) {
    if(sorted_y)
      nmatch = bserchmult(comy[icom], comx[icom], refy, refx, nref,
                          searchrad, refrows, refsepsq);
    else
      nmatch = bserchmult(comx[icom], comy[icom], refx, refy, nref,
                          searchrad, refrows, refsepsq);

    if(nmatch > 0) {
      best_iref = -1;
      best_lr = 0;

      for(imatch = 0; imatch < nmatch; imatch++) {
        iref = refrows[imatch];
        sepsq = refsepsq[imatch];

        rsq = sepsq;
        if(comerr) {
          if(referr)
            rsq /= (comerr[icom]*comerr[icom] + referr[iref]*referr[iref]);
          else
            rsq /= (comerr[icom]*comerr[icom]);
        }
        else {
          if(referr)
            rsq /= (referr[iref]*referr[iref]);
        }
        
        lr = -0.5 * rsq;
        if(comlogrank)
          lr -= comlogrank[icom];
        if(reflogrank)
          lr -= reflogrank[iref];

        if(best_iref < 0 || lr > best_lr) {
          best_iref = iref;
          best_lr = lr;
        }
      }

      prev_icom = best_com_for_ref[best_iref];
      if(prev_icom < 0 || best_lr >= best_lr_for_ref[best_iref]) {
        best_com_for_ref[best_iref] = icom;
        best_lr_for_ref[best_iref] = best_lr;

        if(prev_icom >= 0)
          best_ref_for_com[prev_icom] = -1;
        
        best_ref_for_com[icom] = best_iref;
      }
    }
  }

  free((void *) best_lr_for_ref);
  best_lr_for_ref = (double *) NULL;
  free((void *) refrows);
  refrows = (int *) NULL;
  free((void *) refsepsq);
  refsepsq = (double *) NULL;
  
  return(0);

 error:
  if(best_lr_for_ref)
    free((void *) best_lr_for_ref);
  if(refrows)
    free((void *) refrows);
  if(refsepsq)
    free((void *) refsepsq);

  return(-1);
}
