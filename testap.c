#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fitsio.h>

#include "ap.h"
#include "lfa.h"
#include "util.h"
#include "fitsutil.h"

/* testap: a simple test program for the object detector, produces a
   DS9 ellipse overlay file on stdout from a FITS image.

   Takes the same first four arguments as imcore_conf, except we use
   boolean bad pixel masks here rather than confidence maps. */

static int do_file (char *filename, char *maskfile,
                    int minpix, float thresh, int nbsize,
                    char *errstr);

static void usage (char *av) {
  fprintf(stderr,
          "Usage:\t%s infile maskfile|nomask minpix thresh [nbsize]\n",
          av);
  exit(1);
}

int main (int argc, char *argv[]) {
  char *pn = (char *) NULL, errstr[ERRSTR_LEN];
  char *ep;
  int minpix;
  float thresh;

  int nbsize = 64;

  /* Set the program name for error reporting */
  if(argv[0])
    pn = basename(argv[0]);
  
  if(!pn)
    pn = "testap";

  setprogname(pn);

  if(argc < 5)
    usage(argv[0]);

  minpix = strtol(argv[3], &ep, 10);
  if(*ep)
    fatal(1, "invalid minpix: %s", argv[3]);

  thresh = strtod(argv[4], &ep);
  if(*ep)
    fatal(1, "invalid threshold: %s", argv[4]);

  if(argc > 5) {
    nbsize = strtol(argv[5], &ep, 10);
    if(*ep || nbsize < 1)
      fatal(1, "invalid nbsize: %s", argv[5]);
  }

  if(do_file(argv[1], argv[2], minpix, thresh, nbsize, errstr))
      fatal(1, "%s", errstr);

  return(0);
}

static int test_output (struct ap *ap, struct ap_source *obj, void *data) {
  FILE *fp;
  double gfwhm, ell, spa, cpa, axr, a, b, pa;

  fp = (FILE *) data;

  /* Produce DS9 ellipse overlay */
  ap_ellipse(obj, &gfwhm, &ell, &spa, &cpa);

  axr = 1.0 - ell;
  
  a = sqrt(obj->areal[0] / (M_PI * axr));
  b = a * axr;
  
  pa = atan2(spa, cpa);
  
  printf("image; ellipse %9.3f %9.3f %8.3f %8.3f %8.3f\n",
         obj->x,
         obj->y,
         a,
         b,
         pa * 180.0 / M_PI);

  return(0);
}

static int do_file (char *filename, char *maskfile,
                    int minpix, float thresh, int nbsize,
                    char *errstr) {
  fitsfile *inf, *mskf;
  int status = 0;

  int bitpix, naxis, npix;
  long naxes[2];

  int msk_bitpix, msk_naxis, msk_npix;
  long msk_naxes[2];

  float *buf = (float *) NULL;
  unsigned char *mask = (unsigned char *) NULL;

  float skylev, skynoise;

  struct ap ap;

  /* Open file */
  ffopen(&inf, filename, READONLY, &status);
  if(status) {
    fitsio_err(errstr, status, "ffopen: %s", filename);
    goto error;
  }

  /* Get parameters */
  ffgipr(inf, 2, &bitpix, &naxis, naxes, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgipr");
    goto error;
  }

  if(naxis < 2) {
    report_err(errstr, "image is less than 2-dimensional!");
    goto error;
  }

  /* Read data */
  npix = naxes[0] * naxes[1];
  
  buf = (float *) malloc(npix * sizeof(float));
  if(!buf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  ffgpve(inf, 1, 1, npix, 0, buf, (int *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgpve");
    goto error;
  }

  ffclos(inf, &status);
  if(status) {
    fitsio_err(errstr, status, "ffclos");
    goto error;
  }

  /* Handle mask, if given */
  if(maskfile && strcmp(maskfile, "nomask")) {
    /* Open file */
    ffopen(&mskf, maskfile, READONLY, &status);
    if(status) {
      fitsio_err(errstr, status, "ffopen: %s", maskfile);
      goto error;
    }

    /* Get parameters */
    ffgipr(mskf, 2, &msk_bitpix, &msk_naxis, msk_naxes, &status);
    if(status) {
      fitsio_err(errstr, status, "ffgipr");
      goto error;
    }

    if(msk_naxis < 2) {
      report_err(errstr, "mask is less than 2-dimensional!");
      goto error;
    }
    
    msk_npix = msk_naxes[0]*msk_naxes[1];

    if(msk_npix < npix) {
      report_err(errstr,
                 "%ldx%ld mask smaller than %ldx%ld image",
                 msk_naxes[0], msk_naxes[1],
                 naxes[0], naxes[1]);
      goto error;
    }

    /* Read data */
    mask = (unsigned char *) malloc(npix * sizeof(unsigned char));
    if(!mask) {
      report_syserr(errstr, "malloc");
      goto error;
    }
    
    ffgpvb(mskf, 1, 1, npix, 0, mask, (int *) NULL, &status);
    if(status) {
      fitsio_err(errstr, status, "ffgpvb");
      goto error;
    }

    ffclos(mskf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffclos");
      goto error;
    }
  }

  /* Take off background variations.  Ignores mask for now due to
     issues if whole background bins are masked. */
  if(backremove(buf, NULL, buf, naxes[0], naxes[1], nbsize)) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* Compute sky level of flattened map */
  skylevel_image(buf, mask, npix,
                 -FLT_MAX, 3,
                 &skylev, &skynoise);

  fprintf(stderr, "Sky level = %.3f noise = %.3f threshold = %.3f\n",
          skylev, skynoise, thresh*skynoise);

  /* Detect images */
  if(ap_init(&ap, naxes[0], naxes[1])) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  ap.output = test_output;
  ap.output_user_data = stdout;

  if(ap_image(&ap, buf, NULL, mask, minpix, skylev, thresh*skynoise)) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  fprintf(stderr, "Detected %d sources\n", ap.nsources);

  ap_free(&ap);

  free((void *) buf);
  buf = (float *) NULL;

  if(mask) {
    free((void *) mask);
    mask = (unsigned char *) NULL;
  }

  return(0);
  
 error:
  if(buf)
    free((void *) buf);
  if(mask)
    free((void *) mask);

  return(1);
}
