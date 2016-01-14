#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include <fitsio.h>

#include "ap.h"
#include "lfa.h"
#include "util.h"
#include "fitsutil.h"

static int do_file (char *infile, char *maskfile, char *outfile, int nbsize,
                    char *errstr);

static void usage (char *av) {
  fprintf(stderr, "Usage:\t%s infile maskfile|nomask outfile nbsize\n", av);
  exit(1);
}

int main (int argc, char *argv[]) {
  char *pn = (char *) NULL, errstr[ERRSTR_LEN];
  char *ep;
  int nbsize;

  /* Set the program name for error reporting */
  if(argv[0])
    pn = basename(argv[0]);
  
  if(!pn)
    pn = "testbackground";

  setprogname(pn);

  if(argc != 5)
    usage(argv[0]);

  nbsize = strtol(argv[4], &ep, 10);
  if(*ep || nbsize <= 1)
    fatal(1, "invalid nbsize: %s", argv[4]);

  if(do_file(argv[1], argv[2], argv[3], nbsize, errstr))
      fatal(1, "%s", errstr);

  return(0);
}

static int do_file (char *infile, char *maskfile, char *outfile, int nbsize,
                    char *errstr) {
  fitsfile *inf, *mskf, *outf;
  char fnbuf[FLEN_FILENAME];

  int status = 0;

  int bitpix, naxis, npix;
  long naxes[2];

  int msk_bitpix, msk_naxis, msk_npix;
  long msk_naxes[2];

  int i, nkeys, kclass;
  char card[FLEN_CARD];

  float pedestal;
  int p;

  float *buf = (float *) NULL;
  unsigned char *mask = (unsigned char *) NULL;

  float skylev, skynoise;

  /* Open input file */
  ffopen(&inf, infile, READONLY, &status);
  if(status) {
    fitsio_err(errstr, status, "ffopen: %s", infile);
    goto error;
  }

  /* Create output file */
  fnbuf[0] = '!';
  strncpy(&(fnbuf[1]), outfile, sizeof(fnbuf)-2);
  fnbuf[sizeof(fnbuf)-1] = '\0';

  ffinit(&outf, fnbuf, &status);
  if(status) {
    fitsio_err(errstr, status, "ffinit: %s", outfile);
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

  /* Create output HDU, change to float */
  ffcrim(outf, FLOAT_IMG, naxis, naxes, &status);
  if(status) {
    fitsio_err(errstr, status, "ffcrim");
    goto error;
  }

  /* Copy rest of header */
  ffghsp(inf, &nkeys, (int *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "ffghsp");
    goto error;
  }

  for(i = 1; i <= nkeys; i++) {
    ffgrec(inf, i, card, &status);
    if(status) {
      fitsio_err(errstr, status, "ffgrec: card %d", i);
      goto error;
    }
    
    kclass = ffgkcl(card);
    switch(kclass) {
    case TYP_STRUC_KEY:
    case TYP_CMPRS_KEY:
    case TYP_CKSUM_KEY:
    case TYP_SCAL_KEY:
      break;  /* Skip */
    default:
      if(strncmp(card, "PEDESTAL", 8)) {  /* remove PEDESTAL if there */
        ffprec(outf, card, &status);
        if(status) {
          fitsio_err(errstr, status, "could not write key %d", i);
          goto error;
        }
      }
    }
  }

  /* Read PEDESTAL */
  ffgkye(inf, "PEDESTAL", &pedestal, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    pedestal = 0.0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye");
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

  /* Fix PEDESTAL if present */
  for(p = 0; p < npix; p++)
    buf[p] -= pedestal;

  /* Close file */
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

  /* Compute sky level before */
  skylevel_image(buf, mask, npix,
                 -FLT_MAX, 3,
                 &skylev, &skynoise);

  fprintf(stderr, "Sky level before = %.3f noise = %.3f\n",
          skylev, skynoise);

  /* Remove background */
  if(backremove(buf, mask, buf, naxes[0], naxes[1], nbsize)) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* Compute sky level after */
  skylevel_image(buf, mask, npix,
                 -FLT_MAX, 3,
                 &skylev, &skynoise);

  fprintf(stderr, "Sky level after = %.3f noise = %.3f\n",
          skylev, skynoise);

  /* Write map */
  ffppre(outf, 1, 1, npix, buf, &status);
  if(status) {
    fitsio_err(errstr, status, "ffppre");
    goto error;
  }

  free((void *) buf);
  buf = (float *) NULL;

  if(mask) {
    free((void *) mask);
    mask = (unsigned char *) NULL;
  }
  
  ffclos(outf, &status);
  if(status) {
    fitsio_err(errstr, status, "ffclos");
    goto error;
  }

  return(0);
  
 error:
  if(buf)
    free((void *) buf);
  if(mask)
    free((void *) mask);

  return(1);
}
