#ifndef AP_H
#define AP_H

/* Number of areal profile levels */
#define AP_NAREAL 8

struct ap_parent {
  int first;
  int last;
  int refcount;
  int npixels;
};

struct ap_pixel {
  int next;
  int xa;
  int xb;
  int y;
};

struct ap_source {
  /* 1st moments, centre of first pixel is 1,1 as imcore */
  double x;
  double y;

  /* 2nd moments */
  double sxx;
  double syy;
  double sxy;

  /* Isophotal flux */
  double ziso;

  /* Peak height above sky */
  double zpk;

  /* Areal profiles at 1, 2, 4, ... times threshold above sky */
  int areal[AP_NAREAL];

  /* Bounding box */
  int bb[4];

  /* Flags */
  unsigned char flags;
#define AP_SOURCE_FLAG_BADPIX 0x01  /* contains bad pixels */
};

struct ap {
  int nx;
  int ny;

  float *line;

  struct ap_parent *parents;
  int *free_parents;
  int nalloc_parents;
  int nfree_parents;

  struct ap_pixel *pixels;
  int *free_pixels;
  int nalloc_pixels;
  int nfree_pixels;

  int *line_parents;

  struct ap_source *sources;
  int nsources;
  int nalloc_sources;

  int (*output) (struct ap *ap, struct ap_source *obj, void *user_data);
  void *output_user_data;
};

/* A very simple object detector for autoguiding, based on the APM
   algorithm (see Irwin 1985).  Intended for repeated use on a large
   number of images of the same size and similar source density, so
   re-uses memory to avoid malloc overhead at each call.  There is
   no deblending.

   Image arrays:

   map      The image itself, with sky left on.

   filtmap  Optional: image map convolved with matched detection
                      filter.  If given, this argument is used rather
                      than applying the internal (Hanning) matched
                      filter.  The map needs to be multiplied by
                      the bad pixel mask before convolution, if using
                      one.

   mask     Optional: boolean bad pixel mask, 1 if good, 0 if bad

   Parameters:

   minpix   Minimum number of pixels above threshold for a detection.
   sky      Sky level (e.g. as computed by "skylevel").
   thresh   Threshold in counts above sky.

   The threshold is usually set to some constant value times the
   sky noise, e.g. as in imcore_conf.

   By default, parameters of the detected sources are stored in the
   structure member "sources" and can be accessed directly after a
   successful call.  This memory is re-used at the next call.

   Optionally, by setting structure members "output" and
   "output_user_data", a callback function can be used for output
   rather than storing to memory.  This might be used to write the
   source list directly to a file rather than storing it in RAM. */

int ap_init (struct ap *ap, int nx, int ny);
void ap_free (struct ap *ap);
int ap_image (struct ap *ap,
              float *map, float *filtmap, unsigned char *mask,
              int minpix, float sky, float thresh);

/* Ellipse parameters from second moments */
void ap_ellipse (struct ap_source *obj,
                 double *gfwhm, double *ell, double *spa, double *cpa);

#endif  /* AP_H */
