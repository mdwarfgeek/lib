#ifndef FIND_H
#define FIND_H

struct specfind_line {
  /* 1st moment, centre of first pixel is x=1 */
  double x;
  
  /* Isophotal flux */
  double ziso;

  /* Peak height above sky */
  double zpk;

  /* Bounding box */
  int xl;
  int xh;

  /* Flags */
  unsigned char flags;
};

int specfind (float *line, unsigned char *mask, int nx,
              int minpix, float sky, float thresh,
              struct specfind_line **list, int *nlist);

#endif  /* FIND_H */
