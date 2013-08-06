#ifndef TCUTIL_H
#define TCUTIL_H

int tcutil_init (void);
void tcutil_winsize (int *ncols, int *nrows);
void tcutil_back (void);
void tcutil_up (void);

/* Simplified attribute interface */
enum {
  ATTR_NORM,
  ATTR_BLINK,
  ATTR_BOLD,
  ATTR_HALF,
  ATTR_REVERSE,
  ATTR_STANDOUT,
  ATTR_UNDERLINE,
  NUM_ATTR
};

void tcutil_attr (int a);

#endif  /* TCUTIL_H */
