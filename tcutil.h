#ifndef TCUTIL_H
#define TCUTIL_H

int tcutil_init (void);
void tcutil_winsize (int *ncols, int *nrows);
void tcutil_goto (int x, int y);
void tcutil_back (int n);
void tcutil_up (int n);

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
