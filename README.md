lib
===

Library of miscellaneous astrophysics-related (and some not) C
routines.  Much of my other software depends on it.  The static
library includes all routines that can be built on a standard C99 
compiler without extra library dependencies.  There are a few other
extra routines that do have external dependencies and are not included
in the archive.  These are organized into one file for each external
library.

This is still under active development, and has been for some time.
It has got a bit messy and could do with a tidy.

Some example programs (test*.c) for the astrometric routines are
included.  These need the JPL ephemerides and IERS tables, which are
located using environment variables.  Documentation (including
installation instructions) is still quite lacking, but I'm working on
it, the main reason this is up here in such an unfinished state is for
my convenience.
