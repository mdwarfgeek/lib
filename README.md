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
It has got a bit messy and could do with a tidy.  The Makefile
contains targets for a number of test programs, but these (and their
external depenencies) are not needed to build the library itself.

Some example programs (test*.c) for the astrometric routines are
included.  Many of the astrometric routines need the JPL ephemerides
and IERS tables, which are located using environment variables.

I recommend using the JPL DE430t ephemeris for most purposes, which
can be found in a binary form suitable for use with the library here:
ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/de430t/
Download the file linux_p1550p2650.430t and set JPLEPH_DATA to the
path to this file.

The IERS tables needed are finals2000A.data and tai-utc.dat.  They can
be obtained from: ftp://maia.usno.navy.mil/ser7/
Set IERS_DATA to the directory containing these files.  If running a
telescope or processing data regularly you may want to set up a cron
job to keep these files up to date.

Documentation is still quite lacking, but I'm working on it, the main
reason this is up here in such an unfinished state is because programs
I'm pushing to the repositories depend on the library (usually the more
finished parts).
