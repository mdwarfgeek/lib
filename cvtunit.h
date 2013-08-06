#ifndef __CVTUNIT_H__
#define __CVTUNIT_H__

/* --Useful conversion constants-- */

#define MAS_TO_RAD (M_PI / 648000000.0)
#define DEG_TO_RAD (M_PI / 180.0)
#define AM_TO_RAD (M_PI / 10800.0)
#define AS_TO_RAD (M_PI / 648000.0)
#define HR_TO_RAD (M_PI / 12.0)
#define SEC_TO_RAD (M_PI / 43200.0)

#define RAD_TO_DEG (180.0 / M_PI)
#define RAD_TO_AM  (10800.0 / M_PI)
#define RAD_TO_AS  (648000.0 / M_PI)
#define RAD_TO_HR  (12.0 / M_PI)

#define DEG_TO_AS 3600

/* --Units-- */

/* Radians */
#define UNIT_RAD  0

/* Hours, minutes, seconds, milliseconds */
#define UNIT_HR   1
#define UNIT_M    2
#define UNIT_S    3
#define UNIT_MS   4

/* Degrees and minutes, seconds and milliseconds of arc */
#define UNIT_DEG  5
#define UNIT_AM   6
#define UNIT_AS   7
#define UNIT_MAS  8

/* --Functions-- */

/* Convert base 60 to base 10 */
int base60_to_10 (char *, char **, char *, short, double *, short);

/* Convert base 10 to base 60 */
int base10_to_60 (double, short, char *, size_t,
		  char *, char *, short, short);

/* Return conversion factor between two units */
double get_cfactor (short, short);

/* Extract components of a date */
int parsedate (char *, char *, int *, int *, int *);

#endif  /* __CVTUNIT_H__ */
