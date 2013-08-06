#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "astro.h"
#include "cvtunit.h"
#include "util.h"

/* Barycentering routines dealing with time and radial velocity.
   Follows the IAU 2000 resolutions, as described and elaborated
   in Lindegren & Dravins (2003).  Uses the "observer" structures
   for Earth attitude.

   General notes on accuracy: these routines are intended (in
   conjunction with the time ephemeris support and JPL ephemeris
   interpolation) to achieve microsecond accuracy.  However, this
   depends on the availability of an accurate TT - i.e. on clock
   corrections to UTC and TT(BIPM)-TT(TAI).  If using the simple
   TAI+DTT treatment, the results should be good to about 30 us,
   the rough size of TT(BIPM)-TT(TAI).

   The results have been tested against other publicly available
   implementations.  The results agree at the few tens of ns level
   with Eastman et al. (2010).  The remaining error was tracked
   down to three sources: their value of the AU is different; the
   JPL ephemeris interpolation uses full JD, rather than MJD, and
   thus loses a small amount of precision; and the Fairhead &
   Bretagnon (1990) approximation to TDB is only good to a few
   tens of ns. */

/* Returns time difference SSB - observer (seconds), i.e. add this
   to your observed TDB to convert it to BJD.  Source vector should
   be BCRS, i.e. corrected for space motion. */

double bary_delay (struct observer *obs, double s[3], double pr) {
  double vt[3];
  double dr, ds;

  /* Romer delay for infinite source distance: light path length
     between observer and SSB in AU. */
  dr = v_d_v(s, obs->bop);

  /* First order correction for wavefront curvature if source
     parallax is known.  */
  if(pr > 0) {
    v_x_v(s, obs->bop, vt);
    dr -= 0.5 * pr * v_d_v(vt, vt);
  }

  /* Shapiro delay due to the Sun.  Planets are neglected. */
  ds = 2.0 * log(1.0 + v_d_v(s, obs->hop) / obs->hdist) * GMSUN / (LIGHT * LIGHT);

#if 0
  printf("%.7e %.7e %.7e\n", obs->itgop[0], obs->itgop[1], obs->itgop[2]);
  printf("%.7e %.7e %.7e\n", obs->gop[0], obs->gop[1], obs->gop[2]);

  printf("%25.15f\n", obs->tt+obs->dtdb/DAY);
  printf("%.15f %.15f %.15f\n", obs->bep[0], obs->bep[1], obs->bep[2]);
  printf("%.15f %.15f %.15f\n", s[0], s[1], s[2]);

  printf("%.10f %.10f %.10f %.10f\n",
	 v_d_v(s, obs->bep)*AU/LIGHT,
	 v_d_v(s, obs->gop)*AU/LIGHT,
	 dr*AU/LIGHT, ds/LIGHT);
#endif

  return((dr * AU + ds) / LIGHT);
}

/* Returns (1 + z_b) / (1 + z_obs) - 1.  The IAU recommendation
   is to report c z_b as "Barycentric radial velocity measure". */

double bary_doppler (struct observer *obs, double s[3]) {
  /* Taylor expansion of Eq. 41 in the paper, keeping only
     first order terms.  Next order corrections are < 1mm/s. */
  return(obs->gab * v_d_v(s, obs->vab) * (1.0 + obs->gpearth + obs->gpsun));
}
