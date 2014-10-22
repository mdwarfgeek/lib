#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lfa.h"

/* Barycentering routines dealing with time and radial velocity.
   Follows the IAU 2000 resolutions, as described and elaborated
   in Lindegren & Dravins (2003).  Uses the "observer" structures
   for Earth attitude. */

/* Returns time difference SSB - observer (seconds), i.e. add this
   to your observed TDB to convert it to BJD.  Source vector should
   be BCRS, i.e. corrected for space motion.

   General notes on accuracy: this routine is intended (in
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

double bary_delay (struct observer *obs, double s[3], double pr) {
  double p[3], vt[3];
  double dr, ds, stmp;

  /* Romer delay for infinite source distance: light path length
     between observer and SSB in AU. */
  dr = v_d_v(s, obs->bop);

  if(pr > 0) {
    /* First order correction for wavefront curvature if source
       parallax is known.  */
    v_x_v(s, obs->bop, vt);
    dr -= 0.5 * pr * v_d_v(vt, vt);

    /* Compute normalized observer to source vector */
    memcpy(p, s, sizeof(p));
    v_p_sv(p, -pr, obs->bop);
    v_renorm(p);

    /* Scalar product in Shapiro delay */
    stmp = v_d_v(p, obs->hop);
  }
  else {
    stmp = v_d_v(s, obs->hop);
  }

  /* Shapiro delay due to the Sun.  Planets are neglected. */
  ds = 2.0 * log(1.0 + stmp / obs->hdist) * GMSUN / (LIGHT * LIGHT);

  return((dr * AU + ds) / LIGHT);
}

/* Returns (1 + z_b) / (1 + z_obs) - 1.  The IAU recommendation
   is to report c z_b as "Barycentric radial velocity measure".

   Source vectors s and dsdt should be corrected for space motion.
   The source vector sref is the source position at the catalogue
   epoch.  s and sref should both have unit norm.

   Calculations are based on:
   Wright & Eastman (2014, arXiv:1409.4774v1).
   Compared to this work, we use a slightly more rigorous treatment
   of the Shapiro delay (but only for the Sun, they also include
   the planets) and light travel.  The latter is already included
   in the vectors s and dsdt inside source.c rather than appearing
   here as an approximate term. */

double bary_doppler (struct observer *obs, double sref[3],
                     double s[3], double dsdt[3], double pr) {
  double p[3], dpdt[3], nf;
  double zgr, zppm0, zppm, zsd;
  int i;

  /* Compute normalized observer to source vector (p) and time
     derivative (dpdt). */
  memcpy(p, s, sizeof(p));
  v_p_sv(p, -pr, obs->bop);
  nf = v_renorm(p);

  memcpy(dpdt, dsdt, sizeof(dpdt));
  v_p_sv(dpdt, -pr, obs->bev);
  v_p_sv(dpdt, -pr, obs->gov);
  
  for(i = 0; i < 3; i++)
    dpdt[i] *= nf;

  /* Lorentz transformations of frequency (with GR terms added):

     nu_ssb = nu_emit         (1 + z_gr*)
                       ------------------------
                       gamma_* (1 + beta_* . s)

     at the SSB, and

     nu_meas = nu_ssbo  gamma_o (1 + beta_o . p)
                        ------------------------
                         (1 + z_gr) (1 + z_sd)

     at the observer, where nu_ssbo denotes frequency measured at
     the velocity of the SSB, but the position of the observer.
     The quantity nu_ssbo can be obtained simply by replacing "s"
     with "p" in the nu_ssb equation.  The GR terms at the SSB
     are omitted by convention - i.e. the "Barycentric frequency"
     we want is reckoned as if there was no gravity there.

     For this purpose, we're interested only in the time-variable
     terms as measured at the SSB, and we can only measure
     frequency (we don't usually know what z_gr* is).  Therefore
     the desired quantity for the star is nu_ssbo / nu_ssb0 at
     some reference epoch.  Assuming beta_* is (approximately)
     constant, then:

     nu_ssbo = nu_ssb0  1 + beta_* . s0
                        ---------------
                        1 + beta_* . p

     and therefore:

     nu_meas = nu_ssb0  gamma_o (1 + beta_o . p) (1 + beta_* . s0)
                        ------------------------------------------
                          (1 + z_gr) (1 + z_sd) (1 + beta_* . p)

     Result is (1+z_b) = nu_meas / nu_cat */

  /* beta_star . s and p terms */
  zppm0 = v_d_v(dsdt, sref) * AU / (DAY * LIGHT * pr);
  zppm  = v_d_v(dsdt, p) * AU / (DAY * LIGHT * pr);

  /* z_gr, sum of gravitational redshift terms.
     Sun and Earth only, planets can be neglected at 1 mm/s level. */
  zgr = -obs->gpearth - obs->gpsun;

  /* z_sd, time derivative of the Shapiro delay due to the Sun.
     Planets are neglected. */
  zsd = -2.0 * ((v_d_v(p, obs->hev) +
                 v_d_v(p, obs->gov) +
                 v_d_v(dpdt, obs->hop)) /
                (obs->hdist + v_d_v(p, obs->hop)))
      * GMSUN / (DAY * LIGHT * LIGHT * LIGHT);

  return(obs->gab * (1.0 + v_d_v(p, obs->vab)) * (1.0 + zppm0)
                  / ((1.0 + zgr) * (1.0 + zsd) * (1.0 + zppm)) - 1.0);
}
