#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "astro.h"
#include "cvtunit.h"
#include "util.h"

/* Creates source structure from star coordinates */

void source_star (struct source *src,
		  double ra, double de,      /* radians */
		  double pmra, double pmde,  /* sky projected, arcsec/yr */
		  double plx, double vrad,   /* arcsec, km/s */
		  double epoch) {
  double sa, sd, ca, cd, df, vf;

  /* Set source type */
  src->type = SOURCE_STAR;

  /* Precompute these */
  dsincos(ra, &sa, &ca);
  dsincos(de, &sd, &cd);

  /* Unit vector to star at catalogue epoch */
  src->ref_n[0] = ca*cd;
  src->ref_n[1] = sa*cd;
  src->ref_n[2] = sd;

  /* Convert from catalogue units to radians, rad/day, m/s */
  pmra *= AS_TO_RAD / JYR;  /* da/dt * cos(delta) */
  pmde *= AS_TO_RAD / JYR;  /* dd/dt */
  vrad *= 1000;

  src->pr = plx * AS_TO_RAD;  /* parallax */

  /* Space motion: if star is at x = r * n = n / pr, then
     dx/dt = dx/da * da/dt + dx/dd * dd/dt + dx/dr * dr/dt
           = r dn/da * da/dt + r dn/dd * dd/dt + n * dr/dt
     so
     dn/dt = dn/da * da/dt + dn/dd * dd/dt + n * pr * dr/dt
  */

  /* Doppler shift factor in apparent velocities, accounts for changing
     light travel time to star.  e.g. Klioner 2003. */
  df = 1.0 / (1.0 - vrad / LIGHT);

  /* Rate of change of vector due to radial velocity */
  vf = src->pr * vrad * DAY / AU;   /* pr * dr/dt */

  /* Space motion / day */
  src->ref_dndt[0] = df * (-pmra * sa - pmde * ca*sd + vf * src->ref_n[0]);
  src->ref_dndt[1] = df * ( pmra * ca - pmde * sa*sd + vf * src->ref_n[1]);
  src->ref_dndt[2] = df * (             pmde *    cd + vf * src->ref_n[2]);

  /* Stash epoch */
  src->ref_tdb = J2K + (epoch - 2000.0) * JYR;
}

/* Creates source structure from Keplerian osculating elements */

void source_elem (struct source *src,
		  unsigned char eltype,
		  double epoch,
		  double incl,
		  double anode,
		  double longperi,
		  double aq,
		  double ecc,
		  double lm,
		  double nn) {
  double aa;

  /* Set source type */
  src->type = eltype;

  /* Generate rotation matrix */
  m_identity(src->rot);
  euler_rotate(src->rot, 3, -longperi);
  euler_rotate(src->rot, 1, -incl);
  euler_rotate(src->rot, 3, -anode);

  /* Stash the rest */
  src->ref_tdb = J2K + (epoch - 2000.0) * JYR;
  src->aq = aq;
  src->ecc = ecc;

  switch(eltype) {
  case SOURCE_ELEM_MAJOR:
    /* Mean anomaly from mean longitude */
    src->ma = lm - anode - longperi;

    /* Mean motion */
    src->nn = nn;

    break;
  case SOURCE_ELEM_MINOR:
    /* Mean anomaly */
    src->ma = lm;

    /* Compute mean motion */
    aa = aq * AU;
    src->nn = sqrt(GMSUN / (aa*aa*aa)) / DAY;

    break;
  case SOURCE_ELEM_COMET:
    /* Everything is relative to perihelion for comets */

    break;
  }
}

void source_place (struct observer *obs,
		   struct source *src,
		   unsigned char mask,    /* parts we want */
		   double *s,             /* result */
		   double *dsdt,          /* src motion only */
		   double *pr_r) {
  double vb[3], tmp[3];
  double pr, nf, dt;
  int i;

  if(src->type == SOURCE_STAR) {
    v_copy(s, src->ref_n);
    pr = src->pr;

    /* Space motion */
    if(mask & TR_MOTION) {
      /* Time difference */
      dt = obs->tdb - src->ref_tdb;

      /* Add in Romer delay */
      dt += v_d_v(s, obs->bop) * AU / (LIGHT*DAY);

      v_p_sv(s, dt, src->ref_dndt);
    }

    /* Parallax */
    if(mask & TR_PLX) {
      v_p_sv(s, -pr, obs->bop);
    }
    
    /* Renormalize */
    nf = v_renorm(s);

    pr *= nf;

    if(dsdt) {
      if(mask & TR_PLX) {
	for(i = 0; i < 3; i++)
	  dsdt[i] = nf*src->ref_dndt[i] - pr*(obs->bev[i] + obs->gov[i]);
      }
      else {
	for(i = 0; i < 3; i++)
	  dsdt[i] = nf*src->ref_dndt[i];
      }
    }
  }
  else {
    if(src->type < SOURCE_STAR) {
      /* From JPL, first figure out if we already got it.  The vector
	 "vb" is the motion of the source relative to the SSB only. */
      if(src->type == JPLEPH_SUN) {
	if(mask & TR_PLX) {
	  for(i = 0; i < 3; i++) {
	    s[i]  = -obs->hop[i];
	    vb[i] = -obs->bsv[i];
	  }
	}
	else {
	  for(i = 0; i < 3; i++) {
	    s[i]  = -obs->bsp[i];
	    vb[i] = -obs->bsv[i];
	  }
	}
      }
      else if(src->type == JPLEPH_MOON) {
	if(mask & TR_PLX) {
	  for(i = 0; i < 3; i++) {
	    s[i]  = obs->emp[i] - obs->gop[i];
	    vb[i] = obs->bev[i] + obs->emv[i];
	  }
	}
	else {
	  for(i = 0; i < 3; i++) {
	    s[i]  = obs->bep[i] + obs->emp[i];
	    vb[i] = obs->bev[i] + obs->emv[i];
	  }
	}
      }
      else {
	/* Fetch from ephemeris.  It should already be loaded and
	   interpolated in observer_update, so this call cannot fail. */
	(void) jpleph_fetch(obs->jpltab, obs->tdb, src->type, s, vb);

	if(mask & TR_PLX) {
	  for(i = 0; i < 3; i++)
	    s[i] -= obs->bop[i];
	}
      }
    }
    else {  /* Orbital elements */

      /* Move to separate orbit subroutine, generalise for e = 1 and e > 1 */

#if 0
      /* Mean anomaly, restricted to [0, TWOPI) */
      ma = fmod(src->ma + src->nn * (obs->tdb - src->ref_tdb), TWOPI);
      if(ma < 0)
	ma += TWOPI;

      /* Eccentric anomaly: solve Kepler's equation */
      ea = ma;

      for(iter = 0; iter < 3; iter++) {
	f = ea - src->ecc * sin(ea) - ma;
	df = 1.0 - src->ecc * cos(ea);
	delta = f / df;  /* KABOOM? */
	ea -= delta;
    
	if(fabs(delta) < 1.0e-8) {
	  /* I think that's enough... */
	  break;
	}
      }

      /* True anomaly */
      
#endif

      /* Orbital plane to ecliptic, then ecliptic to GCRS */
      m_x_v(src->rot, s, tmp);
      mt_x_v(obs->eclm, tmp, s);

      m_x_v(src->rot, vb, tmp);
      mt_x_v(obs->eclm, tmp, vb);
    }

    /* Romer delay */
    dt = -sqrt(v_d_v(s, s)) * AU / (DAY*LIGHT);

    /* Apply first order correction for light travel time from source
       ("planetary aberration").  There is no need to iterate at the
       desired mas accuracy level. ?? */
    for(i = 0; i < 3; i++)
      s[i] += dt * vb[i];

    /* Normalize */
    pr = v_renorm(s);

    /* Form full time derivative of source vector */
    if(dsdt) {
      if(mask & TR_PLX) {
	for(i = 0; i < 3; i++)
	  dsdt[i] = pr*(vb[i] - obs->bev[i] - obs->gov[i]);
      }
      else {
	for(i = 0; i < 3; i++)
	  dsdt[i] = pr*vb[i];
      }
    }
  }

  if(mask > (TR_MOTION | TR_PLX))
    observer_ast2obs(obs, s, pr, mask);

  if(pr_r)
    *pr_r = pr;
}

