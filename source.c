#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lfa.h"
#include "cvtunit.h"
#include "util.h"

#define KEPLER_PREC    1.0e-13
#define KEPLER_MAXITER 20

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
  inline_sincos(ra, sa, ca);
  inline_sincos(de, sd, cd);

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

void source_star_vec (struct source *src,
		      double *n,
		      double *dndt,
		      double pr,
		      double epoch) {

  /* Set source type */
  src->type = SOURCE_STAR;

  /* Copy in parameters */
  memcpy(&(src->ref_n[0]), n, sizeof(src->ref_n));
  if(dndt)
    memcpy(&(src->ref_dndt[0]), dndt, sizeof(src->ref_dndt));
  else {
    src->ref_dndt[0] = 0;
    src->ref_dndt[1] = 0;
    src->ref_dndt[2] = 0;
  }

  src->pr = pr;

  /* Stash epoch */
  src->ref_tdb = J2K + (epoch - 2000.0) * JYR;
}

/* Creates source structure from Keplerian osculating elements */

#ifdef TEST
void sla_el2ue_ (double *date, int *jform,
		 double *epoch, double *orbinc, double *anode, double *perih,
		 double *aorq, double *e, double *aorl, double *dm,
		 double *u, int *jstat);
#endif

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
  double tmp[3][3], rot[3][3], r[3], v[3];
  double ma, ea, se, ce, st, ct, f, df, delta;
  double sef, fac, raq;
  int i;

#ifdef TEST
  double u[13];
  int jform, jstat;

  jform = eltype;

  sla_el2ue_(&epoch, &jform,
	     &epoch, &incl, &anode, &longperi,
	     &aq, &ecc, &lm, &nn,
	     u, &jstat);

  fprintf(stderr,
	  "SLALIB universal elements: r_0=(%le, %le, %le) v_0=(%le, %le, %le) alpha=%le r_0=%le\n",
	  u[3], u[4], u[5],
	  u[6]/58.1324409, u[7]/58.1324409, u[8]/58.1324409,  /* canonical to normal */
	  -u[1]/(58.1324409*58.1324409),
	  u[9]);
#endif

  /* Set source type */
  src->type = SOURCE_ELEM;

  /* Generate rotation matrix */
  m_identity(tmp);
  euler_rotate(tmp, 3, -longperi);
  euler_rotate(tmp, 1, -incl);
  euler_rotate(tmp, 3, -anode);
  m_x_m(gcrs2ecl, tmp, rot);

  /* Reference epoch */
  src->ref_tdb = epoch;

  if(eltype == SOURCE_ELEM_COMET) {
    /* Everything is relative to perihelion for comets.
       All of the anomalies are zero here, so conversion to
       universal elements is straightforward. */

    src->mu = GMSUN * DAY*DAY / (AU*AU*AU);  /* AU^3 / day^2 */

    /* Position vector */
    r[0] = aq;
    r[1] = 0;
    r[2] = 0;

    /* Velocity vector */
    raq = 1.0 / aq;

    v[0] = 0;
    v[1] = sqrt(src->mu * (1.0 + ecc)*raq);
    v[2] = 0;

    src->alpha = src->mu * (1.0 - ecc)*raq;
    src->sqrtalpha = sqrt(fabs(src->alpha));
    src->nn = src->sqrtalpha * (1.0 - ecc)*raq;
    src->rref = aq;
  }
  else {
    if(eltype == SOURCE_ELEM_MAJOR) {
      /* Mean anomaly at reference epoch */
      ma = lm - anode - longperi;
      
      /* Compute alpha and mu from mean motion */
      src->alpha = aq*aq * nn*nn;
      src->sqrtalpha = aq * nn;
      src->mu = aq * src->alpha;
    }
    else {  /* SOURCE_ELEM_MINOR */
      /* Mean anomaly at reference epoch */
      ma = lm;
      
      /* Compute mu and alpha */
      raq = 1.0 / aq;

      src->mu = GMSUN * DAY*DAY / (AU*AU*AU);  /* AU^3 / day^2 */
      src->alpha = src->mu * raq;
      src->sqrtalpha = sqrt(src->alpha);

      /* Compute mean motion */
      nn = src->sqrtalpha * raq;
    }

    /* Reduce mean anomaly */
    ma = fmod(ma, TWOPI);

    /* Eccentric anomaly at reference epoch: solve Kepler's equation */
    if(ecc > 0.8)
      ea = (M_PI * ecc + ma) / (1.0 + ecc);  /* for stability */
    else
      ea = ma;
    
    for(i = 0; i < KEPLER_MAXITER; i++) {
      inline_bare_sincos(ea, se, ce);

      f = ea - ecc * se - ma;
      df = 1.0 - ecc * ce;
      delta = f / df;
      
      if(fabs(delta) < KEPLER_PREC) {
	/* I think that's enough... */
	break;
      }

      ea -= delta;
    }

    /* Compute (r/a) times sin and cos of true anomaly at reference epoch */
    sef = sqrt(1.0 - ecc*ecc);

    st = se * sef;
    ct = ce - ecc;

    /* r = a (1 - e cos E) */
    src->rref = aq * df;

    /* Position vector */
    r[0] = aq * ct;
    r[1] = aq * st;
    r[2] = 0;
    
    /* Velocity vector */
    fac = nn * aq / df;
    
    v[0] = -fac * se;
    v[1] = fac * sef * ce;
    v[2] = 0;

    src->nn = nn;
  }

  src->muorref = src->mu / src->rref;
  src->rvref = v_d_v(r, v);

  m_x_v(rot, r, src->ref_n);
  m_x_v(rot, v, src->ref_dndt);

  src->psi = 0;
}

void source_place (struct observer *obs,
		   struct source *src,
		   struct jpleph_table *jpltab,
		   unsigned char mask,    /* parts we want */
		   double *s,             /* result */
		   double *dsdt,          /* src motion only */
		   double *pr_r) {
  double vb[3];
  double pr, nf, dt;
  double c[4];
  double psi, k, dk, delta, f, g, rr;
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
	    vb[i] = obs->bsv[i];
	  }
	}
	else {
	  for(i = 0; i < 3; i++) {
	    s[i]  = obs->bsp[i];
	    vb[i] = obs->bsv[i];
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
	(void) jpleph_fetch(jpltab, obs->tdb, src->type, s, vb);

	if(mask & TR_PLX) {
	  for(i = 0; i < 3; i++)
	    s[i] -= obs->bop[i];
	}
      }
    }
    else {  /* Orbital elements */

      dt = obs->tdb - src->ref_tdb;

      /* Iterative solution to universal Kepler's equation */
      psi = src->psi;  /* initial guess is previous answer */

      for(i = 0; i < KEPLER_MAXITER; i++) {
	/* Evaluate Stumpff functions: c receives s^k c_k(x) */
	stumpff(psi, src->alpha, src->sqrtalpha, c);
	
	/* Newton-Raphson iteration.  Derivative dk/ds uses
	   d/ds (s^k c_k(x)) = c_k-1(x) */
	k  = src->rref * c[1] + src->rvref * c[2] + src->mu * c[3] - dt;
	dk = src->rref * c[0] + src->rvref * c[1] + src->mu * c[2];
	
	delta = k / dk;

#ifdef TEST
	fprintf(stderr, "iter %d: %lf %lf %lf %lf (%lf %lf %lf %lf)\n",
		i+1, psi, k, dk, delta, c[0], c[1], c[2], c[3]);
#endif	

	if(fabs(delta) < KEPLER_PREC) {
	  /* I think that's enough... */
	  break;
	}
	
	psi -= delta;
      }
      
      /* Compute result */
      f = 1.0 - c[2] * src->muorref;
      g = dt - c[3] * src->mu;
      
      for(i = 0; i < 3; i++)
	s[i] = f * src->ref_n[i] + g * src->ref_dndt[i];
      
      rr = 1.0 / sqrt(v_d_v(s, s));
      
      f = -c[1] * src->muorref * rr;
      g = 1.0 - c[2] * src->mu * rr;
      
      for(i = 0; i < 3; i++)
	vb[i] = f * src->ref_n[i] + g * src->ref_dndt[i];

#ifdef TEST
      fprintf(stderr,
	      "Heliocentric state vector: (%le, %le, %le)\n",
	      s[0], s[1], s[2]);
#endif

      /* Convert to barycentric or topocentric as requested */
      if(mask & TR_PLX) {
	for(i = 0; i < 3; i++) {
	  s[i]  -= obs->hop[i];
	  vb[i] += obs->bsv[i];
	}
      }
      else {
	for(i = 0; i < 3; i++) {
	  s[i]  += obs->bsp[i];
	  vb[i] += obs->bsv[i];
	}
      }
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
    observer_ast2obs(obs, s, dsdt, pr, mask);

  if(pr_r)
    *pr_r = pr;
}

