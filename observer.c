#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "astro.h"
#include "cvtunit.h"
#include "util.h"

void observer_init (struct observer *obs,
		    struct jpleph_table *jpltab,
		    struct jpleph_table *tetab,
		    struct iers_table *iertab,
		    double longitude,  /* WGS84 longitude, E positive */
		    double latitude,   /* WGS84 latitude, N positive */
		    double height) {   /* Height above geoid, m */
  double sinphi, cosphi, sinlam, coslam;
  double u, z;

  double ang[NPNANG];

  obs->jpltab = jpltab;
  obs->tetab = tetab;
  obs->iertab = iertab;

  obs->latitude = latitude;
  obs->longitude = longitude;
  obs->height = height;

  /* Geocentric location of observer */
  dsincos(latitude, &sinphi, &cosphi);
  dsincos(longitude, &sinlam, &coslam);

  geoc(sinphi, cosphi, height, &u, &z);

  obs->itgop[0] = u*coslam / AU;
  obs->itgop[1] = u*sinlam / AU;
  obs->itgop[2] = z / AU;

  /* Earth gravitational potential at observer / c^2 */
  obs->gpearth = GMEARTH / (LIGHT*LIGHT * sqrt(u*u + z*z));

  /* Observer longitude and latitude rotation matrices */
  m_identity(obs->lmm);
  euler_rotate_sc(obs->lmm, 3, sinlam, coslam);

  m_identity(obs->phm);
  euler_rotate_sc(obs->phm, 2, cosphi, sinphi);

  /* Minimum impact parameter for light deflection: solar radius in AU */
  obs->minbdefl = RSUN / AU;

  /* GCRS to mean ecliptic of J2000, for lazy people */
  pfb06ang(0.0, ang);

  m_identity(obs->eclm);
  euler_rotate(obs->eclm, 3,  ang[PNANG_GAM]);
  euler_rotate(obs->eclm, 1,  ang[PNANG_PHI]);
  euler_rotate(obs->eclm, 3, -ang[PNANG_PSI]);

  /* These are optional, so set default values that do nothing */
  obs->dut1 = 0;
  obs->xp = 0;
  obs->yp = 0;

  obs->dtdb = 0;
}

int observer_update (struct observer *obs,
		     double tt,
		     unsigned char mask) {
  int rv, i;
  double jdtk, jctk, sp;

  double tmp[3][3];

  double tei, dteidt, tev[3], dtevdt[3];

  /* Time arguments: TT Julian days and centuries since 2000.0 */
  jdtk = (tt-J2K);
  jctk = jdtk / (100*JYR);

  if(mask & OBSERVER_UPDATE_IERS) {
    if(obs->iertab) {
      /* Update Earth orientation parameters from IERS tabulation.
	 Strictly, the argument should be UTC, but TT is used here. 
	 dX and dY are ignored, because we use IAU 2000B rather
         than 2000A nutation, so strictly they don't apply. */
      rv = iers_fetch(obs->iertab, tt,
		      &(obs->dut1), &(obs->xp), &(obs->yp),
		      (double *) NULL, (double *) NULL);
      if(rv < 0)
	return(rv);
      
      /* Convert units */
      obs->xp *= AS_TO_RAD;
      obs->yp *= AS_TO_RAD;
    }
    
    /* TIO locator, kind of pointless if we are doing nutation,
       which is only accurate to about 1 mas with 2000B, but
       it is included anyway. */
    sp = -47e-6 * jctk;
    
    /* Polar motion = R_1(-yp) * R_2(-xp) * R_3(s') */
    m_identity(obs->pmm);
    euler_rotate(obs->pmm, 3, sp*AS_TO_RAD);
    euler_rotate(obs->pmm, 2, -obs->xp);
    euler_rotate(obs->pmm, 1, -obs->yp);

    /* Longitude and polar motion matrix */
    m_x_m(obs->lmm, obs->pmm, obs->lpm);

    /* Geocentre-Observer position vector in TIRS */
    mt_x_v(obs->pmm, obs->itgop, obs->tigop);

    /* Geocentre-Observer velocity vector in TIRS */
    obs->tigov[0] = -obs->tigop[1] * TWOPI * (1.0 + ERADAY);
    obs->tigov[1] =  obs->tigop[0] * TWOPI * (1.0 + ERADAY);
    obs->tigov[2] = 0;
  }

  if(mask & OBSERVER_UPDATE_NUT)
    /* Nutation and CIO locator */
    nut00b(jctk, &(obs->dpsi), &(obs->deps), &(obs->sxy), (double *) NULL);

  if(mask & OBSERVER_UPDATE_PFB) {
     /* Precession */
     pfb06ang(jctk, obs->ang_pfb);

     /* Add in (stored) nutation */
     obs->ang_pnfb[PNANG_GAM]  = obs->ang_pfb[PNANG_GAM];
     obs->ang_pnfb[PNANG_PHI]  = obs->ang_pfb[PNANG_PHI];
     obs->ang_pnfb[PNANG_PSI]  = obs->ang_pfb[PNANG_PSI]  + obs->dpsi;
     obs->ang_pnfb[PNANG_EPSA] = obs->ang_pfb[PNANG_EPSA] + obs->deps;

     /* Compute matrix */
     makecim(obs->ang_pnfb, obs->sxy, 0.0, 0.0, obs->cim);
  }

  if(mask & OBSERVER_UPDATE_ERA) {
    /* Earth rotation angle */
    obs->era = TWOPI*fmod(ERA2K +
			  fmod(jdtk, 1.0) +
			  fmod(ERADAY*jdtk, 1.0) -
			  obs->dut1*(1.0+ERADAY)/DAY,
			  1.0);
    
    /* Earth rotation matrix = R_3(era) */
    m_identity(obs->erm);
    euler_rotate(obs->erm, 3, obs->era);

    /* GCRS - TIRS */
    m_x_m(obs->erm, obs->cim, tmp);

    /* Geocentre-Observer position and velocity in GCRS */
    mt_x_v(tmp, obs->tigop, obs->gop);
    mt_x_v(tmp, obs->tigov, obs->gov);

    /* GCRS - Topocentric (-h, delta) */
    m_x_m(obs->lpm, tmp, obs->ctm);
  }

  if(mask & OBSERVER_UPDATE_TDB && obs->tetab) {
    /* Fetch time ephemeris integral and Earth velocity vector */
    rv = jpleph_fetch(obs->tetab, tt, TIMEEPH_TEI, &tei, &dteidt);
    if(rv < 0)
      return(rv);

    rv = jpleph_fetch(obs->tetab, tt, TIMEEPH_TEV, tev, dtevdt);
    if(rv < 0)
      return(rv);

    /* Compute delta(TDB) in seconds */
    obs->dtdb = ZTDB
      + v_d_v(tev, obs->gop)*AU*AU / (LIGHT*LIGHT*DAY)  /* topocentric */
      + tei*DAY / (1.0 - obs->tetab->lc);               /* geocentric */
  }

  obs->tdb = tt + obs->dtdb / DAY;

  if(mask & OBSERVER_UPDATE_SOLSYS) {
    /* Components of Earth position and velocity relative to SSB and
       heliocentre from ephemerides.  TT is used as an approximation
       to TDB if dtdb is not available. */
    rv = jpleph_fetch(obs->jpltab, obs->tdb, JPLEPH_EMB, obs->bep, obs->bev);
    if(rv < 0)
      return(rv);
    
    rv = jpleph_fetch(obs->jpltab, obs->tdb, JPLEPH_MOON, obs->emp, obs->emv);
    if(rv < 0)
      return(rv);
    
    rv = jpleph_fetch(obs->jpltab, obs->tdb, JPLEPH_SUN, obs->bsp, obs->bsv);
    if(rv < 0)
      return(rv);
    
    for(i = 0; i < 3; i++) {
      /* SSB->Earth = SSB->EMB - Earth->EMB */
      obs->bep[i] -= obs->jpltab->emfac * obs->emp[i];
      obs->bev[i] -= obs->jpltab->emfac * obs->emv[i];
      
      /* Sun->Earth = SSB->Earth - SSB->Sun */
      obs->hep[i] = obs->bep[i] - obs->bsp[i];
      obs->hev[i] = obs->bev[i] - obs->bsv[i];
    }
  }

  if(mask & OBSERVER_UPDATE_ERA || mask & OBSERVER_UPDATE_SOLSYS) {
    for(i = 0; i < 3; i++) {
      /* Barycentric and heliocentric position vectors of observer */
      obs->bop[i] = obs->bep[i] + obs->gop[i];
      obs->hop[i] = obs->hep[i] + obs->gop[i];

      /* Aberration vector (also used for Doppler shift) */
      obs->vab[i] = (obs->bev[i] + obs->gov[i]) * AU / (DAY*LIGHT);
    }

    /* Lorentz factor */
    obs->gab = 1.0 / sqrt(1.0 - v_d_v(obs->vab, obs->vab));

    /* Distance from observer to Sun, AU */
    obs->hdist = sqrt(v_d_v(obs->hop, obs->hop));

    /* Gravitational potential at observer due to Sun / c^2 */
    obs->gpsun = GMSUN / (LIGHT*LIGHT * obs->hdist * AU);
  }

  obs->tt = tt;

  return(0);
}

