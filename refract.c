#include <stdio.h>
#include <math.h>

#include "lfa.h"
#include "cvtunit.h"
#include "util.h"

#define TROPO_HEIGHT  11000  /* m */
#define TLR          0.0065  /* K/m */

/* Computes refraction constants following Saastamoinen (1972).
 * This treatment is better than 0.2 arcsec for ZD < 80 degrees
 * (below 0.1 arcsec by 75 degrees), so should be adequate for most
 * telescope pointing purposes, and it's a lot easier to work with
 * than numerically integrating an atmospheric model.  The formula
 * blows up for ZD > 85 degrees or so, and should not be used there.
 */

void refract_const (double temperat,   /* K */
		    double humidity,   /* 0-1 */
		    double pressure,   /* mbar */
		    double wavelength, /* um */
		    double height,     /* m */
		    double *refco) {
  double ttropo, ptropo, ph2o, scl, q, pt;

  /* Compute temperature and pressure at tropopause */
  ttropo = temperat + TLR * (height-TROPO_HEIGHT);
  ptropo = pressure * pow(ttropo / temperat, 5.26);

  /* Compute partial pressure of water vapour */
  ph2o = humidity * pow(temperat / 247.1, 18.36);

  /* Compute scale factor - arcsec to radians, and
   * correction to refractivity for wavelength dependence, Eq. 38 */
  scl = AS_TO_RAD * 170.2649 / (173.3 - 1.0/(wavelength*wavelength));

  /* Compute refraction constants, Eq. 31a and 30a */
  q = (pressure - 0.156*ph2o) / temperat;
  pt = 0.000288 * 1.0e-6 * (pressure*temperat + 0.190 * ptropo*ttropo);

  refco[0] = scl*(16.271 * q - 0.0000749 * pressure);
  refco[1] = scl*(16.271 * 0.0000394 * q*q - 0.0000749 * pressure + 5*pt);
  refco[2] = scl*(3*pt - 0.013 * 1.0e-6 * pressure*pressure / temperat);
  refco[3] = scl*(-0.014 * 1.0e-12 * (pressure*temperat*temperat +
				      0.64*ptropo*ttropo*ttropo));
  refco[4] = scl*(0.0003 * 1.0e-15 * (pressure*temperat*temperat*temperat +
				      2*ptropo*ttropo*ttropo*ttropo));
}

void refract_corr (double *refco, double tanz, double *refr, double *deriv) {
  double incr;
  int i;

  incr = tanz*tanz;

  *refr  = refco[NREFCO-1];
  *deriv = (2*NREFCO-1)*refco[NREFCO-1];

  for(i = NREFCO-2; i >= 0; i--) {
    *refr  = refco[i] + incr * (*refr);
    *deriv = (2*i+1) * refco[i] + incr * (*deriv);
  }

  *refr  *= tanz;
  *deriv *= (1.0+incr);  /* sec^2 z */
}

/* This treatment is good to about 1 arcsec at 85 degrees, i.e. better than
   the refraction model itself.  It avoids use of transcendental functions
   by Taylor expanding. */

void refract_vec (double *refco, unsigned char unref,
		  double *vi, double *vo, double *dvidt, double *dvodt) {
  double z, dzdt;
  double zsq, rezsq, rsq, tansqzd, secsqzd;
  int i;
  double fac, zfac, refsum, dersum, wt, wtfac;

  double drefsum, sdersum, dtansqzddt, dfacdt, dwtdt;

  /* Hold constant below about 3 degrees */
  z = vi[2];
  dzdt = dvidt ? dvidt[2] : 0;
  if(z < 0.05) {
    z = 0.05;
    dzdt = 0;
  }

  /* Compute tan^2 and sec^2 input zenith distance */
  zsq = z*z;
  rsq = vi[0]*vi[0] + vi[1]*vi[1];
  rezsq = 1.0 / zsq;
  tansqzd = rsq * rezsq;
  secsqzd = (rsq + zsq) * rezsq;

  /* Use it to compute the polynomial (in tan^2 zd_in) parts of the
     refraction formula. */
  refsum = refco[NREFCO-1];
  dersum = (2*NREFCO-1)*refco[NREFCO-1];

  for(i = NREFCO-2; i >= 0; i--) {
    refsum = refco[i] + tansqzd * refsum;
    dersum = (2*i+1) * refco[i] + tansqzd * dersum;
  }

  if(unref) {
    /* Unrefracting is easy, formula is in observed ZD. */
    wtfac = 1.0;
    wt = refsum;
  }
  else {
    /* Compute -(f/f')/tan(zd_in) ~= -delta/tan(zd_in)
       where delta=(zd_out-zd_in), for one Newton-Raphson iteration.
       Velocity is approximated using just leading order terms in refsum. */
    wtfac = 1.0 / (1.0 + dersum * secsqzd);
    wt = -refsum * wtfac;
  }

  /* The vector of direction cosines we want out is:
     cos(az)*sin(zd_out)
     sin(az)*sin(zd_out)
             cos(zd_out)
     Taylor expanding the sin and cos about zd_in to second order yields:
     sin(zd_out) ~= sin(zd_in) + cos(zd_in)*(zd_out-zd_in) - 0.5*sin(zd_in)*(zd_out-zd_in)**2
                  = sin(zd_in) (1 + delta * cot(zd_in) - 0.5 * delta**2)
                  = sin(zd_in) (1 + wt - 0.5 * wt**2 * tan(zd_in)**2)
     cos(zd_out) ~= cos(zd_in) - sin(zd_in)*(zd_out-zd_in) - 0.5*cos(zd_in)*(zd_out-zd_in)**2
                  = cos(zd_in) (1 - delta * tan(zd_in) - 0.5 * delta**2)
                  = cos(zd_in) (1 - (wt + 0.5 * wt**2) * tan(zd_in)**2)
  */

  /* Compute scale factor for x,y */
  fac = 1.0 + wt - 0.5*wt*wt*tansqzd;
  
  /* Compute scale factor for z */
  zfac = wt * (1.0 + 0.5*wt);

  /* Compute refracted vector */
  vo[0] = vi[0] * fac;
  vo[1] = vi[1] * fac;
  vo[2] = vi[2] - z * zfac * tansqzd;

  /* Compute velocity vector if requested */
  if(dvidt && dvodt) {
    /* First derivative of "refsum" with respect to tan^2 zd 
       and derivative of "dersum" with respect to tan^2 zd */
    drefsum = (NREFCO-1)*refco[NREFCO-1];
    sdersum = (2*NREFCO-1)*(NREFCO-1)*refco[NREFCO-1];

    for(i = NREFCO-2; i > 0; i--) {
      drefsum = i * refco[i] + tansqzd * drefsum;
      sdersum = (2*i+1)*i * refco[i] + tansqzd * sdersum;
    }

    /* Derivatives of expressions above */
    dtansqzddt = 2.0*((vi[0]*dvidt[0] + vi[1]*dvidt[1]) * rezsq - tansqzd * dzdt / z);

    if(unref)
      dwtdt = drefsum * dtansqzddt;
    else {
      dwtdt = -wtfac * (drefsum +
			wt * (dersum + sdersum*secsqzd)) * dtansqzddt;
    }

    dfacdt = dwtdt - wt*(dwtdt*tansqzd + 0.5*wt*dtansqzddt);

    dvodt[0] = dvidt[0] * fac + vi[0] * dfacdt;
    dvodt[1] = dvidt[1] * fac + vi[1] * dfacdt;
    dvodt[2] = dvidt[2] - dzdt *               zfac * tansqzd
                        -    z * dwtdt * (1.0 + wt) * tansqzd
                        -    z *               zfac * dtansqzddt;
  }
}
