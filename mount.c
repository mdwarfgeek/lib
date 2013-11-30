#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lfa.h"

/* Aim and boresight vectors to roll and pitch angles by decomposing
   posture matrix.  Follows Wallace (2002) equations.  Returns false
   if there were no solutions, otherwise true. */

int mount_ab2rp (double *aim, double *daimdt,
		 double *bore,
		 double snp, double cnp,
		 unsigned char flip,
		 double *r, double *p, double *drdt, double *dpdt) {
  double salph, calphsq, calph, sp, cp, x, y, sr, cr;
  double dsalphdt, dcalphdt, dspdt, dcpdt, dxdt, dydt, dsrdt, dcrdt;

  /* From Wallace, p = atan(salph/calph) - atan(bore_z / bore_x)
                     = atan((salph*bore_x - calph*bore_z) /
		            (salph*bore_z + calph*bore_x))
                     = atan(sp/cp)
  */

  /* sin and cos of first part */
  salph = aim[2] + snp*bore[1];
  calphsq = aim[0]*aim[0] + aim[1]*aim[1] -
    bore[1]*(2*aim[2]*snp + bore[1]) - snp*snp;

  if(calphsq < 0)  /* no solutions */
    return(0);

  calph = sqrt(calphsq);

  if(flip)
    calph = -calph;

  /* sin and cos of p */
  sp = salph*bore[0] - calph*bore[2];
  cp = salph*bore[2] + calph*bore[0];

  /* From Wallace, r = atan((aim_y*x - aim_x*y) / (aim_x*x + aim_y*y))
     x = bore_x*cp + bore_z*sp
     y = bore_y*cnp + (bore_x*sp - bore_z*cp) * snp
  */

  x = cp*bore[0] + sp*bore[2];
  y = cnp*bore[1] + (sp*bore[0] - cp*bore[2])*snp;

  /* sin and cos of r */
  sr = x*aim[1] - y*aim[0];
  cr = x*aim[0] + y*aim[1];

  /* Resulting roll and pitch angles */
  *r = atan2(sr, cr);
  *p = atan2(sp, cp);

  /* Time derivatives, if requested */
  if(drdt || dpdt) {
    dsalphdt = daimdt[2];

    if(calph > 0)
      dcalphdt = (aim[0]*daimdt[0] + aim[1]*daimdt[1] - 
		  2*bore[1]*daimdt[2]*snp) / calph;
    else
      dcalphdt = 0;  /* pole */

    dspdt = dsalphdt*bore[0] - dcalphdt*bore[2];
    dcpdt = dsalphdt*bore[2] + dcalphdt*bore[0];

    dxdt = dcpdt*bore[0] + dspdt*bore[2];
    dydt = (dspdt*bore[0] - dcpdt*bore[2])*snp;

    dsrdt = dxdt*aim[1] + x*daimdt[1] - dydt*aim[0] - y*daimdt[0];
    dcrdt = dxdt*aim[0] + x*daimdt[0] + dydt*aim[1] + y*daimdt[1];

    if(drdt)
      *drdt = (dsrdt*cr - dcrdt*sr) / (sr*sr + cr*cr);
    if(dpdt)
      *dpdt = (dspdt*cp - dcpdt*sp) / (sp*sp + cp*cp);
  }

  return(1);
}

/* Roll and pitch angles to posture matrix, defined so A = P B */

void mount_rp2pos (double r, double p,
		   double snp, double cnp,
		   double pos[3][3]) {

  m_identity(pos);
  euler_rotate(pos, 2, p);
  euler_rotate_sc(pos, 1, snp, cnp);
  euler_rotate(pos, 3, -r);
}

