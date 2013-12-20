#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lfa.h"

/* Aim and boresight vectors to roll and pitch angles by decomposing
   posture matrix.  Follows Wallace (2002) equations.  Returns false
   if there were no solutions, otherwise true.  Intended for
   generating demands to send to a telescope drive system. */

int mount_ab2rp (double *aim, double *daimdt,
		 double *bore,
		 double snp, double cnp,
		 unsigned char flip,
		 double pos[3][3], double *r, double *p,
		 double dposdt[3][3], double *drdt, double *dpdt) {
  double salph, calphsq, calph, sp, cp, x, y, sr, cr, npsq, nrsq, np, nr;
  double dposdt_a[3][3], dposdt_b[3][3];
  double dsalphdt, dcalphdt, dspdt, dcpdt, dxdt, dydt, dsrdt, dcrdt;
  int i, j;

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

  /* Normalization constants */
  npsq = 1.0/(sp*sp+cp*cp);
  nrsq = 1.0/(sr*sr+cr*cr);

  np = sqrt(npsq);
  nr = sqrt(nrsq);

  /* Resulting posture matrix */
  m_identity(pos);
  euler_rotate_sc(pos, 2, sp*np, cp*np);
  euler_rotate_sc(pos, 1, snp, cnp);
  m_copy(pos, dposdt_b);  /* for later */
  euler_rotate_sc(pos, 3, -sr*nr, cr*nr);

  /* Resulting roll and pitch angles */
  *r = atan2(sr, cr);
  *p = atan2(sp, cp);

  /* Time derivatives, if requested */
  if(drdt || dpdt) {
    dsalphdt = daimdt[2];

    if(calph > 0)
      dcalphdt = (aim[0]*daimdt[0] + aim[1]*daimdt[1] - 
		  bore[1]*daimdt[2]*snp) / calph;
    else
      dcalphdt = 0;  /* pole */

    dspdt = dsalphdt*bore[0] - dcalphdt*bore[2];
    dcpdt = dsalphdt*bore[2] + dcalphdt*bore[0];

    dxdt = dcpdt*bore[0] + dspdt*bore[2];
    dydt = (dspdt*bore[0] - dcpdt*bore[2])*snp;

    dsrdt = dxdt*aim[1] + x*daimdt[1] - dydt*aim[0] - y*daimdt[0];
    dcrdt = dxdt*aim[0] + x*daimdt[0] + dydt*aim[1] + y*daimdt[1];

    /* Product rule for resulting posture matrix */
    m_identity(dposdt_a);
    euler_drot_sc(dposdt_a, 2,
		   (dspdt*cp - dcpdt*sp)*cp*np*npsq,
		   (dcpdt*sp - dspdt*cp)*sp*np*npsq);
    euler_rotate_sc(dposdt_a, 1, snp, cnp);
    euler_rotate_sc(dposdt_a, 3, -sr*nr, cr*nr);

    euler_drot_sc(dposdt_b, 3,
		  -(dsrdt*cr - dcrdt*sr)*cr*nr*nrsq,
		   (dcrdt*sr - dsrdt*cr)*sr*nr*nrsq);

    for(j = 0; j < 3; j++)
      for(i = 0; i < 3; i++)
	dposdt[j][i] = dposdt_a[j][i] + dposdt_b[j][i];

    /* Resulting roll and pitch angles */
    if(drdt)
      *drdt = (dsrdt*cr - dcrdt*sr) * nrsq;
    if(dpdt)
      *dpdt = (dspdt*cp - dcpdt*sp) * npsq;
  }

  return(1);
}

void mount_pa (double *aimp, double *daimpdt,
	       double *bore,
	       double pos[3][3], double dposdt[3][3],
	       double *a, double *dadt) {
  double borep[3], dborepdt_a[3], dborepdt_b[3], dborepdt[3];
  double sa, ca, dsadt, dcadt;
  int i;

  /* Transform PA unit vector */
  mt_x_v(pos, aimp, borep);

  /* Angle, from Wallace Sect. 4.3 */
  sa = borep[1]*bore[0] - borep[0]*bore[1];
  ca = borep[2]*(bore[0]*bore[0] + bore[1]*bore[1]) -
       bore[2]*(borep[0]*bore[0] + borep[1]*bore[1]);
    
  *a = atan2(sa, ca);

  if(dadt) {
    /* Product rule for time derivative of transformed PA unit vector */
    mt_x_v(pos, daimpdt, dborepdt_a);
    mt_x_v(dposdt, aimp, dborepdt_b);

    for(i = 0; i < 3; i++)
      dborepdt[i] = dborepdt_a[i] + dborepdt_b[i];

    dsadt = dborepdt[1]*bore[0] - dborepdt[0]*bore[1];
    dcadt = dborepdt[2]*(bore[0]*bore[0] + bore[1]*bore[1]) -
            bore[2]*(dborepdt[0]*bore[0] + dborepdt[1]*bore[1]);
    
    *dadt = (dsadt*ca - dcadt*sa) / (sa*sa + ca*ca);
  }
}

/* Roll and pitch angles to posture matrix, defined so A = P B.
   Intended for processing encoder readings. */

void mount_rp2pos (double r, double p,
		   double snp, double cnp,
		   double pos[3][3]) {

  m_identity(pos);
  euler_rotate(pos, 2, p);
  euler_rotate_sc(pos, 1, snp, cnp);
  euler_rotate(pos, 3, -r);
}

