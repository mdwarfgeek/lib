#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lfa.h"
#include "cvtunit.h"
#include "util.h"

void euler_rotate (double m[3][3], int axis, double angle) {
  double s, c;

  inline_sincos(angle, s, c);
  euler_rotate_sc(m, axis, s, c);
}

void euler_rotate_sc (double m[3][3], int axis, double sa, double ca) {
  int i = 0, j = 1, k;
  double atmp, btmp;

  switch(axis) {
  case 1:
    i = 1;
    j = 2;
    break;
  case 2:
    i = 2;
    j = 0;
    break;
  case 3:
    i = 0;
    j = 1;
    break;
  }

  for(k = 0; k < 3; k++) {
    atmp =  ca * m[i][k] + sa * m[j][k];
    btmp = -sa * m[i][k] + ca * m[j][k];
    m[i][k] = atmp;
    m[j][k] = btmp;
  }
}

void euler_drot_sc (double m[3][3], int axis, double dsa, double dca) {
  int h = 2, i = 0, j = 1, k;
  double atmp, btmp;

  switch(axis) {
  case 1:
    h = 0;
    i = 1;
    j = 2;
    break;
  case 2:
    h = 1;
    i = 2;
    j = 0;
    break;
  case 3:
    h = 2;
    i = 0;
    j = 1;
    break;
  }

  for(k = 0; k < 3; k++) {
    atmp =  dca * m[i][k] + dsa * m[j][k];
    btmp = -dsa * m[i][k] + dca * m[j][k];
    m[i][k] = atmp;
    m[j][k] = btmp;
    m[h][k] = 0;
  }
}

void m_x_v (double m[3][3], double vi[3], double vo[3]) {
  int i, j;

  for(j = 0; j < 3; j++) {
    vo[j] = 0;

    for(i = 0; i < 3; i++)
      vo[j] += m[j][i] * vi[i];
  }
}

void mt_x_v (double m[3][3], double vi[3], double vo[3]) {
  int i, j;

  for(j = 0; j < 3; j++) {
    vo[j] = 0;

    for(i = 0; i < 3; i++)
      vo[j] += m[i][j] * vi[i];
  }
}

void m_x_m (double a[3][3], double b[3][3], double c[3][3]) {
  int i, j, k;

  for(j = 0; j < 3; j++)
    for(i = 0; i < 3; i++) {
      c[j][i] = 0;

      for(k = 0; k < 3; k++)
	c[j][i] += a[j][k] * b[k][i];
    }
}

double v_renorm (double v[3]) {
  double nf;

  nf = v_d_v(v, v);
  if(nf > 0) {
    nf = 1.0 / sqrt(nf);

    v[0] *= nf;
    v[1] *= nf;
    v[2] *= nf;
  }
  else
    nf = 1.0;

  return nf;
}

void ad_to_v (double a, double d, double v[3]) {
  double sa, ca, sd, cd;

  /* Precompute these */
  inline_sincos(a, sa, ca);
  inline_sincos(d, sd, cd);

  /* Unit vector of direction cosines */
  v[0] = ca*cd;
  v[1] = sa*cd;
  v[2] = sd;
}

void v_to_ad (double v[3], unsigned char flip, double *a, double *d) {
  if(flip) {
    *a = atan2(-v[1], -v[0]);
    *d = atan2(v[2], -sqrt(v[0]*v[0]+v[1]*v[1]));
  }
  else {
    *a = atan2(v[1], v[0]);
    *d = atan2(v[2], sqrt(v[0]*v[0]+v[1]*v[1]));
  }
}

void v_to_ad_dt (double v[3], double dvdt[3], unsigned char flip,
		 double *dadt, double *dddt) {
  double rsq;
  
  rsq = v[0]*v[0] + v[1]*v[1];
  if(rsq > 0) {
    *dadt = (v[0]*dvdt[1] - v[1]*dvdt[0]) / rsq;
    *dddt = (dvdt[2]*rsq - v[2]*(v[0]*dvdt[0] + v[1]*dvdt[1]))
      / ((rsq + v[2]*v[2])*sqrt(rsq));

    if(flip)
      *dddt *= -1;
  }
  else {
    /* Pole */
    *dadt = 0;
    *dddt = 0;
  }
}

void v_to_az (double v[3], unsigned char flip, double *a, double *z) {
  if(flip) {
    *a = atan2(-v[1], -v[0]);
    *z = atan2(-sqrt(v[0]*v[0]+v[1]*v[1]), v[2]);
  }
  else {
    *a = atan2(v[1], v[0]);
    *z = atan2(sqrt(v[0]*v[0]+v[1]*v[1]), v[2]);
  }
}

double v_angle_v (double u[3], double v[3]) {
  double t[3];

  /* Robust method, using ratio of vector and scalar products:
     |u x v| = |u||v| sin theta
      u . v  = |u||v| cos theta */
  v_x_v(u, v, t);

  return atan2(sqrt(v_d_v(t, t)),
	            v_d_v(u, v));
}

