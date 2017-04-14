#include <math.h>

#include "lfa.h"

/* Area of circular segment A(r^2, r^2 sin theta, r^2 cos theta)
   Compute as area of circular sector 0.5 * r^2 * theta minus
   area of triangle 0.5 * r^2 * sin theta.  The angle can be
   larger than 90 degrees so we use atan2. */

static inline double circseg (double rsq, double rsqst, double rsqct) {
  return(0.5*(rsq * atan2(rsqst, rsqct) - rsqst));
}

/* Fraction of pixel inside circle of radius r, general case.
   x, y are coordinates of the centre of the pixel relative to
   the centre of the circle.  */

double pixovcirc (double x, double y, double r) {
  double rsq;
  double t;
  double x1, x2, y1, y2;
  double x1sq, x2sq, y1sq, y2sq;
  double hx1, hx2, hy1, hy2;
  double a, b, st, ct, area;

  rsq = r*r;

  /* Use symmetry to reduce to first octant with x > 0, y > 0 and x > y */
  x = fabs(x);
  y = fabs(y);
  r = fabs(r);  /* just in case */
  if(y > x) {
    t = x;
    x = y;
    y = t;
  }

  /* Evaluate the most straightforward case first: if the vertex
     furthest from the origin is inside the circle, then the whole
     pixel must be inside the circle. */
  x2 = x+0.5;
  y2 = y+0.5;

  x2sq = x2*x2;
  y2sq = y2*y2;

  if(x2sq+y2sq <= rsq)  /* entirely inside */
    return(1.0);
  else {
    /* We'll need the other ones too */
    x1 = x-0.5;
    y1 = y-0.5;

    x1sq = x1*x1;
    y1sq = y1*y1;
    
    /* The distinction here is whether intersection points are in the
       same side of the square or different sides.  A side can have
       0, 1 or 2 intersections.  If 1, a vertex must be included, if
       2, no vertex is included. */

    /* Handle the 2-intersection cases first because they allow us to
       precompute quantities we'll need later.  These are circular
       segments that are outside the pixel (so their contribution to
       the area is negative).  In each case we can just use the double
       angle formula because both angles are the same. */
    area = 0;

    if(x1sq < rsq) {
      hx1 = sqrt(rsq - x1sq);
      if(-hx1 > y1 && hx1 < y2) {
        /* r cos t1 = x1 */
        st = -2 * x1 * hx1;
        ct = 2 * x1sq - rsq;
        area -= circseg(rsq, st, ct);
      }
    }
    else
      hx1 = 0;

    if(y1sq < rsq) {
      hy1 = sqrt(rsq - y1sq);
      if(-hy1 > x1 && hy1 < x2) {
        /* r cos t1 = y1 */
        st = -2 * y1 * hy1;
        ct = 2 * y1sq - rsq;
        area -= circseg(rsq, st, ct);
      }
    }
    else
      hy1 = 0;
    
    if(x2sq < rsq) {
      hx2 = sqrt(rsq - x2sq);
      if(-hx2 > y1 && hx2 < y2) {
        /* r cos t1 = x2 */
        st = 2 * x2 * hx2;
        ct = 2 * x2sq - rsq;
        area -= circseg(rsq, st, ct);
      }
    }
    else
      hx2 = 0;
    
    if(y2sq < rsq) {
      hy2 = sqrt(rsq - y2sq);
      if(-hy2 > x1 && hy2 < x2) {
        /* r cos t1 = y2 */
        st = 2 * y2 * hy2;
        ct = 2 * y2sq - rsq;
        area -= circseg(rsq, st, ct);
      }
    }
    else
      hy2 = 0;

    /* Now handle 1-intersection cases: how many vertices are inside
       the circle? */
    if(x1sq+y1sq >= rsq) {  /* vertices entirely outside */
      /* Does pixel contain the origin? */
      if(x1 < 0 && y1 < 0)  /* circle is inside pixel */
        area += rsq*M_PI;

      /* otherwise pixel is outside circle, area is zero */
    }
    else if(x1sq+y2sq >= rsq &&
            x2sq+y1sq >= rsq) {  /* one vertex x1,y1 inside */
      /* Overlap area is circular segment plus triangle.
         r sin t1 = y1, r cos t1 = b
         r cos t2 = x1, r sin t1 = a */
      a = hx1;
      b = hy1;
      st = a*b - x1*y1;  /* sin(t2-t1) */
      ct = b*x1 + a*y1;  /* cos(t2-t1) */
      area += circseg(rsq, st, ct) + 0.5*(a-y1)*(b-x1);
    }
    else if(x2sq+y1sq >= rsq) {  /* two vertices x1,y1 and x1,y2 inside */
      /* Overlap area is circular segment plus rectangle.
         r sin t1 = y1, r cos t1 = a
         r sin t2 = y2, r cos t2 = b */
      a = hy1;
      b = hy2;
      st = y2*a - y1*b;
      ct = a*b + y1*y2;
      area += circseg(rsq, st, ct) + 0.5*(a-x1+b-x1);
    }
    else {  /* three vertices, all except x2,y2 inside */
      /* Overlap area is circular segment plus square minus triangle.
         r cos t1 = x2, r sin t1 = a
         r sin t2 = y2, r cos t2 = b */
      a = hx2;
      b = hy2;
      st = y2*x2 - a*b;
      ct = b*x2 + a*y2;
      area += circseg(rsq, st, ct) + 1.0 - 0.5*(y2-a)*(x2-b);
    }

    return(area);
  }
}
