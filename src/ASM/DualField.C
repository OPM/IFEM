// $Id$
//==============================================================================
//!
//! \file DualField.C
//!
//! \date Oct 22 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of dual fields for goal-oriented error estimation.
//!
//==============================================================================

#include "DualField.h"
#include "Vec3Oper.h"


DualRealFunc::DualRealFunc (const Vec3& o, const Vec3& n, const Vec3& XZp,
                            double d, double w, size_t p) : X0(o), normal(n)
{
  tangent.cross(XZp-X0,normal);
  normal.normalize();
  tangent.normalize();
  depth = d;
  width = w;
  patch = p;
}


DualRealFunc::DualRealFunc (const Vec3& o, const Vec3& n,
                            double d, double w, size_t p) : X0(o)
{
  normal.x = n.x;
  normal.y = n.y;
  normal.normalize();
  tangent.x = -normal.y;
  tangent.y =  normal.x;
  depth = d;
  width = w;
  patch = p;
}


DualRealFunc::DualRealFunc (const Vec3& o, const Vec3Pair& d,
                            size_t p) : X0(o), Xll(d.first), Xur(d.second)
{
  normal.x = tangent.y = 1.0;
  depth = width = 0.0;
  patch = p;
}


bool DualRealFunc::inDomain (const Vec3& X) const
{
  if (depth > 0.0)
  {
    double x = (X-X0)*normal;
    if (x > 0.0 || -x >= depth)
      return false;
    else if (width <= 0.0)
      return true;

    double y = (X-X0)*tangent;
    return y <= 0.5*width && y >= -0.5*width;
  }

  // Rectangular point domain
  return (X.x >= Xll.x && X.x <= Xur.x) &&
         (X.y >= Xll.y && X.y <= Xur.y) &&
         (X.z >= Xll.z && X.z <= Xur.z);
}


double DualRealFunc::value (const Vec3& X, bool ignoreDomain) const
{
  if (ignoreDomain || this->inDomain(X))
    return depth > 0.0 ? 1.0 + (X-X0)*normal/depth : 1.0;

  return 0.0; // outside the function domain
}


DualVecFunc::DualVecFunc (int c, const Vec3& o, const Vec3& n, const Vec3& XZp,
                          double d, double w, size_t p)
  : VecFunc(3), W(o,n,XZp,d,w,p)
{
  comp = c;
}


DualVecFunc::DualVecFunc (int c, const Vec3& o, const Vec3& n,
                          double d, double w, size_t p)
  : VecFunc(2), W(o,n,d,w,p)
{
  comp = c > 2 ? 6 : c;
}


DualVecFunc::DualVecFunc (int c, const Vec3& o, const Vec3Pair& d, size_t p)
  : VecFunc(2), W(o,d,p)
{
  comp = c;
}


/*!
  The field does not evaluate to zero outside the defined function domain, since
  it is used to find the control point values of a spline representation of it.
*/

Vec3 DualVecFunc::evaluate (const Vec3& X) const
{
  double w = W.value(X,true);
  if (fabs(w) < 1.0e-12)
    return Vec3();

  if (W.isPointExtraction())
    switch (comp) {
    case 1: return Vec3(W.ecc(X,1),0.0,0.0);
    case 2: return Vec3(0.0,W.ecc(X,2),0.0);
    case 3: return Vec3(0.0,0.0,W.ecc(X,3));
    case 4: return Vec3(0.0,0.5*W.ecc(X,3),0.5*W.ecc(X,2));
    case 5: return Vec3(0.5*W.ecc(X,3),0.0,0.5*W.ecc(X,1));
    case 6: return Vec3(0.5*W.ecc(X,2),0.5*W.ecc(X,1),0.0);
    }
  else
    switch (comp) {
    case 1: return  w*W.x();
    case 2: return  w*W.y();
    case 3: return  w*W.z();
    case 4: return -w*W.ecc(X,3)*W.y() + w*W.ecc(X,2)*W.z();
    case 5: return  w*W.ecc(X,3)*W.x() - w*W.ecc(X,1)*W.z();
    case 6: return -w*W.ecc(X,2)*W.x() + w*W.ecc(X,1)*W.y();
    }

  return Vec3();
}
