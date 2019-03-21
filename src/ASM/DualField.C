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
                            double d, double w)
{
  X0 = o;
  normal = n;
  tangent.cross(XZp-X0,normal);
  normal.normalize();
  tangent.normalize();
  depth = d;
  width = w;
}


DualRealFunc::DualRealFunc (const Vec3& o, const Vec3& n, double d, double w)
{
  X0 = o;
  normal = n;
  normal.z = 0.0;
  normal.normalize();
  tangent.x = -normal.y;
  tangent.y = -normal.x;
  tangent.z = 0.0;
  depth = d;
  width = w;
}


bool DualRealFunc::inDomain (const Vec3& X) const
{
  double x = (X-X0)*normal;
  if (x > 0.0 || -x >= depth)
    return false;

  if (width > 0.0)
  {
    double y = (X-X0)*tangent;
    return y <= 0.5*width && y >= -0.5*width;
  }

  return true;
}


double DualRealFunc::value (const Vec3& X, bool ignoreDomain) const
{
  if (ignoreDomain || this->inDomain(X))
    return 1.0 + (X-X0)*normal/depth;

  return 0.0; // outside the function domain
}


DualVecFunc::DualVecFunc (int c, const Vec3& o, const Vec3& n, const Vec3& XZp,
                          double d, double w) : VecFunc(3), W(o,n,XZp,d,w)
{
  comp = c;
}


DualVecFunc::DualVecFunc (int c, const Vec3& o, const Vec3& n,
                          double d, double w) : VecFunc(2), W(o,n,d,w)
{
  comp = c > 2 ? 6 : c;
}


/*!
  The field does not evaluate to zero outside the defined function domain, since
  it is used to find the control point values of a spline representation of it.
*/

Vec3 DualVecFunc::evaluate (const Vec3& X) const
{
  double w = W.value(X,true);
  if (w != 0.0)
    switch (comp) {
    case 1: return Vec3( w           , 0.0         , 0.0         );
    case 2: return Vec3( 0.0         , w           , 0.0         );
    case 3: return Vec3( 0.0         , 0.0         , w           );
    case 4: return Vec3( 0.0         ,-w*W.ecc(X,3), w*W.ecc(X,2));
    case 5: return Vec3( w*W.ecc(X,3), 0.0         ,-w*W.ecc(X,1));
    case 6: return Vec3(-w*W.ecc(X,2), w*W.ecc(X,1), 0.0         );
    }

  return Vec3();
}
