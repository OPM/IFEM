// $Id$
//==============================================================================
//!
//! \file Vec3Oper.C
//!
//! \date Dec 19 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Global algebraic operations involving Vec3 objects.
//!
//==============================================================================

#include "Vec3Oper.h"
#include "Vec3.h"


double Vec3::comparisonTolerance = 1.0e-4;


/*!
  \brief Multiplication of a matrix and a vector.
  \return \f$ {\bf y} = {\bf A} {\bf x} \f$

  Special version for computation of 3D point vectors.
  The number of rows in \b A must be (at least) 3.
  Extra rows (if any) in \b A are ignored.
*/

Vec3 operator* (const utl::matrix<real>& A, const std::vector<real>& x)
{
  std::vector<real> y(A.rows(),0.0);
  A.multiply(x,y);
  return Vec3(y);
}


Vec3 operator* (const Vec3& a, real value)
{
  return Vec3(a.x*value, a.y*value, a.z*value);
}


Vec3 operator* (real value, const Vec3& a)
{
  return Vec3(a.x*value, a.y*value, a.z*value);
}


Vec3 operator/ (const Vec3& a, real value)
{
  return Vec3(a.x/value, a.y/value, a.z/value);
}


real operator* (const Vec3& a, const Vec3& b)
{
  return a.x*b.x + a.y*b.y + a.z*b.z;
}


Vec3 operator+ (const Vec3& a, const Vec3& b)
{
  return Vec3(a.x+b.x, a.y+b.y, a.z+b.z);
}


Vec3 operator- (const Vec3& a, const Vec3& b)
{
  return Vec3(a.x-b.x, a.y-b.y, a.z-b.z);
}


bool operator== (const Vec3& a, const Vec3& b)
{
  return a.equal(b,Vec3::comparisonTolerance);
}


bool operator!= (const Vec3& a, const Vec3& b)
{
  return !a.equal(b,Vec3::comparisonTolerance);
}


bool operator< (const Vec3& a, const Vec3& b)
{
  return a.lessThan(b,Vec3::comparisonTolerance);
}


std::ostream& operator<< (std::ostream& os, const Vec3& a)
{
  return a.print(os);
}
