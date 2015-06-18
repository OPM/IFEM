// $Id$
//==============================================================================
//!
//! \file Vec3.h
//!
//! \date Dec 10 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of a point in 3D space with some basic operations.
//!
//==============================================================================

#ifndef _VEC3_H
#define _VEC3_H

#include <vector>
#include <iostream>
#include <cmath>


/*!
  \brief Simple class for representing a point in 3D space.
*/

class Vec3
{
  Real v[3]; //!< The actual point vector

public:
  Real& x; //!< Reference to X-component
  Real& y; //!< Reference to Y-component
  Real& z; //!< Reference to Z-component

  //! \brief Default constructor creating a point at origin.
  Vec3() : x(v[0]), y(v[1]), z(v[2]) { x = y = z = 0.0; }
  //! \brief Constructor creating a point at the specified location.
  Vec3(Real X, Real Y, Real Z) : x(v[0]), y(v[1]), z(v[2])
  { x = X; y = Y; z = Z; }
  //! \brief Constructor creating a point at the specified location.
  Vec3(const Real* pos, size_t n = 3) : x(v[0]), y(v[1]), z(v[2])
  { for (size_t i = 0; i < 3; i++) v[i] = i < n ? pos[i] : 0.0; }
  //! \brief Constructor creating a point from the given \a std::vector.
  Vec3(const std::vector<Real>& X) : x(v[0]), y(v[1]), z(v[2])
  { for (size_t i = 0; i < 3; i++) v[i] = i < X.size() ? X[i] : 0.0; }
  //! \brief Copy constructor.
  Vec3(const Vec3& X) : x(v[0]), y(v[1]), z(v[2])
  { x = X.x; y = X.y; z = X.z; }
  //! \brief Constructor creating a cross product of two other vectors.
  Vec3(const Vec3& X, const Vec3& Y) : x(v[0]), y(v[1]), z(v[2])
  { this->cross(X,Y); }

  //! \brief Assignment operator.
  Vec3& operator=(const Vec3& X) { x = X.x; y = X.y; z = X.z; return *this; }
  //! \brief Overloaded assignment operator.
  Vec3& operator=(Real val) { x = y = z = val; return *this; }

  //! \brief Indexing operator for component reference (1-based).
  const Real& operator()(int i) const { return v[i-1]; }
  //! \brief Indexing operator for component reference (0-based).
  const Real& operator[](int i) const { return v[i]; }
  //! \brief Indexing operator for component assignment.
  Real& operator()(int i) { return v[i-1]; }
  //! \brief Indexing operator for component assignment.
  Real& operator[](int i) { return v[i]; }

  //! \brief Reference through a pointer.
  const Real* ptr() const { return v; }

  //! \brief Conversion to std::vector.
  std::vector<Real> vec(size_t n = 3) const
  { return std::vector<Real>(v, v + (n < 3 ? n : 3)); }

  //! \brief Multiplication with a scalar.
  Vec3& operator*=(Real c) { x *= c; y *= c; z *= c; return *this; }
  //! \brief Division by a scalar.
  Vec3& operator/=(Real d) { x /= d; y /= d; z /= d; return *this; }

  //! \brief Add the given vector \b X to \a *this.
  Vec3& operator+=(const Vec3& X)
  {
    x += X.x;
    y += X.y;
    z += X.z;
    return *this;
  }

  //! \brief Subtract the given vector \b X from \a *this.
  Vec3& operator-=(const Vec3& X)
  {
    x -= X.x;
    y -= X.y;
    z -= X.z;
    return *this;
  }

  //! \brief Return the sum of the vector.
  Real sum() const { return x+y+z; }

  //! \brief Return the length of the vector.
  Real length() const { return sqrt(length2()); }

  //! \brief Return the square of the length of the vector.
  Real length2() const { return x*x+y*y+z*z; }

  //! \brief Normalize the vector and return its length.
  Real normalize(double tol = 1.0e-16)
  {
    double len = sqrt(x*x+y*y+z*z);
    if (len <= tol) return len;
    x /= len; y /= len; z /= len;
    return len;
  }

  //! \brief Cross product between two vectors.
  Vec3& cross(const Vec3& a, const Vec3& b)
  {
    x = a.y*b.z - a.z*b.y;
    y = a.z*b.x - a.x*b.z;
    z = a.x*b.y - a.y*b.x;
    return *this;
  }

  //! \brief Equality check between two vectors.
  bool equal(const Vec3& a, double tol = 1.0e-6) const
  {
    if (fabs(x-a.x) <= tol)
      if (fabs(y-a.y) <= tol)
        if (fabs(z-a.z) <= tol)
          return true;

    return false;
  }

  //! \brief Check if the vector is zero with the given tolerance.
  bool isZero(double tol = 1.0e-6) const
  {
    if (fabs(x) <= tol)
      if (fabs(y) <= tol)
        if (fabs(z) <= tol)
          return true;

    return false;
  }

  //! \brief Check if one vector is less than the other.
  //! \details First we compare the z-coordinates. Only if the z-coordinates are
  //! equal, we compare the y-coordinates, and if they are equal too,
  //! we finally compare the x-coordinates.
  bool lessThan(const Vec3& a, double tol = 1.0e-6) const
  {
    for (int i = 2; i >= 0; i--)
      if (v[i]+tol < a.v[i])
        return true;
      else if (v[i] > a.v[i]+tol)
        return false;

    return false;
  }

  //! \brief Print out a vector.
  virtual std::ostream& print(std::ostream& os, double tol = 1.0e-12) const
  {
    if (fabs(x) > tol)
      os << x;
    else
      os <<"0";

    if (fabs(y) > tol)
      os <<" "<< y;
    else
      os <<" 0";

    if (fabs(z) > tol)
      os <<" "<< z;
    else
      os <<" 0";

    return os;
  }

  static double comparisonTolerance; //!< Coordinate comparison tolerance
};


/*!
  \brief Simple class for representing a point in 3D space and time.
*/

class Vec4 : public Vec3
{
public:
  Real t;   //!< The time coordinate
  int  idx; //!< Nodal point index

  //! \brief Default constructor creating a point at origin.
  Vec4() { t = 0.0; idx = -1; }
  //! \brief Constructor creating a point at the specified location.
  Vec4(Real X, Real Y, Real Z, Real T = 0.0) : Vec3(X,Y,Z) { t = T; idx = -1; }
  //! \brief Constructor creating a point at the specified location.
  Vec4(const Vec3& X, Real T = 0.0, int id = -1) : Vec3(X) { t = T; idx = id; }
  //! \brief Constructor creating a point from the given \a std::vector.
  Vec4(const std::vector<Real>& X) : Vec3(X)
  { t = X.size() > 3 ? X[3] : 0.0; idx = -1; }
  //! \brief Copy constructor.
  Vec4(const Vec4& X) : Vec3(X) { t = X.t; idx = X.idx; }

  //! \brief Assignment operator.
  Vec4& operator=(const Vec4& X)
  { x = X.x; y = X.y; z = X.z; t = X.t; idx = X.idx; return *this; }
  //! \brief Overloaded assignment operator.
  Vec4& operator=(Real val) { x = y = z = val; return *this; }

  //! \brief Assignment method.
  Vec4& assign(const Vec3& X)
  {
    x = X.x; y = X.y; z = X.z;
    const Vec4* x4 = dynamic_cast<const Vec4*>(&X);
    if (x4) { t = x4->t; idx = x4->idx; }
    return *this;
  }

  //! \brief Print out a vector.
  virtual std::ostream& print(std::ostream& os, double tol = 1.0e-12) const
  {
    if (idx >= 0) os <<"(idx="<< idx <<") ";
    this->Vec3::print(os,tol);
    if (t > 0.0) os <<" "<< t;
    return os;
  }
};


typedef std::pair<Vec3,Vec3> Vec3Pair; //!< A pair of two point vectors
typedef std::vector<Vec3>    Vec3Vec;  //!< An array of point vectors

#endif
