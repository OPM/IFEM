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
  Vec3() : x(v[0]), y(v[1]), z(v[2]) { x = y = z = Real(0); }
  //! \brief Constructor creating a point at the specified location.
  Vec3(Real X, Real Y, Real Z = Real(0)) : x(v[0]), y(v[1]), z(v[2])
  { x = X; y = Y; z = Z; }
  //! \brief Constructor creating a point at the specified location.
  Vec3(const Real* pos, size_t n = 3) : x(v[0]), y(v[1]), z(v[2])
  { for (size_t i = 0; i < 3; i++) v[i] = i < n ? pos[i] : Real(0); }
  //! \brief Constructor creating a point from the given \a std::vector.
  Vec3(const std::vector<Real>& X) : x(v[0]), y(v[1]), z(v[2])
  { for (size_t i = 0; i < 3; i++) v[i] = i < X.size() ? X[i] : Real(0); }
  //! \brief Copy constructor.
  Vec3(const Vec3& X) : x(v[0]), y(v[1]), z(v[2])
  { x = X.x; y = X.y; z = X.z; }
  //! \brief Constructor creating a cross product of two other vectors.
  Vec3(const Vec3& X, const Vec3& Y) : x(v[0]), y(v[1]), z(v[2])
  { this->cross(X,Y); }

  //! \brief Empty destructor.
  virtual ~Vec3() {}

  //! \brief Assignment operator.
  Vec3& operator=(const Vec3& X) { x = X.x; y = X.y; z = X.z; return *this; }
  //! \brief Overloaded assignment operator.
  Vec3& operator=(Real val) { x = y = z = val; return *this; }

  //! \brief Indexing operator for component reference (1-based).
  const Real& operator()(int i) const { return v[i-1]; }
  //! \brief Indexing operator for component reference (0-based).
  const Real& operator[](int i) const { return v[i]; }
  //! \brief Indexing operator for component assignment (1-based).
  Real& operator()(int i) { return v[i-1]; }
  //! \brief Indexing operator for component assignment (0-based).
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

  //! \brief Return the square of the length of the vector.
  Real length2() const { return x*x+y*y+z*z; }

  //! \brief Return the length of the vector.
  double length() const { return sqrt(this->length2()); }

  //! \brief Normalize the vector and return its length.
  Real normalize(double tol = 1.0e-16)
  {
    double len = this->length();
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
      else if (v[i]-tol > a.v[i])
        return false;

    return false;
  }

  //! \brief Check if a vector is inside the box defined by two other vectors.
  bool inside(const Vec3& a, const Vec3& b, double tol = 1.0e-6) const
  {
    for (int i = 0; i < 3; i++)
      if (v[i]+tol < a.v[i] || v[i]-tol > b.v[i])
        return false;

    return true;
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
  const Real* u; //!< Spline parameters of point

  Real t;   //!< The time coordinate
  int  idx; //!< Nodal point index

  //! \brief Default constructor creating a point at origin.
  Vec4(const Real* par = nullptr, Real T = Real(0))
  {
    u = par;
    t = T;
    idx = -1;
  }

  //! \brief Constructor creating a point at the specified location.
  Vec4(Real X, Real Y, Real Z, Real T = Real(0)) : Vec3(X,Y,Z)
  {
    u = nullptr;
    t = T;
    idx = -1;
  }

  //! \brief Constructor creating a point at the specified location.
  Vec4(const Vec3& X, Real T = Real(0), int id = -1) : Vec3(X)
  {
    u = nullptr;
    t = T;
    idx = id;
  }

  //! \brief Constructor creating a point at the specified location.
  Vec4(const Vec3& X, Real T, const Real* par) : Vec3(X)
  {
    u = par;
    t = T;
    idx = -1;
  }

  //! \brief Constructor creating a point from the given \a std::vector.
  Vec4(const std::vector<Real>& X, const Real* par = nullptr) : Vec3(X)
  {
    u = par;
    t = X.size() > 3 ? X[3] : Real(0);
    idx = -1;
  }

  //! \brief Copy constructor.
  Vec4(const Vec4& X) : Vec3(X)
  {
    u = X.u;
    t = X.t;
    idx = X.idx;
  }

  //! \brief Assignment operator.
  Vec4& operator=(const Vec4& X)
  {
    x = X.x;
    y = X.y;
    z = X.z;
    t = X.t;
    idx = X.idx;
    u = X.u;

    return *this;
  }

  //! \brief Overloaded assignment operator.
  Vec4& operator=(Real val) { x = y = z = val; return *this; }

  //! \brief Assignment method.
  Vec4& assign(const Vec3& X)
  {
    const Vec4* x4 = dynamic_cast<const Vec4*>(&X);
    if (x4)
      this->operator=(*x4);
    else
    {
      x = X.x;
      y = X.y;
      z = X.z;
    }

    return *this;
  }

  //! \brief Print out a vector.
  virtual std::ostream& print(std::ostream& os, double tol = 1.0e-12) const
  {
    if (idx >= 0) os <<"(idx="<< idx <<") ";
    this->Vec3::print(os,tol);
    if (u) os <<" ("<< u[0] <<" "<< u[1] <<" "<< u[2] <<")";
    if (t > Real(0)) os <<" "<< t;
    return os;
  }
};


using Vec3Pair    = std::pair<Vec3,Vec3>;    //!< A pair of two point vectors
using PointValue  = std::pair<Vec3,double>;  //!< A point with associated value
using Vec3Vec     = std::vector<Vec3>;       //!< An array of point vectors
using PointValues = std::vector<PointValue>; //!< An array of point values

#endif
