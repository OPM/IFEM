// $Id$
//==============================================================================
//!
//! \file Point.h
//!
//! \date Sep 19 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of points in 3D space with coordinates and parameters.
//!
//==============================================================================

#ifndef UTL_POINT_H
#define UTL_POINT_H

#include "Vec3.h"


namespace utl
{
  /*!
    \brief Class for representing points in 3D space.
    \details The only difference between this class and the parent class Vec4
    is that it keeps the spline parameters in an internal array.
  */

  class Point : public Vec4
  {
    Real par[3]; //!< Internal storage of spline parameters

  public:
    //! \brief Default constructor.
    Point() : Vec4(par) { par[0] = par[1] = par[2] = Real(0); }

    //! \brief Constructor creating a point at the specified location.
    Point(const Vec3& X, const std::vector<Real>& U) : Vec4(X,Real(0),par)
    {
      for (size_t i = 0; i < 3; i++)
        par[i] = i < U.size() ? U[i] : 0.0;
    }

    //! \brief Copy constructor.
    Point(const Point& X) : Vec4(X,X.t,par)
    {
      for (size_t i = 0; i < 3; i++)
        par[i] = X.par[i];
    }

    //! \brief Constructor creating a point from the given \a std::vector.
    Point(const std::vector<Real>& X) : Vec4(X,par)
    {
      for (size_t i = 0; i < 3; i++)
        par[i] = 3+i < X.size() ? X[3+i] : 0.0;
    }

    //! \brief Assignment operator.
    Point& operator=(const Vec4& X)
    {
      x = X.x;
      y = X.y;
      z = X.z;
      t = X.t;
      idx = X.idx;
      if (X.u)
        for (size_t i = 0; i < 3; i++)
          par[i] = X.u[i];

      return *this;
    }
  };
}


//! \brief Multiplication of a point and a scalar.
inline utl::Point operator*(const utl::Point& a, Real value)
{
  return utl::Point({ value*a.x   , value*a.y   , value*a.z,
                      value*a.u[0], value*a.u[1], value*a.u[2] });
}


//! \brief Multiplication of a scalar and a point.
inline utl::Point operator*(Real value, const utl::Point& a)
{
  return utl::Point({ value*a.x   , value*a.y   , value*a.z,
                      value*a.u[0], value*a.u[1], value*a.u[2] });
}


//! \brief Summation of two points.
inline utl::Point operator+(const utl::Point& a, const utl::Point& b)
{
  return utl::Point({ a.x   +b.x   , a.y   +b.y   , a.z   +b.z,
                      a.u[0]+b.u[0], a.u[1]+b.u[1], a.u[2]+b.u[2] });
}


typedef std::vector<utl::Point> PointVec; //!< An array of points

#endif
