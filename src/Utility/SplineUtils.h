// $Id$
//==============================================================================
//!
//! \file SplineUtils.h
//!
//! \date Mar 29 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Various utility functions on spline objects.
//!
//==============================================================================

#ifndef SPLINE_UTILS_H
#define SPLINE_UTILS_H

#include "Function.h"

class Vec4;

namespace Go {
  class Point;
  class SplineCurve;
  class SplineSurface;
}


namespace SplineUtils //! Various utility functions on spline objects.
{
  //! \brief Helper method for casting a \a Go::Point object to Vec3.
  Vec3 toVec3(const Go::Point& X, int nsd = 3);
  //! \brief Helper method for casting a \a Go::Point and time object to Vec4.
  Vec4 toVec4(const Go::Point& X, real time);

  //! \brief Projects a scalar-valued function onto a spline curve.
  Go::SplineCurve* project(const Go::SplineCurve* curve,
			   const RealFunc& f, real time = real(0));
  //! \brief Projects a vector-valued function onto a spline curve.
  Go::SplineCurve* project(const Go::SplineCurve* curve,
			   const VecFunc& f, int nComp, real time = real(0));

  //! \brief Projects a scalar-valued function onto a spline surface.
  Go::SplineSurface* project(const Go::SplineSurface* surface,
			     const RealFunc& f, real time = real(0));
  //! \brief Projects a vector-valued function onto a spline surface.
  Go::SplineSurface* project(const Go::SplineSurface* surface,
			     const VecFunc& f, int nComp, real time = real(0));
}

#endif
