// $Id$
//==============================================================================
//!
//! \file LagrangeInterpolator.h
//!
//! \date Feb 12 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Utilities for Lagrange interpolation.
//!
//==============================================================================

#ifndef LAGRANGE_INTERPOLATOR_H_
#define LAGRANGE_INTERPOLATOR_H_

#include "matrix.h"


/*!
  \brief Lagrange interpolation in one dimension.
*/

class LagrangeInterpolator
{
public:
  //! \brief The constructor initializes the grid point array.
  explicit LagrangeInterpolator(const std::vector<Real>& g) : grid(g) {}

  //! \brief Evaluate the interpolator at the specified point.
  Real evaluate(Real x, size_t nr) const;
  //! \brief Interpolates the given data at the specified point.
  Real interpolate(Real x, const std::vector<Real>& data) const;
  //! \brief Returns an interpolation matrix for a grid of data points.
  utl::matrix<Real> get(const std::vector<Real>& xg) const;

private:
  std::vector<Real> grid; //!< Grid points
};

#endif
