// $Id$
//==============================================================================
//!
//! \file LagrangeInterpolator.C
//!
//! \date Feb 12 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Utilities for Lagrange interpolation.
//!
//==============================================================================

#include "LagrangeInterpolator.h"


Real LagrangeInterpolator::evaluate (Real x, size_t nr) const
{
  Real result = Real(1);
  for (size_t j = 0; j < grid.size(); j++)
    if (j != nr)
      result *= (x-grid[j])/(grid[nr]-grid[j]);

  return result;
}


Real LagrangeInterpolator::interpolate (Real x,
                                        const std::vector<Real>& data) const
{
  Real result = Real(0);
  for (size_t i = 0; i < data.size(); i++)
    result += data[i]*this->evaluate(x,i);

  return result;
}


utl::matrix<Real> LagrangeInterpolator::get (const std::vector<Real>& xg) const
{
  utl::matrix<Real> result(xg.size(),grid.size());
  for (size_t i = 0; i < grid.size(); i++)
    for (size_t j = 0; j < xg.size(); j++)
      result(j+1,i+1) = this->evaluate(xg[j],i);

  return result;
}
