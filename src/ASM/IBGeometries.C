
// $Id$
//==============================================================================
//!
//! \file IBGeometries.C
//!
//! \date Jan 24 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Physical geometries for immersed boundary simulations.
//!
//==============================================================================

#include "IBGeometries.h"
#include <cstddef>


double Hole2D::Alpha (double X, double Y, double) const
{
  // Determine radius of the point to be checked
  double r2 = (X-Xc)*(X-Xc) + (Y-Yc)*(Y-Yc);

  // Determine if point is located within the hole or not
  double alpha = r2 > R*R ? 1.0 : 0.0;

  // Alpha is the penalization parameter in the sense of equation (29) in the
  // immersed boundary paper. It can be used either to indicate where a point
  // is located, or can be used as a penalization value for the local stiffness
  // contribution at the current integration point.
  // In the latter case, it needs to have data type double.
  return alpha;
}


double PerforatedPlate2D::Alpha (double X, double Y, double) const
{
  // Determine if point is located within any of the holes or not
  double alpha = 1.0;
  for (size_t i = 0; i < holes.size() && alpha > 0.0; i++)
    alpha = holes[i].Alpha(X,Y);

  return alpha;
}
