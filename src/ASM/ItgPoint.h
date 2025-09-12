// $Id$
//==============================================================================
//!
//! \file ItgPoint.h
//!
//! \date Sep 30 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integration point representation.
//!
//==============================================================================

#ifndef _ITG_POINT_H
#define _ITG_POINT_H

#include <cstddef>


/*!
  \brief Class representing an integration point.
*/

class ItgPoint
{
public:
  //! \brief Default constructor.
  explicit ItgPoint(size_t i = 0)
  {
    iGP = i;
    iel = -1;
    idx = 0;
    u = v = w = xi = eta = zeta = 0.0;
  }

  //! \brief Constructor initializing the spline domain parameters.
  explicit ItgPoint(double a, double b = 0.0, double c = 0.0, size_t i = 0)
  {
    u = a;
    v = b;
    w = c;

    iGP = i;
    idx = 0;
    iel = -1;
    xi = eta = zeta = 0.0;
  }

  //! \brief Alternative constructor initializing the spline domain parameters.
  explicit ItgPoint(const double* par, size_t i = 0)
  {
    u = par[0];
    v = par[1];
    w = par[2];

    iGP = i;
    idx = 0;
    iel = -1;
    xi = eta = zeta = 0.0;
  }

  //! \brief Empty destructor.
  virtual ~ItgPoint() {}

  size_t iGP; //!< Global integration point counter

  double u; //!< First spline parameter of the point
  double v; //!< Second spline parameter of the point
  double w; //!< Third spline parameter of the point

  int    iel;  //!< Identifier of the element containing this point
  size_t idx;  //!< Global index (0-based) of the element containing this point
  double xi;   //!< First local coordinate within current element
  double eta;  //!< Second local coordinate within current element
  double zeta; //!< Third local coordinate within current element
};

#endif
