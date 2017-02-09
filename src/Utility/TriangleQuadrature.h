// $Id$
//==============================================================================
//!
//! \file TriangleQuadrature.h
//!
//! \date Feb 07 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Triangle quadrature rules.
//!
//==============================================================================

#ifndef _TRIANGLE_QUADRATURE_H
#define _TRIANGLE_QUADRATURE_H


/*!
  \brief Triangle quadrature rules.
*/

namespace TriangleQuadrature
{
  //! \brief Get quadrature point (area-) coordinates in the domain [0,1].
  //! \param[in] n Number of quadrature points
  const double* getCoord(int n);
  //! \brief Get quadrature point weights.
  //! \param[in] n Number of quadrature points
  const double* getWeight(int n);
}

#endif
