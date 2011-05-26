// $Id$
//==============================================================================
//!
//! \file GaussQuadrature.h
//!
//! \date Oct 27 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Gaussian quadrature rules in one dimension.
//!
//==============================================================================

#ifndef _GAUSS_QUADRATURE_H
#define _GAUSS_QUADRATURE_H


/*!
  \brief Gaussian quadrature rules in one dimension.
*/

class GaussQuadrature
{
public:
  //! \brief Get gauss point coordinates in the domain [-1,1].
  //! \param[in] n Number of gauss points
  static const double* getCoord (int n) { return getGauss(n,0); }
  //! \brief Get gauss point weights.
  //! \param[in] n Number of gauss points
  static const double* getWeight(int n) { return getGauss(n,1); }

private:
  //! \brief Returns the gauss point coordinates or weights.
  //! \param[in] n Number of gauss points
  //! \param[in] i Option telling what to return (0=coordinates, 1=weights)
  static const double* getGauss(int n, int i);
};

#endif
