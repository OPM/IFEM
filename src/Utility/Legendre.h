// $Id$
//==============================================================================
//!
//! \file Legendre.h
//!
//! \date Mar 19 2009
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Various utility methods for Spectral elements.
//!
//==============================================================================

#ifndef _LEGENDRE_H
#define _LEGENDRE_H

#include "MatVec.h"


namespace Legendre
{
  //! \brief Get Gauss-Legendre points and weights in the domain [-1,1].
  //! \param[out] weights Computed Gauss-Legendre weight
  //! \param[out] points Computed Gauss-Legendre points	
  //! \param[in] n Number of Gauss points
  bool GL(RealArray& weights, RealArray& points, int n);

  //! \brief Get Gauss-Lobatto-Legendre points and weights in the domain [-1,1].
  //! \param[out] weights Computed Gauss-Legendre weight
  //! \param[out] points Computed Gauss-Legendre points	
  //! \param[in] n Number of Gauss points
  bool GLL(Vector& weights, Vector& points, int n);

  //! \brief Evaluates the \a n-th Legendre polynomial.
  //! \param[in] n Polynomial degree
  //! \param[in] x Evaluation point
  //! \param[out] retval Polynomial value at point \a x
  bool LegendreEval(int n, Real x, Real& retval);

  //! \brief Evaluates the first derivative of the \a n-th Legendre polynomial.
  //! \param[in] n Polynomial degree
  //! \param[in] x Evaluation point
  //! \param[out] retval Value of the first derivative at point \a x
  bool LegendreDerEval(int n, Real x, Real& retval);

  //! \brief Evaluates first derivatives of the \a n one-dimensional
  //! Lagrange interpolation polynomials through \a n GLL-points.
  //! \param[in] n Number of GLL points\polynomials
  //! \param[out] der Evaluated values
  bool basisDerivatives(int n, Matrix& der);
}

#endif
