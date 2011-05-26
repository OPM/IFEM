// $Id$
//==============================================================================
//!
//! \file Lagrange.h
//!
//! \date Feb 10 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Evaluation of Lagrange basis functions.
//!
//==============================================================================

#ifndef _LAGRANGE_H
#define _LAGRANGE_H

#include "MatVec.h"


/*!
  \brief Evaluation of Lagrange basis functions.
*/

class Lagrange
{
public:
  //! \brief Constructor initializing the reference to natural coordinates.
  //! \param[in] p Natural interpolation point coordinates in range [-1,1]
  Lagrange(const RealArray& p) : points(p) {}

  //! \brief Evaluates a 1D Lagrange polynomial.
  //! \param[in] polnum Which polynomial of the basis to evaluate
  //! \param[in] xi Natural coordinate of the evaluation point
  //! \param[out] retval The computed polynomial value
  bool evalPol(int polnum, double xi, double& retval) const;
  //! \brief Evaluates the first derivative of a 1D Lagrange polynomial.
  //! \param[in] polnum Which polynomial of the basis to evaluate
  //! \param[in] xi Natural coordinate of the evaluation point
  //! \param[out] retval The computed polynomial derivative
  bool evalDer(int polnum, double xi, double& retval) const;

  //! \brief Evaluates a 1D, 2D or 3D Lagrangian basis at a given point.
  //! \param[out] val Values of all basis functions
  //! \param[out] derval Derivatives of all basis functions
  //! \param[in] p1 Polynomial degree in first parameter direction
  //! \param[in] x1 Natural coordinate in first parameter direction
  //! \param[in] p2 Polynomial degree in second parameter direction
  //! \param[in] x2 Natural coordinate in second parameter direction
  //! \param[in] p3 Polynomial degree in third parameter direction
  //! \param[in] x3 Natural coordinate in third parameter direction
  //!
  //! \details If derval is a 0-pointer, the derivatives are not computed.
  //! If \a p2 is zero, a 1D basis is assumed. Otherwise,
  //! if \a p3 is zero, a 2D basis is assumed. If \a p1, \a p2, and \a p3
  //! all are non-zero, a 3D basis is assumed.
  static bool computeBasis (Vector& val, Matrix& derval,
			    int p1 = 0, double x1 = 0.0,
			    int p2 = 0, double x2 = 0.0,
			    int p3 = 0, double x3 = 0.0);

  //! \brief Evaluates a 1D, 2D or 3D Lagrangian basis at a given point.
  //! \param[out] val Values of all basis functions
  //! \param[out] derval Derivatives of all basis functions
  //! \param[in] p1 Natural point coordinates in first parameter direction
  //! \param[in] x1 Natural coordinate in first parameter direction
  //! \param[in] p2 Natural point coordinates in second parameter direction
  //! \param[in] x2 Natural coordinate in second parameter direction
  //! \param[in] p3 Natural point coordinates in third parameter direction
  //! \param[in] x3 Natural coordinate in third parameter direction
  //!
  //! \details If derval is a 0-pointer, the derivatives are not computed.
  //! If \a p2 is empty, a 1D basis is assumed. Otherwise,
  //! if \a p3 is empty, a 2D basis is assumed. If \a p1, \a p2, and \a p3
  //! all are non-empty, a 3D basis is assumed.
  static bool computeBasis (Vector& val, Matrix& derval,
			    const RealArray& p1, double x1,
			    const RealArray& p2 = RealArray(), double x2 = 0.0,
			    const RealArray& p3 = RealArray(), double x3 = 0.0);

private:
  const RealArray& points; //!< Natural coordinates of the interpolation points
};

#endif
