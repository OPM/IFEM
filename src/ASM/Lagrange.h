// $Id$
//==============================================================================
//!
//! \file Lagrange.h
//!
//! \date Feb 10 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Evaluation of %Lagrange basis functions.
//!
//==============================================================================

#ifndef _LAGRANGE_H
#define _LAGRANGE_H

#include "MatVec.h"


namespace Lagrange //! Evaluation of %Lagrange basis functions.
{
  //! \brief Evaluates a 1D, 2D or 3D Lagrangian basis at a given point.
  //! \param[out] val Values of all basis functions
  //! \param[in] p1 Polynomial degree in first parameter direction
  //! \param[in] x1 Natural coordinate in first parameter direction
  //! \param[in] p2 Polynomial degree in second parameter direction
  //! \param[in] x2 Natural coordinate in second parameter direction
  //! \param[in] p3 Polynomial degree in third parameter direction
  //! \param[in] x3 Natural coordinate in third parameter direction
  //!
  //! \details If \a p2 is zero, a 1D basis is assumed. Otherwise,
  //! if \a p3 is zero, a 2D basis is assumed. If \a p1, \a p2, and \a p3
  //! all are non-zero, a 3D basis is assumed.
  bool computeBasis(RealArray& val, int p1, double x1,
                    int p2 = 0, double x2 = 0.0,
                    int p3 = 0, double x3 = 0.0);

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
  //! \details If \a p2 is zero, a 1D basis is assumed. Otherwise,
  //! if \a p3 is zero, a 2D basis is assumed. If \a p1, \a p2, and \a p3
  //! all are non-zero, a 3D basis is assumed.
  bool computeBasis(RealArray& val, Matrix& derval, int p1, double x1,
                    int p2 = 0, double x2 = 0.0,
                    int p3 = 0, double x3 = 0.0);

  //! \brief Evaluates a 1D, 2D or 3D Lagrangian basis at a given point.
  //! \param[out] val Values of all basis functions
  //! \param[out] derval Pointer to derivatives of all basis functions
  //! \param[in] p1 Natural point coordinates in first parameter direction
  //! \param[in] x1 Natural coordinate in first parameter direction
  //! \param[in] p2 Natural point coordinates in second parameter direction
  //! \param[in] x2 Natural coordinate in second parameter direction
  //! \param[in] p3 Natural point coordinates in third parameter direction
  //! \param[in] x3 Natural coordinate in third parameter direction
  //!
  //! \details If \a derval is a null pointer, the derivatives are not computed.
  //! If \a p2 is empty, a 1D basis is assumed. Otherwise,
  //! if \a p3 is empty, a 2D basis is assumed. If \a p1, \a p2, and \a p3
  //! all are non-empty, a 3D basis is assumed.
  bool computeBasis(RealArray& val, Matrix* derval,
                    const RealArray& p1, double x1,
                    const RealArray& p2 = {}, double x2 = 0.0,
                    const RealArray& p3 = {}, double x3 = 0.0);
};

#endif
