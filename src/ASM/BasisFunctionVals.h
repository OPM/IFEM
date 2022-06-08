// $Id$
//==============================================================================
//!
//! \file BasisFunctionVals.h
//!
//! \date Jun 1 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Basis function values container.
//!
//==============================================================================

#ifndef _BASIS_FUNCTION_VALS_H
#define _BASIS_FUNCTION_VALS_H

#include "MatVec.h"


/*!
  \brief Struct holding basis function values and derivatives.
*/

struct BasisFunctionVals
{
  Vector     N;    //!< Basis function values
  Matrix    dNdu;  //!< Basis function derivatives
  Matrix3D d2Ndu2; //!< Second order derivatives of basis functions
  Matrix4D d3Ndu3; //!< Third order derivatives of basis functions
};

#endif
