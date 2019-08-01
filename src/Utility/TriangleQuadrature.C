// $Id$
//==============================================================================
//!
//! \file TriangleQuadrature.C
//!
//! \date Feb 07 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Triangle quadrature rules.
//!
//==============================================================================

#include "TriangleQuadrature.h"
#include <iostream>


//! \brief 1-point rule coordinates.
static const double T1[2] = { 0.3333333333333333, 0.3333333333333333 };
//! \brief 1-point rule weights.
static const double W1[1] = { 1.0 };

//! \brief 3-point rule coordinates.
static const double T3[6] = {
  0.5, 0.5,
  0.5, 0.5,
  0.0, 0.5
};
//! \brief 3-point rule weights.
static const double W3[3] = {
  0.3333333333333333,
  0.3333333333333333,
  0.3333333333333333
};

//! \brief 4-point rule coordinates.
static const double T4[8] = {
  0.3333333333333333, 0.3333333333333333,
  0.6, 0.2,
  0.2, 0.6,
  0.2, 0.2
};
//! \brief 4-point rule weights.
static const double W4[4] = {
 -0.5625,
  0.5208333333333333,
  0.5208333333333333,
  0.5208333333333333,
};

//! \brief Prints an error message for non-supported quadrature rules.
const double* error (int n)
{
  std::cerr <<" *** TriangleQuadrature: "<< n <<"-point rule is not available."
            << std::endl;
  return nullptr;
}


const double* TriangleQuadrature::getCoord (int n)
{
  switch (n) {
  case 0: return T1; // dummy
  case 1: return T1;
  case 3: return T3;
  case 4: return T4;
  }

  return error(n);
}


const double* TriangleQuadrature::getWeight (int n)
{
  switch (n) {
  case 0: return W1; // dummy
  case 1: return W1;
  case 3: return W3;
  case 4: return W4;
  }

  return error(n);
}
