// $Id$
//==============================================================================
//!
//! \file FiniteElement.h
//!
//! \date Mar 30 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Finite element quantities at an integration point.
//!
//==============================================================================

#ifndef _FINITE_ELEMENT_H
#define _FINITE_ELEMENT_H

#include "MatVec.h"
#include "Vec3.h"


/*!
  \brief Class representing a finite element.
*/

class FiniteElement
{
public:
  //! \brief Default constructor.
  FiniteElement(size_t nb = 0) : iel(0), iGP(0), N(nb), detJxW(1.0) {}

  int      iel;         //!< Element identifier
  size_t   iGP;         //!< Global integration point counter
  double   u;           //!< First parameter of current point
  double   v;           //!< Second parameter of current point
  double   w;           //!< Third parameter of current point
  double   xi;          //!< First local coordinate of current integration point
  double   eta;         //!< Second local coordinate of current integration point
  double   zeta;        //!< Third local coordinate of current integration point
  Vector   N;           //!< Basis function values
  Vector   Navg;        //!< Volume-averaged basis function values
  Matrix   dNdX;        //!< First derivatives (gradient) of the basis functions
  Matrix3D d2NdX2;      //!< Second derivatives of the basis functions
  Matrix   G;           //!< Matrix used for stabilized methods
  std::vector<Vec3> XC; //!< Matrix with element corner coordinates
  double   detJxW;      //!< Weighted determinant of the coordinate mapping
};


/*!
  \brief Class representing a mixed finite element.
*/

class MxFiniteElement : public FiniteElement
{
public:
  //! \brief Default constructor.
  MxFiniteElement(size_t nb1 = 0, size_t nb2 = 0) : FiniteElement(nb1),
    N1(N), N2(nb2), dN1dX(dNdX) {}

  Vector&  N1;     //!< Basis function values for the first basis
  Vector   N2;     //!< Basis function values for the second basis
  Matrix&  dN1dX;  //!< Basis function gradients for the first basis
  Matrix   dN2dX;  //!< Basis function gradients for the second basis
};

#endif
