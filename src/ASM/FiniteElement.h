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
#include "Tensor.h"


/*!
  \brief Class representing a finite element.
*/

class FiniteElement
{
public:
  //! \brief Default constructor.
  FiniteElement(size_t n = 0, size_t i = 0) : iel(0), iGP(i), N(n), p(0), Te(3)
  { u = v = w = xi = eta = zeta = 0.0; detJxW = 1.0; }

  //! \brief Empty destructor.
  virtual ~FiniteElement() {}

  int iel; //!< Element identifier

  // Gauss point quantities
  size_t   iGP;    //!< Global integration point counter
  double   u;      //!< First parameter of current point
  double   v;      //!< Second parameter of current point
  double   w;      //!< Third parameter of current point
  double   xi;     //!< First local coordinate of current integration point
  double   eta;    //!< Second local coordinate of current integration point
  double   zeta;   //!< Third local coordinate of current integration point
  double   detJxW; //!< Weighted determinant of the coordinate mapping
  Vector   N;      //!< Basis function values
  Matrix   dNdX;   //!< First derivatives (gradient) of the basis functions
  Matrix3D d2NdX2; //!< Second derivatives of the basis functions
  Matrix   G;      //!< Matrix used for stabilized methods

  // Element quantities
  short int           p;    //!< Polynomial order of the basis functions
  Vec3Vec             XC;   //!< Array with element corner coordinate vectors
  Vector              Navg; //!< Volume-averaged basis function values
  Matrix              Xn;   //!< Matrix of element nodal coordinates
  Tensor              Te;   //!< Local-to-global element transformation matrix
  std::vector<Tensor> Tn;   //!< Array of element nodal rotation matrices
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

  //! \brief Empty destructor.
  virtual ~MxFiniteElement() {}

  Vector&  N1;     //!< Basis function values for the first basis
  Vector   N2;     //!< Basis function values for the second basis
  Matrix&  dN1dX;  //!< Basis function gradients for the first basis
  Matrix   dN2dX;  //!< Basis function gradients for the second basis
};

#endif
