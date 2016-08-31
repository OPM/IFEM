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

#pragma once

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
  FiniteElement(size_t n = 0, size_t i = 0) : iGP(i), N(n), iel(0), p(0), Te(3)
  { u = v = w = xi = eta = zeta = 0.0; detJxW = 1.0; }

  //! \brief Empty destructor.
  virtual ~FiniteElement() {}

  //! \brief Returns a const reference to the basis function values.
  virtual const Vector& basis(char) const { return N; }
  //! \brief Returns a reference to the basis function values.
  virtual Vector& basis(char) { return N; }

  //! \brief Returns a const reference to the basis function derivatives.
  virtual const Matrix& grad(char) const { return dNdX; }
  //! \brief Returns a reference to the basis function derivatives.
  virtual Matrix& grad(char) { return dNdX; }

  //! \brief Returns a const reference to the basis function 2nd-derivatives.
  virtual const Matrix3D& hess(char) const { return d2NdX2; }
  //! \brief Returns a reference to the basis function 2nd-derivatives.
  virtual Matrix3D& hess(char) { return d2NdX2; }

  // Gauss point quantities
  size_t   iGP;    //!< Global integration point counter
  double   u;      //!< First parameter of current point
  double   v;      //!< Second parameter of current point
  double   w;      //!< Third parameter of current point
  double   xi;     //!< First local coordinate of current integration point
  double   eta;    //!< Second local coordinate of current integration point
  double   zeta;   //!< Third local coordinate of current integration point
  double   detJxW; //!< Weighted determinant of the coordinate mapping
  Vector     N;    //!< Basis function values
  Matrix    dNdX;  //!< First derivatives (gradient) of the basis functions
  Matrix3D d2NdX2; //!< Second derivatives of the basis functions
  Matrix   G;      //!< Matrix used for stabilized methods

  // Element quantities
  int                 iel;  //!< Element identifier
  short int           p;    //!< Polynomial order of the basis functions
  Vec3Vec             XC;   //!< Array with element corner coordinate vectors
  Vector              Navg; //!< Volume-averaged basis function values
  Matrix              Xn;   //!< Matrix of element nodal coordinates
  Tensor              Te;   //!< Local-to-global element transformation matrix
  std::vector<Tensor> Tn;   //!< Array of element nodal rotation matrices
};

std::ostream& operator<<(std::ostream& os, const FiniteElement& obj);

/*!
  \brief Class representing a mixed finite element.
*/

class MxFiniteElement : public FiniteElement
{
public:
  //! \brief Default constructor.
  MxFiniteElement(const std::vector<size_t>& n, size_t i = 0) :
    FiniteElement(n.front(),i), Nx(n.size()-1), dNxdX(n.size()-1), d2NxdX2(n.size()-1)
  {}

  //! \brief Empty destructor.
  virtual ~MxFiniteElement() {}

  //! \brief Returns a const reference to the basis function values.
  virtual const Vector& basis(char b) const { return b == 1 ? N : Nx[b-2]; }
  //! \brief Returns a reference to the basis function values.
  Vector& basis(char b) { return b == 1 ? N : Nx[b-2]; }

  //! \brief Returns a const reference to the basis function derivatives.
  virtual const Matrix& grad(char b) const { return b == 1 ? dNdX : dNxdX[b-2]; }
  //! \brief Returns a reference to the basis function derivatives.
  virtual Matrix& grad(char b) { return b == 1 ? dNdX : dNxdX[b-2]; }

  //! \brief Returns a const reference to the basis function 2nd-derivatives.
  virtual const Matrix3D& hess(char b) const { return b == 1 ? d2NdX2 : d2NxdX2[b-2]; }
  //! \brief Returns a reference to the basis function 2nd-derivatives.
  virtual Matrix3D& hess(char b) { return b == 1 ? d2NdX2 : d2NxdX2[b-2]; }

protected:
  Vectors                 Nx;
  std::vector<Matrix>    dNxdX;
  std::vector<Matrix3D> d2NxdX2;
};

