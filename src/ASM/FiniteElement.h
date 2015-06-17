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

  //! \brief Returns a const reference to a vector of basis function values
  virtual const Vector& basis(char basis) const { return N; }
  //! \brief Returns a reference to a vector of basis function values
  virtual Vector& basis(char basis) { return N; }

  //! \brief Returns a const reference to a matrix of basis function derivatives
  virtual const Matrix& grad(char basis) const { return dNdX; }
  //! \brief Returns a reference to a matrix of basis function derivatives
  virtual Matrix& grad(char basis) { return dNdX; }

  //! \brief Returns a const reference to a matrix3d of basis function derivatives
  virtual const Matrix3D& hess(char basis) const { return d2NdX2; }
  //! \brief Returns a reference to a matrix3d of basis function derivatives
  virtual Matrix3D& hess(char basis) { return d2NdX2; }

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
  MxFiniteElement(const std::vector<size_t>& n) :
    FiniteElement(n[0]), Nx(n.size()-1), dNxdX(n.size()-1), d2NxdX2(n.size()-1)
  {}

  //! \brief Empty destructor.
  virtual ~MxFiniteElement() {}

  //! \brief Returns a const reference to a vector of basis function values
  const Vector& basis(char basis) const { return basis==1?N:Nx[basis-2]; }

  //! \brief Returns a reference to a vector of basis function values
  Vector& basis(char basis) { return basis==1?N:Nx[basis-2]; }

  //! \brief Returns a const reference to a matrix of basis function derivatives
  const Matrix& grad(char basis) const { return basis==1?dNdX:dNxdX[basis-2]; }

  //! \brief Returns a reference to a matrix of basis function derivatives
  Matrix& grad(char basis) { return basis==1?dNdX:dNxdX[basis-2]; }

  //! \brief Returns a const reference to a matrix3d of basis function derivatives
  const Matrix3D& hess(char basis) const { return basis==1?d2NdX2:d2NxdX2[basis-2]; }

  //! \brief Returns a reference to a matrix3d of basis function derivatives
  Matrix3D& hess(char basis) { return basis==1?d2NdX2:d2NxdX2[basis-2]; }

protected:
  Vectors Nx;
  std::vector<Matrix> dNxdX;
  std::vector<Matrix3D> d2NxdX2;
};

#endif
