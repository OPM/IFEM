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

#include "ItgPoint.h"
#include "MatVec.h"
#include "Vec3.h"
#include "Tensor.h"


/*!
  \brief Class representing a finite element.
*/

class FiniteElement : public ItgPoint
{
public:
  //! \brief Default constructor.
  explicit FiniteElement(size_t n = 0, size_t i = 0) : ItgPoint(i), N(n), Te(3)
  { p = q = r = 0; h = 0.0; detJxW = 1.0; }

  //! \brief Empty destructor.
  virtual ~FiniteElement() {}

  //! \brief Returns the number of bases.
  virtual size_t getNoBasis() const { return 1; }

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

  //! \brief Returns a const reference to the basis function 3nd-derivatives.
  virtual const Matrix4D& hess2(char) const { return d3NdX3; }
  //! \brief Returns a reference to the basis function 3nd-derivatives.
  virtual Matrix4D& hess2(char) { return d3NdX3; }

protected:
  //! \brief Writes the finite element object to the given output stream.
  virtual std::ostream& write(std::ostream& os) const;

  //! \brief Global output stream operator.
  friend std::ostream& operator<<(std::ostream& os, const FiniteElement& fe)
  {
    return fe.write(os);
  }

public:
  // Gauss point quantities
  double   detJxW; //!< Weighted determinant of the coordinate mapping
  Vector     N;    //!< Basis function values
  Matrix    dNdX;  //!< First derivatives (gradient) of the basis functions
  Matrix3D d2NdX2; //!< Second derivatives of the basis functions
  Matrix4D d3NdX3; //!< Third derivatives of the basis functions
  Matrix     G;    //!< Covariant basis / Matrix used for stabilized methods
  Matrix     H;    //!< Hessian

  // Element quantities
  short int           p;    //!< Polynomial order of the basis in u-direction
  short int           q;    //!< Polynomial order of the basis in v-direction
  short int           r;    //!< Polynomial order of the basis in r-direction
  double              h;    //!< Characteristic element size/diameter
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
  //! \brief The constructor initializes the size of each basis.
  explicit MxFiniteElement(const std::vector<size_t>& n, size_t ip = 0);

  //! \brief Empty destructor.
  virtual ~MxFiniteElement() {}

  //! \brief Returns the number of bases.
  virtual size_t getNoBasis() const { return 1+M.size(); }

  //! \brief Returns a const reference to the basis function values.
  virtual const Vector& basis(char b) const { return b < 2 ? N : M[b-2]; }
  //! \brief Returns a reference to the basis function values.
  virtual Vector& basis(char b) { return b < 2 ? N : M[b-2]; }

  //! \brief Returns a const reference to the basis function derivatives.
  virtual const Matrix& grad(char b) const { return b < 2 ? dNdX : dMdX[b-2]; }
  //! \brief Returns a reference to the basis function derivatives.
  virtual Matrix& grad(char b) { return b < 2 ? dNdX : dMdX[b-2]; }

  //! \brief Returns a const reference to the basis function 2nd-derivatives.
  virtual const Matrix3D& hess(char b) const { return b < 2 ? d2NdX2 : d2MdX2[b-2]; }
  //! \brief Returns a reference to the basis function 2nd-derivatives.
  virtual Matrix3D& hess(char b) { return b < 2 ? d2NdX2 : d2MdX2[b-2]; }

  //! \brief Returns a const reference to the basis function 3rd-derivatives.
  virtual const Matrix4D& hess2(char b) const { return b < 2 ? d3NdX3 : d3MdX3[b-2]; }
  //! \brief Returns a reference to the basis function 3rd-derivatives.
  virtual Matrix4D& hess2(char b) { return b < 2 ? d3NdX3 : d3MdX3[b-2]; }

protected:
  //! \brief Writes the finite element object to the given output stream.
  virtual std::ostream& write(std::ostream& os) const;

private:
  std::vector<Vector>     M;    //!< Basis function values
  std::vector<Matrix>    dMdX;  //!< First derivatives of the basis functions
  std::vector<Matrix3D> d2MdX2; //!< Second derivatives of the basis functions
  std::vector<Matrix4D> d3MdX3; //!< Third derivatives of the basis functions
};

#endif
