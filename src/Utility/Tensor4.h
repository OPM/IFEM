// $Id$
//==============================================================================
//!
//! \file Tensor4.h
//!
//! \date Oct 27 2015
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of fourth-order tensors with some basic operations.
//!
//==============================================================================

#ifndef _TENSOR4_H
#define _TENSOR4_H

#include <vector>
#include <iostream>


/*!
  \brief Simple class for representing a non-symmetric fourth-order tensor.
*/

class Tensor4
{
protected:
  typedef unsigned short int t_ind; //!< Tensor index type (for convenience)

  t_ind             n; //!< Number of spatial dimensions for the tensor
  t_ind             m; //!< Dimension of the matrix representation
  std::vector<Real> v; //!< The actual tensor component values

  //! \brief Auxilliary method used by the constructors.
  virtual void redim(t_ind nsd) { n = nsd; m = n*n; }

  //! \brief Prints out the tensor to an output stream.
  virtual std::ostream& print(std::ostream& os) const;

private:
  //! \brief Returns a 0-based array index for the given tensor indices.
  inline t_ind index(t_ind i, t_ind j, t_ind k, t_ind l) const
  {
    return i-1 + n*(j-1 + n*(k-1 + n*(l-1)));
  }

public:
  //! \brief The default constructor creates a (scaled) identity tensor.
  Tensor4(t_ind nsd = 3, Real scale = Real(1), bool makeJ = false);
  //! \brief Constructor creating a tensor from a vector of components.
  //! \details The provided vector is assumed to contain the components of the
  //! matrix representation of the tensor, stored in a one-dimensional array.
  Tensor4(const std::vector<Real>& x, t_ind nsd = 3);

  //! \brief Sets \a this to the 0-tensor.
  void zero() { std::fill(v.begin(),v.end(),Real(0)); }

  //! \brief Type casting to a one-dimensional vector, for referencing.
  operator const std::vector<Real>&() const { return v; }
  //! \brief Type casting to a one-dimensional vector, for assignment.
  operator std::vector<Real>&() { return v; }

  //! \brief Reference through a pointer.
  const Real* ptr() const { return &v.front(); }

  //! \brief Index-1 based component reference.
  const Real& operator()(t_ind i, t_ind j, t_ind d, t_ind l) const;
  //! \brief Index-1 based component access.
  Real& operator()(t_ind i, t_ind j, t_ind k, t_ind l);

  //! \brief Assignment operator.
  Tensor4& operator=(const Tensor4& T);
  //! \brief Overloaded assignment operator.
  Tensor4& operator=(Real val);

  //! \brief Incrementation operator.
  Tensor4& operator+=(Real val);

  //! \brief Output stream operator.
  friend std::ostream& operator<<(std::ostream& os, const Tensor4& T)
  {
    return T.print(os);
  }
};


/*!
  \brief Simple class for representing a symmetric fourth-order tensor.
*/

class SymmTensor4 : public Tensor4
{
protected:
  //! \brief Auxilliary method used by the constructors.
  virtual void redim(t_ind nsd);

  //! \brief Prints out the tensor to an output stream.
  virtual std::ostream& print(std::ostream& os) const;

private:
  //! \brief Returns a 0-based array index for the given row and column indices.
  //! \details Symmetric second-order tensors in 3D are assumed stored with the
  //! following order in a one-dimensional array: s11, s22, s33, s12, s23, s13.
  inline t_ind index(t_ind i, t_ind j) const
  {
    if (i == j)
      return i-1; // diagonal term
    else if (n == 2)
      return 2;   // off-diagonal term (2D)
    else if (n == 1)
      return 0;

    if (i == j+1 || i+2 == j) std::swap(i,j);
    return i+2; // upper triangular term (3D)
  }

public:
  //! \brief The default constructor creates an identity tensor.
  SymmTensor4(t_ind nsd = 3, bool makeJ = false);
  //! \brief Constructor creating a tensor from a vector of components.
  //! \details The provided vector is assumed to contain the components of the
  //! matrix representation of the tensor, stored in a one-dimensional array.
  SymmTensor4(const std::vector<Real>& x, t_ind nsd = 3);

  //! \brief Index-1 based component reference.
  const Real& operator()(t_ind i, t_ind j, t_ind d, t_ind l) const;
  //! \brief Index-1 based component access.
  Real& operator()(t_ind i, t_ind j, t_ind k, t_ind l);
};

#endif
