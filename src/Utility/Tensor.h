// $Id$
//==============================================================================
//!
//! \file Tensor.h
//!
//! \date Dec 17 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of second-order tensors with some basic operations.
//!
//==============================================================================

#ifndef _TENSOR_H
#define _TENSOR_H

#include <vector>
#include <iostream>

class Vec3;


/*!
  \brief Simple class for representing a non-symmetric second-order tensor.
*/

class Tensor
{
protected:
  typedef unsigned short int t_ind; //!< Tensor index type

  const t_ind       n; //!< Number of spatial dimensions for the tensor
  std::vector<real> v; //!< The actual tensor component values

  //! \brief Return a 0-based array index for the given tensor indices.
  //! \details Assuming column-wise storage for non-symmetric tensors.
  virtual t_ind index(t_ind i, t_ind j) const { return i-1 + n*(j-1); }

  //! \brief Print out the tensor to an output stream.
  virtual std::ostream& print(std::ostream& os) const;

public:
  //! \brief Constructor creating a zero tensor.
  Tensor(const t_ind nsd) : n(nsd) { v.resize(n*n,real(0)); }
  //! \brief Constructor creating a transformation from two tangent vectors.
  Tensor(const std::vector<real>& t1, const std::vector<real>& t2);
  //! \brief Copy constructor.
  Tensor(const Tensor& T);

  //! \brief Set to 0-tensor.
  void zero() { std::fill(v.begin(),v.end(),real(0)); }

  //! \brief Type casting to a one-dimensional vector, for referencing.
  operator const std::vector<real>&() const { return v; }
  //! \brief Type casting to a one-dimensional vector, for assignment.
  operator std::vector<real>&() { return v; }

  //! \brief Reference through a pointer.
  const real* ptr() const { return &v.front(); }

  //! \brief Index-1 based component reference.
  const real& operator()(t_ind i, t_ind j) const { return v[this->index(i,j)]; }
  //! \brief Index-1 based component access.
  real& operator()(t_ind i, t_ind j) { return v[this->index(i,j)]; }

  //! \brief Assignment operator.
  Tensor& operator=(const Tensor& T);
  //! \brief Overloaded assignment operator.
  Tensor& operator=(const std::vector<real>& val);
  //! \brief Overloaded assignment operator.
  Tensor& operator=(real val);

  //! \brief Incrementation operator.
  Tensor& operator+=(const Tensor& T);
  //! \brief Incrementation operator.
  Tensor& operator+=(real val);

  //! \brief Scaling operator.
  Tensor& operator*=(real val);

  //! \brief Inner product of two tensors.
  real innerProd(const Tensor& T);

  //! \brief Return the dimension of this tensor.
  t_ind dim() const { return n; }

  //! \brief Query whether this tensor is symmetric or not.
  bool symmetric() const { return v.size() == (size_t)(n*(n+1)/2); }

  //! brief Query whether this tensor is zero within the given tolerance.
  bool isZero(real tol = real(1.0e-6)) const;

  //! \brief Transpose the tensor.
  virtual Tensor& transpose();

  //! \brief Compute the trace of the tensor.
  virtual real trace() const;

  //! \brief Compute the determinant of the tensor.
  virtual real det() const;

  //! \brief Invert the tensor.
  //! \param[in] tol Division by zero tolerance
  //! \return Determinant of the tensor
  virtual real inverse(real tol = real(0));

  // Global operators

  //! \brief Multiplication between a tensor and a point vector.
  friend Vec3 operator*(const Tensor& T, const Vec3& v);
  //! \brief Multiplication between a point vector and transpose of a tensor.
  friend Vec3 operator*(const Vec3& v, const Tensor& T);

  //! \brief Output stream operator.
  friend std::ostream& operator<<(std::ostream& os, const Tensor& T)
  {
    return T.print(os);
  }
};


/*!
  \brief Simple class for representing a symmetric second-order tensor.
*/

class SymmTensor : public Tensor
{
  //! \brief Sets the number of spatial dimensions for the tensor to \a newDim.
  //! \return \e true if the dimension was changed, otherwise \e false
  //!
  //! \details This method is private because the tensor dimension is not
  //! supposed to be changed by the application. It is only for internal use.
  bool redim(const t_ind newDim);

protected:
  //! \brief Return a 0-based array index for the given tensor indices.
  //! \details Symmetric tensors are assumed stored with the following order:
  //! s11, s22, s33, s12, s23, s13.
  virtual t_ind index(t_ind i, t_ind j) const
  {
    if (i == j)
      return i-1; // diagonal term
    else if (n == 2)
      return 2;   // off-diagonal term (2D)

    if (i == j+1 || i+2 == j) std::swap(i,j);
    return i+2; // upper triangular term (3D)
  }

  //! \brief Print out the lower triangle of the tensor to an output stream.
  virtual std::ostream& print(std::ostream& os) const;

public:
  //! \brief Constructor creating a zero tensor.
  SymmTensor(const t_ind nsd) : Tensor(nsd) { v.resize(n*(n+1)/2,real(0)); }
  //! \brief Constructor creating a symmetric tensor from a vector.
  SymmTensor(const std::vector<real>& vec);
  //! \brief Copy constructor.
  SymmTensor(const SymmTensor& T);

  //! \brief Transpose the symmetric tensor (do nothing).
  virtual Tensor& transpose() { return *this; }

  //! \brief Compute the trace of the symmetric tensor.
  virtual real trace() const;

  //! \brief Compute the determinant of the symmetric tensor.
  virtual real det() const;

  //! \brief Invert the symmetric tensor.
  //! \param[in] tol Division by zero tolerance
  //! \return Determinant of the tensor
  virtual real inverse(real tol = real(0));

  //! \brief Congruence transformation of a symmetric tensor.
  SymmTensor& transform(const Tensor& T);

  //! \brief Construct the right Cauchy-Green tensor from a deformation tensor.
  SymmTensor& rightCauchyGreen(const Tensor& F);

  //! \brief Return the von Mises value of the symmetric tensor.
  real vonMises() const;

  // Global operators

  //! \brief Multiplication between a scalar and a symmetric tensor.
  friend SymmTensor operator*(real a, const SymmTensor& T);

  //! \brief Adding a scaled unit tensor to a symmetric tensor.
  friend SymmTensor operator+(const SymmTensor& T, real a);
  //! \brief Subtracting a scaled unit tensor from a symmetric tensor.
  friend SymmTensor operator-(const SymmTensor& T, real a);
};


/*!
  \brief Simple class for representing a symmetric fourth-order tensor.
*/

class SymmTensor4
{
  typedef unsigned short int t_ind; //!< Tensor index type

  t_ind                    n; //!< Number of spatial dimensions for the tensor
  t_ind                    m; //!< Dimension of the matrix representation
  const std::vector<real>& v; //!< The actual tensor component values
  real*                  ptr; //!< Non-const pointer to tensor component values

  //! \brief Return a 0-based array index for the given row and column indices.
  //! \details Symmetric second-order tensors in 3D are assumed stored with the
  //! following order in a one-dimensional array: s11, s22, s33, s12, s23, s13.
  inline t_ind index(t_ind i, t_ind j) const
  {
    if (i == j)
      return i-1; // diagonal term
    else if (n == 2)
      return 2;   // off-diagonal term (2D)

    if (i == j+1 || i+2 == j) std::swap(i,j);
    return i+2; // upper triangular term (3D)
  }

public:
  //! \brief The constructor creates a tensor from a vector of components.
  //! \details The provided vector is assumed to contain the components of the
  //! matrix representation of the tensor, stored in a one-dimensional array.
  SymmTensor4(const std::vector<real>& x, t_ind nsd = 3);

  //! \brief Index-1 based component reference.
  const real& operator()(t_ind i, t_ind j, t_ind d, t_ind l) const;
  //! \brief Index-1 based component access.
  real& operator()(t_ind i, t_ind j, t_ind k, t_ind l);
};


/*!
  \brief Abstract interface to problem-specific local coordinate systems.
*/

class LocalSystem
{
protected:
  //! \brief Protected default constructor since this is an interface class.
  LocalSystem() {}

public:
  //! \brief Empty default destructor.
  virtual ~LocalSystem() {}

  //! \brief Computes the global-to-local transformation at the point \a X.
  virtual const Tensor& getTmat(const Vec3& X) const = 0;

  static int patch; //!< Counter used to establish multi-patch local systems
};

#endif
