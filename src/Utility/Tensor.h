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
public:
  typedef unsigned short int t_ind; //!< Tensor index type
protected:
  const t_ind       n; //!< Number of spatial dimensions for the tensor
  std::vector<Real> v; //!< The actual tensor component values

  //! \brief Returns a 0-based array index for the given tensor indices.
  //! \details Assuming column-wise storage for non-symmetric tensors.
  virtual t_ind index(t_ind i, t_ind j) const { return i-1 + n*(j-1); }

  //! \brief Prints out the tensor to an output stream.
  virtual std::ostream& print(std::ostream& os) const;

private:
  //! \brief Creates a 3D transformation from three unit vectors.
  void define3Dtransform(const Vec3& v1, const Vec3& v2, const Vec3& v3);

public:
  //! \brief Constructor creating a zero tensor.
  Tensor(const t_ind nsd) : n(nsd) { v.resize(n*n,Real(0)); }
  //! \brief Constructor creating a transformation from a face normal vector.
  Tensor(const Vec3& vn);
  //! \brief Constructor creating a transformation from two tangent vectors.
  Tensor(const Vec3& t1, const Vec3& t2, bool t1isZ = false);
  //! \brief Constructor creating a transformation from three unit vectors.
  Tensor(const Vec3& v1, const Vec3& v2, const Vec3& v3);
  //! \brief Copy constructor.
  Tensor(const Tensor& T);

  //! \brief Sets \a this to the 0-tensor.
  void zero() { std::fill(v.begin(),v.end(),Real(0)); }

  //! \brief Type casting to a one-dimensional vector, for referencing.
  operator const std::vector<Real>&() const { return v; }
  //! \brief Type casting to a one-dimensional vector, for assignment.
  operator std::vector<Real>&() { return v; }

  //! \brief Reference through a pointer.
  const Real* ptr() const { return &v.front(); }

  //! \brief Index-1 based component reference.
  const Real& operator()(t_ind i, t_ind j) const { return v[this->index(i,j)]; }
  //! \brief Index-1 based component access.
  Real& operator()(t_ind i, t_ind j) { return v[this->index(i,j)]; }

  //! \brief Assignment operator.
  Tensor& operator=(const Tensor& T);
  //! \brief Overloaded assignment operator.
  Tensor& operator=(const std::vector<Real>& val);
  //! \brief Overloaded assignment operator.
  Tensor& operator=(Real val);

  //! \brief Incrementation operator.
  Tensor& operator+=(const Tensor& T);
  //! \brief Incrementation operator.
  Tensor& operator+=(Real val);

  //! \brief Decrementation operator.
  Tensor& operator-=(const Tensor& T);
  //! \brief Decrementation operator.
  Tensor& operator-=(Real val);

  //! \brief Scaling operator.
  Tensor& operator*=(Real val);

  //! \brief Returns the inner-product of \a *this and the given tensor.
  Real innerProd(const Tensor& T) const;

  //! \brief Returns the dimension of this tensor.
  t_ind dim() const { return n; }

  //! \brief Returns the size of this tensor.
  size_t size() const { return v.size(); }

  //! \brief Query whether this tensor is symmetric or not.
  virtual bool symmetric() const { return false; }

  //! brief Query whether this tensor is zero within the given tolerance.
  bool isZero(Real tol = Real(1.0e-6)) const;

  //! \brief Transposes the tensor.
  virtual Tensor& transpose();
  //! \brief Makes the tensor symmetric.
  virtual Tensor& symmetrize();

  //! \brief Returns the trace of the tensor.
  virtual Real trace() const;

  //! \brief Returns the determinant of the tensor.
  virtual Real det() const;

  //! \brief Inverts the tensor.
  //! \param[in] tol Division by zero tolerance
  //! \return Determinant of the tensor
  virtual Real inverse(Real tol = Real(0));

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
  //! \brief Resets the number of spatial dimensions of the tensor.
  //! \param[in] nsd The new tensor dimension
  //! \param[in] with33 If \e true and \a nsd = 2, include the 33 term also
  //! \return \e true if the dimension was changed, otherwise \e false
  //!
  //! \details This method is private because the tensor dimension is not
  //! supposed to be changed by the application. It is only for internal use.
  bool redim(const t_ind nsd, bool with33 = false);

protected:
  //! \brief Returns a 0-based array index for the given tensor indices.
  //! \details Symmetric 3D tensors are assumed stored with the following order:
  //! s11, s22, s33, s12, s23, s13. Symmetric 2D tensors have the order
  //! s11, s22, s12 and if the 33 component is included: s11, s22, s33, s12.
  virtual t_ind index(t_ind i, t_ind j) const
  {
    if (i == j)
      return i-1; // diagonal term
    else if (n == 2)
      return v.size()-1; // off-diagonal term (2D)

    if (i == j+1 || i+2 == j) std::swap(i,j);
    return i+2; // upper triangular term (3D)
  }

  //! \brief Prints out the lower triangle of the tensor to an output stream.
  virtual std::ostream& print(std::ostream& os) const;

public:
  //! \brief Constructor creating a zero tensor.
  //! \param[in] nsd The tensor dimension (1, 2 or 3)
  //! \param[in] with33 If \e true and \a nsd = 2, include the 33 term also
  //!
  //! \details The combination \a nsd = 2 and \a with33 = \e true results in a
  //! 2D tensor with the 33-term as the forth component. This is typically used
  //! to represent the symmetric stress tensor in 2D plane strain models,
  //! where the \f$\sigma_{zz}\f$ component is nonzero.
  SymmTensor(const t_ind nsd, bool with33 = false);
  //! \brief Constructor creating a symmetric tensor from a vector.
  SymmTensor(const std::vector<Real>& vec);
  //! \brief Copy constructor.
  SymmTensor(const SymmTensor& T) : Tensor(0) { this->copy(T); }

  //! \brief Copies a symmetric tensor, possibly with dimension change.
  void copy(const SymmTensor& T);

  //! \brief Query whether this tensor is symmetric or not.
  virtual bool symmetric() const { return true; }

  //! \brief Transposes the symmetric tensor (i.e., does nothing).
  virtual Tensor& transpose() { return *this; }
  //! \brief Makes the symmetric tensor symmetric (i.e., does nothing).
  virtual Tensor& symmetrize() { return *this; }

  //! \brief Returns the trace of the symmetric tensor.
  virtual Real trace() const;

  //! \brief Returns the determinant of the symmetric tensor.
  virtual Real det() const;

  //! \brief Inverts the symmetric tensor.
  //! \param[in] tol Division by zero tolerance
  //! \return Determinant of the tensor
  virtual Real inverse(Real tol = Real(0));

  //! \brief Congruence transformation of a symmetric tensor.
  SymmTensor& transform(const Tensor& T);

  //! \brief Constructs the right Cauchy-Green tensor from a deformation tensor.
  SymmTensor& rightCauchyGreen(const Tensor& F);

  //! \brief Returns the inner-product (L2-norm) of the symmetric tensor.
  Real L2norm(bool doSqrt = true) const;
  //! \brief Returns the von Mises value of the symmetric tensor.
  Real vonMises(bool doSqrt = true) const;
  //! \brief Computes the principal values of the symmetric tensor.
  void principal(Vec3& p) const;

  // Global operators

  //! \brief Adding a scaled unit tensor to a symmetric tensor.
  friend SymmTensor operator+(const SymmTensor& T, Real a);
  //! \brief Subtracting a scaled unit tensor from a symmetric tensor.
  friend SymmTensor operator-(const SymmTensor& T, Real a);

  //! \brief Multiplication between a scalar and a symmetric tensor.
  friend SymmTensor operator*(Real a, const SymmTensor& T);
};


//! \brief Adding two symmetric tensors.
inline SymmTensor operator+(const SymmTensor& A, const SymmTensor& B)
{
  SymmTensor C(A); C += B; return C;
}

//! \brief Subtracting two symmetric tensors.
inline SymmTensor operator-(const SymmTensor& A, const SymmTensor& B)
{
  SymmTensor C(A); C -= B; return C;
}

//! \brief Inner-product of two symmetric tensors.
inline Real operator*(const SymmTensor& A, const SymmTensor& B)
{
  return A.innerProd(B);
}


/*!
  \brief Simple class for representing a symmetric fourth-order tensor.
*/

class SymmTensor4
{
  typedef unsigned short int t_ind; //!< Tensor index type

  t_ind                    n; //!< Number of spatial dimensions for the tensor
  t_ind                    m; //!< Dimension of the matrix representation
  const std::vector<Real>& v; //!< The actual tensor component values
  Real*                  ptr; //!< Non-const pointer to tensor component values

  //! \brief Returns a 0-based array index for the given row and column indices.
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
  SymmTensor4(const std::vector<Real>& x, t_ind nsd = 3);

  //! \brief Index-1 based component reference.
  const Real& operator()(t_ind i, t_ind j, t_ind d, t_ind l) const;
  //! \brief Index-1 based component access.
  Real& operator()(t_ind i, t_ind j, t_ind k, t_ind l);
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
