// $Id$
//==============================================================================
//!
//! \file MatVec.h
//!
//! \date Oct 1 2007
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Index 1-based matrices and vectors for algebraic operations.
//!
//==============================================================================

#ifndef UTL_MATVEC_H
#define UTL_MATVEC_H

#include "matrix.h"

//! A real-valued vector with algebraic operations
typedef utl::vector<Real>   Vector;
//! A real-valued matrix with algebraic operations
typedef utl::matrix<Real>   Matrix;
//! A real-valued three-dimensional matrix with algebraic operations
typedef utl::matrix3d<Real> Matrix3D;

//! A real-valued array without algebraic operations
typedef std::vector<Real>      RealArray;
//! A real-valued two-dimensional array without algebraic operations
typedef std::vector<RealArray> Real2DMat;
//! A real-valued three-dimensional array without algebraic operations
typedef std::vector<Real2DMat> Real3DMat;
//! An array of real-valued vectors with algebraic operations
typedef std::vector<Vector>    Vectors;


namespace utl
{
  //! \brief Congruence transformation of a symmetric matrix.
  //! \details The following matrix multiplication is performed:
  //! \f[ {\bf A} = {\bf T}^T{\bf A}{\bf T} \f]
  //! where \b A is a full, symmetric matrix and \b T is an identity matrix
  //! with the nodal sub-matrix \b Tn inserted on the diagonal.
  //! \param A The matrix to be transformed
  //! \param[in] Tn Nodal transformation matrix
  //! \param[in] k Index telling where to insert \b Tn on the diagonal of \b T
  bool transform(Matrix& A, const Matrix& Tn, size_t k);

  //! \brief Congruence transformation of a vector.
  //! \details The vector \b V is pre-multiplied with the transformation matrix
  //! \b T which is the identity matrix with the nodal sub-matrix \b Tn
  //! inserted on the diagonal.
  //! \param V The vector to be transformed
  //! \param[in] Tn Nodal transformation matrix
  //! \param[in] k Index telling where to insert \b Tn on the diagonal of \b T
  //! \param[in] transpose If \e true, the transpose of \b Tn is used instead
  bool transform(Vector& V, const Matrix& Tn, size_t k, bool transpose = false);

  //! \brief Inverts the square matrix \b A.
  bool invert(Matrix& A);
}

#endif
