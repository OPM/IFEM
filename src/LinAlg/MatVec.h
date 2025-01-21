// $Id$
//==============================================================================
//!
//! \file MatVec.h
//!
//! \date Oct 1 2007
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Global algebraic operations on index 1-based matrices and vectors.
//!
//==============================================================================

#ifndef UTL_MATVEC_H
#define UTL_MATVEC_H

#include "matrixnd.h"

// Some convenience type definitions:

//! A real-valued vector with algebraic operations
typedef utl::vector<Real>   Vector;
//! A real-valued matrix with algebraic operations
typedef utl::matrix<Real>   Matrix;
//! A real-valued three-dimensional matrix with algebraic operations
typedef utl::matrix3d<Real> Matrix3D;
//! A real-valued four-dimensional matrix with algebraic operations
typedef utl::matrix4d<Real> Matrix4D;

//! A real-valued array without algebraic operations
typedef std::vector<Real>      RealArray;
//! A real-valued two-dimensional array without algebraic operations
typedef std::vector<RealArray> Real2DMat;
//! A real-valued three-dimensional array without algebraic operations
typedef std::vector<Real2DMat> Real3DMat;
//! An array of real-valued vectors with algebraic operations
typedef std::vector<Vector>    Vectors;
//! An array of real-valued matrices with algebraic operations
typedef std::vector<Matrix>    Matrices;


namespace utl
{
  //! \brief Multiplication of a vector and a scalar.
  //! \return \f$ {\bf Y} = c {\bf X} \f$
  Vector operator*(const Vector& X, Real c);
  //! \brief Multiplication of a scalar and a vector.
  //! \return \f$ {\bf Y} = c {\bf X} \f$
  inline Vector operator*(Real c, const Vector& X) { return X*c; }
  //! \brief Division of a vector by a scalar.
  //! \return \f$ {\bf Y} = \frac{1}{d} {\bf X} \f$
  inline Vector operator/(const Vector& X, Real d) { return X*(Real(1)/d); }

  //! \brief Addition of two vectors.
  //! \return \f$ {\bf Z} = {\bf X} + {\bf Y} \f$
  Vector operator+(const Vector& X, const Vector& Y);
  //! \brief Subtraction of two vectors.
  //! \return \f$ {\bf Z} = {\bf X} - {\bf Y} \f$
  Vector operator-(const Vector& X, const Vector& Y);

  //! \brief Multiplication of a matrix and a scalar.
  //! \return \f$ {\bf B} = c {\bf A} \f$
  Matrix operator*(const Matrix& A, Real c);
  //! \brief Multiplication of a scalar and a matrix.
  //! \return \f$ {\bf B} = c {\bf A} \f$
  inline Matrix operator*(Real c, const Matrix& A) { return A*c; }

  //! \brief Dot product of two vectors.
  //! \return \f$ a = {\bf X}^T {\bf Y} \f$
  inline Real operator*(const Vector& X, const Vector& Y) { return X.dot(Y); }

  //! \brief Multiplication of a matrix and a vector.
  //! \return \f$ {\bf Y} = {\bf A} {\bf X} \f$
  RealArray operator*(const Matrix& A, const Vector& X);
  //! \brief Multiplication of a vector and a matrix.
  //! \return \f$ {\bf Y} = {\bf A}^T {\bf X} \f$
  RealArray operator*(const Vector& X, const Matrix& A);

  //! \brief Multiplication of two matrices.
  //! \return \f$ {\bf C} = {\bf A} {\bf B} \f$
  Matrix operator*(const Matrix& A, const Matrix& B);

  //! \brief Congruence transformation of a symmetric matrix.
  //! \param A The matrix to be transformed
  //! \param[in] Tn Nodal transformation matrix
  bool transform(Matrix& A, const Matrix& Tn);

  //! \brief Congruence transformation of a vector.
  //! \param V The vector to be transformed
  //! \param[in] Tn Nodal transformation matrix
  //! \param[in] transpose If \e true, the transpose of \b Tn is used instead
  bool transform(Vector& V, const Matrix& Tn, bool transpose = false);

  //! \brief Inverts the square matrix \b A.
  bool invert(Matrix& A);
}

#endif
