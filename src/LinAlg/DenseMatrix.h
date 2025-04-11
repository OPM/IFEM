// $Id$
//==============================================================================
//!
//! \file DenseMatrix.h
//!
//! \date Jan 4 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Dense system matrix representation.
//!
//==============================================================================

#ifndef _DENSE_MATRIX_H
#define _DENSE_MATRIX_H

#include "SystemMatrix.h"

class SparseMatrix;


/*!
  \brief Class for representing a dense system matrix.
*/

class DenseMatrix : public SystemMatrix
{
public:
  //! \brief Default constructor.
  DenseMatrix(size_t m = 0, size_t n = 0, bool s = false);
  //! \brief Copy constructor.
  DenseMatrix(const DenseMatrix& A);
  //! \brief Special constructor taking data from a one-dimensional array.
  explicit DenseMatrix(const RealArray& data, size_t nrows = 0);
  //! \brief Special constructor, type conversion from Matrix.
  explicit DenseMatrix(const Matrix& A, bool s = false);
  //! \brief The destructor frees the dynamically allocated arrays.
  virtual ~DenseMatrix() { delete[] ipiv; }

  //! \brief Returns the matrix type.
  virtual LinAlg::MatrixType getType() const { return LinAlg::DENSE; }

  //! \brief Creates a copy of the system matrix and returns a pointer to it.
  virtual SystemMatrix* copy() const { return new DenseMatrix(*this); }

  //! \brief Resizes the matrix to dimension \f$r \times c\f$.
  //! \details Will preserve existing matrix content within the new dimension.
  bool redim(size_t r, size_t c);

  //! \brief Marks the matrix as symmetric.
  //! \details If marked as symmetric, Cholesky factorization will be employed.
  void setSymmetric(bool s = true) { symm = s && myMat.rows() == myMat.cols(); }

  //! \brief Returns the dimension of the system matrix.
  //! \param[in] idim Which direction to return the dimension in
  virtual size_t dim(int idim) const
  {
    return idim > 0 && idim < 3 ? myMat.dim(idim) : myMat.size();
  }

  //! \brief Access to the matrix itself.
  Matrix& getMat() { return myMat; }
  //! \brief Index-1 based element access.
  Real& operator()(size_t r, size_t c) { return myMat(r,c); }
  //! \brief Index-1 based element reference.
  const Real& operator()(size_t r, size_t c) const { return myMat(r,c); }

  //! \brief Dumps the system matrix on a specified format.
  virtual void dump(std::ostream& os, LinAlg::StorageFormat format,
                    const char* label);

  //! \brief Initializes the element assembly process.
  //! \param[in] sam Auxiliary data describing the FE model topology, etc.
  virtual void initAssembly(const SAM& sam, bool);

  //! \brief Initializes the matrix to zero assuming it is properly dimensioned.
  virtual void init();

  //! \brief Adds an element matrix into the associated system matrix.
  //! \param[in] eM  The element matrix
  //! \param[in] sam Auxiliary data for FE assembly management
  //! \param[in] e   Identifier for the element that \a eM belongs to
  //! \return \e true on successful assembly, otherwise \e false
  virtual bool assemble(const Matrix& eM, const SAM& sam, int e);
  //! \brief Adds an element matrix into the associated system matrix.
  //! \param[in] eM  The element matrix
  //! \param[in] sam Auxiliary data describing the FE model topology,
  //!                nodal DOF status and constraint equations
  //! \param     B   The system right-hand-side vector
  //! \param[in] e   Identifier for the element that \a eM belongs to
  //! \return \e true on successful assembly, otherwise \e false
  //!
  //! \details When multi-point constraints are present, contributions from
  //! these are also added into the system right-hand-side vector, \a B.
  virtual bool assemble(const Matrix& eM, const SAM& sam,
                        SystemVector& B, int e);
  //! \brief Adds an element matrix into the associated system matrix.
  //! \param[in] eM  The element matrix
  //! \param[in] sam Auxiliary data for FE assembly management
  //! \param     B   The system right-hand-side vector
  //! \param[in] meq Matrix of element equation numbers
  //! \return \e true on successful assembly, otherwise \e false
  //!
  //! \details When multi-point constraints are present, contributions from
  //! these are also added into the system right-hand-side vector, \a B.
  virtual bool assemble(const Matrix& eM, const SAM& sam,
                        SystemVector& B, const std::vector<int>& meq);

  //! \brief Augments a similar matrix symmetrically to the current matrix.
  //! \param[in] B  The matrix to be augmented
  //! \param[in] r0 Row offset for the augmented matrix
  //! \param[in] c0 Column offset for the augmented matrix
  virtual bool augment(const SystemMatrix& B, size_t r0, size_t c0);

  //! \brief Multiplication with a scalar.
  virtual void mult(Real alpha) { myMat.multiply(alpha); }

  //! \brief Adds a matrix with similar dimension to the current matrix.
  //! \param[in] B     The matrix to be added
  //! \param[in] alpha Scale factor for matrix \b B
  virtual bool add(const SystemMatrix& B, Real alpha);

  //! \brief Adds the constant &sigma; to the diagonal of this matrix.
  virtual bool add(Real sigma, int ieq);

  //! \brief Performs the matrix-vector multiplication \b C = \a *this * \b B.
  virtual bool multiply(const SystemVector& B, SystemVector& C) const;

  using SystemMatrix::solve;
  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param B Right-hand-side vector on input, solution vector on output
  //! \param[out] rc Reciprocal condition number of the LHS-matrix (optional)
  virtual bool solve(SystemVector& B, Real* rc = nullptr);
  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param B Right-hand-side matrix on input, solution matrix on output
  bool solve(Matrix& B);

  //! \brief Solves a standard symmetric-definite eigenproblem.
  //! \param[out] eigVal Computed eigenvalues
  //! \param[out] eigVec Computed eigenvectors stored column by column
  //! \param[in] nev The number of eigenvalues and eigenvectors to compute
  //!
  //! \details The eigenproblem is assumed to be on the form
  //! \b A \b x = &lambda; \b x where \b A ( = \a *this )
  //! is assumed to be symmetric and positive definite.
  //! The eigenproblem is solved by the LAPack library subroutine
  //! <a href="https://www.netlib.org/lapack/explore-3.1.1-html/dsyevx.f.html">DSYEVX</a>.
  //! \sa LAPack library documentation https://www.netlib.org/lapack/.
  bool solveEig(RealArray& eigVal, Matrix& eigVec, int nev);

  //! \brief Solves a non-symmetric eigenproblem.
  //! \sa LAPack library documentation.
  //! \param[out] r_val Real part of the computed eigenvalues
  //! \param[out] c_val Complex part of the computed eigenvalues
  //!
  //! \details The eigenproblem is assumed to be on the form
  //! \b A \b x = &lambda; \b x where \b A ( = \a *this )
  //! is a square non-symmetric matrix.
  //! The eigenproblem is solved by the LAPack library subroutine
  //! <a href="https://www.netlib.org/lapack/explore-3.1.1-html/dgeev.f.html">DGEEV</a>.
  //! \sa LAPack library documentation https://www.netlib.org/lapack/.
  bool solveEigNon(RealArray& r_val, RealArray& c_val);

  //! \brief Solves a generalized symmetric-definite eigenproblem.
  //! \param B Symmetric and positive definite mass matrix.
  //! \param[out] eigVal Computed eigenvalues
  //! \param[out] eigVec Computed eigenvectors stored column by column
  //! \param[in] nev The number of eigenvalues and eigenvectors to compute
  //! \param[in] shift Eigenvalue shift (unused)
  //!
  //! \details The eigenproblem is assumed to be on the form
  //! \b A \b x = &lambda; \b B \b x where \b A ( = \a *this ) and \b B
  //! both are assumed to be symmetric and \b B also to be positive definite.
  //! The eigenproblem is solved by the LAPack library subroutine
  //! <a href="https://www.netlib.org/lapack/explore-3.1.1-html/dsygvx.f.html">DSYGVX</a>.
  //! \sa LAPack library documentation https://www.netlib.org/lapack/.
  bool solveEig(DenseMatrix& B, RealArray& eigVal, Matrix& eigVec, int nev,
                Real shift = Real(0));

  //! \brief Returns the L-infinity norm of the matrix.
  virtual Real Linfnorm() const { return myMat.normInf(); }

protected:
  //! \brief Augments a dense matrix symmetrically to the current matrix.
  //! \param[in] B  The matrix to be augmented
  //! \param[in] r0 Row offset for the augmented matrix
  //! \param[in] c0 Column offset for the augmented matrix
  bool augment(const Matrix& B, size_t r0, size_t c0);

  //! \brief Augments a sparse matrix symmetrically to the current matrix.
  //! \param[in] B  The matrix to be augmented
  //! \param[in] r0 Row offset for the augmented matrix
  //! \param[in] c0 Column offset for the augmented matrix
  bool augment(const SparseMatrix& B, size_t r0, size_t c0);

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param B Right-hand-side vectors on input, solution vectors on output
  //! \param[in] nrhs Number of right-hand-side vectors
  //! \param[out] rcond Reciprocal condition number of the LHS-matrix (optional)
  //!
  //! \brief This is the function which actually solves the equation system,
  //! using the LAPack library subroutines. The two public \a solve methods just
  //! forward to this method.
  bool solve(Real* B, size_t nrhs, Real* rcond = nullptr);

  //! \brief Writes the system matrix to the given output stream.
  virtual std::ostream& write(std::ostream& os) const { return os << myMat; }

private:
  Matrix myMat; //!< The actual dense matrix
  int*   ipiv;  //!< Pivot indices used in \a solve
  bool   symm;  //!< Flags whether the matrix is symmetric or not
};

//! \brief Multiply a matrix with a scalar.
//! \param[in] alpha Scalar value
//! \param[in] A The matrix to scale
DenseMatrix operator*(Real alpha, const DenseMatrix& A);

//! \brief Multiply a matrix with a scalar.
//! \param[in] A The matrix to scale
//! \param[in] alpha Scalar value
DenseMatrix operator*(const DenseMatrix& A, Real alpha);

#endif
