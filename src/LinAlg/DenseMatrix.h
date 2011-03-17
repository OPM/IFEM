// $Id: DenseMatrix.h,v 1.16 2010-12-06 09:05:46 rho Exp $
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
  DenseMatrix(size_t m = 0, size_t n = 0) : myMat(m,n) { ipiv = 0; }
  //! \brief Copy constructor.
  DenseMatrix(const DenseMatrix& A);
  //! \brief Special constructor taking data from a one-dimensional array.
  DenseMatrix(const RealArray& data, size_t nrows = 0);
  //! \brief The destructor frees the dynamically allocated arrays.
  virtual ~DenseMatrix() { if (ipiv) delete[] ipiv; }

  //! \brief Returns the matrix type.
  virtual Type getType() const { return DENSE; }

  //! \brief Creates a copy of the system matrix and returns a pointer to it.
  virtual SystemMatrix* copy() const { return new DenseMatrix(*this); }

  //! \brief Resizes the matrix to dimension \f$r \times c\f$.
  //! \details Will preserve existing matrix content within the new dimension.
  bool redim(size_t r, size_t c);

  //! \brief Returns the dimension of the system matrix.
  //! \param[in] idim Which direction to return the dimension in.
  virtual size_t dim(int idim = 1) const;

  //! \brief Index-1 based element access.
  real& operator()(size_t r, size_t c) { return myMat(r,c); }
  //! \brief Index-1 based element reference.
  const real& operator()(size_t r, size_t c) const { return myMat(r,c); }

  //! \brief Initializes the element assembly process.
  //! \details Must be called once before the element assembly loop.
  //! \param[in] sam Auxilliary data describing the FE model topology, etc.
  virtual void initAssembly(const SAM& sam);

  //! \brief Initializes the matrix to zero assuming it is properly dimensioned.
  virtual void init();

  //! \brief Adds an element stiffness matrix into the system stiffness matrix.
  //! \param[in] eM  The element stiffness matrix
  //! \param[in] sam Auxilliary data describing the FE model topology,
  //!                nodal DOF status and constraint equations
  //! \param[in] e   Identifier for the element that \a eM belongs to
  //! \return \e true on successful assembly, otherwise \e false
  virtual bool assemble(const Matrix& eM, const SAM& sam, int e);
  //! \brief Adds an element stiffness matrix into the system stiffness matrix.
  //! \details When multi-point constraints are present, contributions from
  //! these are also added into the system right-hand-side load vector.
  //! \param[in] eM  The element stiffness matrix
  //! \param[in] sam Auxilliary data describing the FE model topology,
  //!                nodal DOF status and constraint equations
  //! \param     B   The system right-hand-side load vector
  //! \param[in] e   Identifier for the element that \a eM belongs to
  //! \return \e true on successful assembly, otherwise \e false
  virtual bool assemble(const Matrix& eM, const SAM& sam,
			SystemVector& B, int e);

  //! \brief Augments a similar matrix symmetrically to the current matrix.
  //! \param[in] B  The matrix to be augmented
  //! \param[in] r0 Row offset for the augmented matrix
  //! \param[in] c0 Column offset for the augmented matrix
  virtual bool augment(const SystemMatrix& B, size_t r0, size_t c0);

  //! \brief Adds a matrix with similar dimension to the current matrix.
  //! \param[in] B     The matrix to be added
  //! \param[in] alpha Scale factor for matrix \b B
  virtual bool add(const SystemMatrix& B, real alpha = 1.0);

  //! \brief Adds the diagonal matrix \f$\sigma\f$\b I to the current matrix.
  virtual bool add(real sigma);

  //! \brief Performs the matrix-vector multiplication \b C = \a *this * \b B.
  virtual bool multiply(const SystemVector& B, SystemVector& C);

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param B Right-hand-side vector on input, solution vector on output
  //! \param newLHS \e true if the left-hand-side matrix is updated
  virtual bool solve(SystemVector& B, bool newLHS = true);

  //! \brief Solves a standard symmetric-definite eigenproblem.
  //! \details The eigenproblem is assumed to be on the form
  //! \b A \b x = \f$\lambda\f$ \b x where \b A ( = \a *this )
  //! is assumed to be symmetric and positive definite.
  //! The eigenproblem is solved by the LAPack library subroutine \a DSYEVX.
  //! \sa LAPack library documentation.
  //! \param[out] eigVal Computed eigenvalues
  //! \param[out] eigVec Computed eigenvectors stored column by column
  //! \param[in] nev The number of eigenvalues and eigenvectors to compute
  bool solveEig(RealArray& eigVal, Matrix& eigVec, int nev);

  //! \brief Solves a non-symmetric eigenproblem.
  //! \details The eigenproblem is assumed to be on the form
  //! \b A \b x = \f$\lambda\f$ \b x where \b A ( = \a *this )
  //! is a square non-symmetric matrix.
  //! The eigenproblem is solved by the LAPack library subroutine \a DGEEV.
  //! \sa LAPack library documentation.
  //! \param[out] r_val Real part of the computed eigenvalues
  //! \param[out] c_val Complex part of the computed eigenvalues
  bool solveEigNon(RealArray& r_val, RealArray& c_val);

  //! \brief Solves a generalized symmetric-definite eigenproblem.
  //! \details The eigenproblem is assumed to be on the form
  //! \b A \b x = \f$\lambda\f$ \b B \b x where \b A ( = \a *this ) and \b B
  //! both are assumed to be symmetric and \b B also to be positive definite.
  //! The eigenproblem is solved by the LAPack library subroutine \a DSYGVX.
  //! \sa LAPack library documentation.
  //! \param B Symmetric and positive definite mass matrix.
  //! \param[out] eigVal Computed eigenvalues
  //! \param[out] eigVec Computed eigenvectors stored column by column
  //! \param[in] nev The number of eigenvalues and eigenvectors to compute
  //! \param[in] shift Eigenvalue shift (unused)
  bool solveEig(DenseMatrix& B, RealArray& eigVal, Matrix& eigVec, int nev,
		real shift = 0.0);

  //! \brief Returns the L-infinity norm of the matrix.
  virtual real Linfnorm() const { return myMat.normInf(); }

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

  //! \brief Writes the system matrix to the given output stream.
  virtual std::ostream& write(std::ostream& os) const { return os << myMat; }

private:
  Matrix myMat; //!< The actual dense matrix
  int*   ipiv;  //!< Pivot indices used in \a solve
};

#endif
