// $Id$
//==============================================================================
//!
//! \file DiagMatrix.h
//!
//! \date Aug 29 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Diagonal system matrix representation.
//!
//==============================================================================

#ifndef _DIAG_MATRIX_H
#define _DIAG_MATRIX_H

#include "SystemMatrix.h"


/*!
  \brief Class for representing a diagonal system matrix.
*/

class DiagMatrix : public SystemMatrix
{
public:
  //! \brief Default constructor.
  DiagMatrix(size_t m = 0) : myMat(m) {}
  //! \brief Copy constructor.
  DiagMatrix(const DiagMatrix& A) : SystemMatrix(A), myMat(A.myMat) {}
  //! \brief Special constructor taking data from a one-dimensional array.
  DiagMatrix(const RealArray& data, size_t nrows = 0);
  //! \brief Empty destructor.
  virtual ~DiagMatrix() {}

  //! \brief Returns the matrix type.
  virtual LinAlg::MatrixType getType() const { return LinAlg::DIAG; }

  //! \brief Creates a copy of the system matrix and returns a pointer to it.
  virtual SystemMatrix* copy() const { return new DiagMatrix(*this); }

  //! \brief Resizes the matrix to dimension \f$r \times r\f$.
  //! \details Will preserve existing matrix content within the new dimension.
  bool redim(size_t r) { return myMat.resize(r,utl::RETAIN); }

  //! \brief Returns the dimension of the system matrix.
  virtual size_t dim(int) const { return myMat.size(); }

  //! \brief Access to the matrix itself.
  Vector& getMat() { return myMat; }
  //! \brief Index-1 based element access.
  Real& operator()(size_t r) { return myMat(r); }
  //! \brief Index-1 based element reference.
  const Real& operator()(size_t r) const { return myMat(r); }

  //! \brief Dumps the system matrix on a specified format.
  virtual void dump(std::ostream& os, LinAlg::StorageFormat format,
                    const char* label);

  //! \brief Initializes the element assembly process.
  //! \details Must be called once before the element assembly loop.
  //! \param[in] sam Auxiliary data describing the FE model topology, etc.
  virtual void initAssembly(const SAM& sam, bool);

  //! \brief Initializes the matrix to zero assuming it is properly dimensioned.
  virtual void init() { myMat.fill(Real(0)); }

  //! \brief Adds an element matrix into the associated system matrix.
  //! \param[in] eM  The element matrix
  //! \param[in] sam Auxiliary data for FE assembly management
  //! \param[in] e   Identifier for the element that \a eM belongs to
  //! \return \e true on successful assembly, otherwise \e false
  virtual bool assemble(const Matrix& eM, const SAM& sam, int e);
  //! \brief Adds an element matrix into the associated system matrix.
  //! \param[in] eM  The element matrix
  //! \param[in] sam Auxiliary data for FE assembly management
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

  //! \brief Multiplication with a scalar.
  virtual void mult(Real alpha) { myMat *= alpha; }

  //! \brief Adds a matrix with similar dimension to the current matrix.
  //! \param[in] B     The matrix to be added
  //! \param[in] alpha Scale factor for matrix \b B
  virtual bool add(const SystemMatrix& B, Real alpha);

  //! \brief Adds the diagonal matrix &sigma;\b I to the current matrix.
  virtual bool add(Real sigma);

  //! \brief Performs the matrix-vector multiplication \b C = \a *this * \b B.
  virtual bool multiply(const SystemVector& B, SystemVector& C) const;

  using SystemMatrix::solve;
  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param B Right-hand-side vector on input, solution vector on output
  virtual bool solve(SystemVector& B, Real*);

  //! \brief Returns the L-infinity norm of the matrix.
  virtual Real Linfnorm() const { return myMat.normInf(); }

protected:
  //! \brief Writes the system matrix to the given output stream.
  virtual std::ostream& write(std::ostream& os) const { return os << myMat; }

private:
  Vector myMat; //!< The actual diagonal matrix
};

#endif
