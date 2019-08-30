// $Id$
//==============================================================================
//!
//! \file SPRMatrix.h
//!
//! \date Jan 4 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of the system matrix on the SPR-format with interface
//! to a Fortran module for assembling and solving linear equation systems.
//!
//==============================================================================

#ifndef _SPR_MATRIX_H
#define _SPR_MATRIX_H

#include "SystemMatrix.h"

//! \brief Size of the MSPAR array.
#define NS 60


/*!
  \brief Class for representing the system matrix on the SPR-format.
  \details It is an interface to a Fortran module for assembling and solving
  linear systems of equations.
*/

class SPRMatrix : public SystemMatrix
{
public:
  //! \brief Default constructor.
  SPRMatrix() {}
  //! \brief Copy constructor.
  SPRMatrix(const SPRMatrix& A);
  //! \brief The destructor frees the dynamically allocated arrays.
  virtual ~SPRMatrix();

  //! \brief Returns the matrix type.
  virtual LinAlg::MatrixType getType() const { return LinAlg::SPR; }

  //! \brief Creates a copy of the system matrix and returns a pointer to it.
  virtual SystemMatrix* copy() const { return new SPRMatrix(*this); }

  //! \brief Returns the dimension of the system matrix.
  virtual size_t dim(int = 1) const { return mpar[7]; }

  //! \brief Initializes the element assembly process.
  //! \details Must be called once before the element assembly loop.
  //! The SPR data structures are initialized and the all symbolic operations
  //! that are need before the actual assembly can start are performed.
  //! \param[in] sam Auxiliary data describing the FE model topology, etc.
  virtual void initAssembly(const SAM& sam, bool);

  //! \brief Initializes the matrix to zero assuming it is properly dimensioned.
  virtual void init();

  using SystemMatrix::assemble;
  //! \brief Adds an element matrix into the associated system matrix.
  //! \param[in] eM  The element matrix
  //! \param[in] sam Auxiliary data describing the FE model topology,
  //!                nodal DOF status and constraint equations
  //! \param[in] e   Identifier for the element that \a eM belongs to
  //! \return \e true on successful assembly, otherwise \e false
  virtual bool assemble(const Matrix& eM, const SAM& sam, int e);
  //! \brief Adds an element matrix into the associated system matrix.
  //! \details When multi-point constraints are present, contributions from
  //! these are also added into the system right-hand-side vector.
  //! \param[in] eM  The element matrix
  //! \param[in] sam Auxiliary data describing the FE model topology,
  //!                nodal DOF status and constraint equations
  //! \param     B   The system right-hand-side vector
  //! \param[in] e   Identifier for the element that \a eM belongs to
  //! \return \e true on successful assembly, otherwise \e false
  virtual bool assemble(const Matrix& eM, const SAM& sam,
			SystemVector& B, int e);

  //! \brief Adds a matrix with similar sparsity pattern to the current matrix.
  //! \param[in] B     The matrix to be added
  //! \param[in] alpha Scale factor for matrix \b B
  virtual bool add(const SystemMatrix& B, Real alpha = 1.0);

  //! \brief Adds the diagonal matrix \f$\sigma\f$\b I to the current matrix.
  virtual bool add(Real sigma);

  //! \brief Performs the matrix-vector multiplication \b C = \a *this * \b B.
  virtual bool multiply(const SystemVector& B, SystemVector& C) const;

  using SystemMatrix::solve;
  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param B Right-hand-side vector on input, solution vector on output
  virtual bool solve(SystemVector& B, bool, Real*);

  //! \brief Solves a generalized symmetric-definite eigenproblem.
  //! \details The eigenproblem is assumed to be on the form
  //! \b A \b x = \f$\lambda\f$ \b B \b x where \b A ( = \a *this ) and \b B
  //! both are assumed to be symmetric and \b B also to be positive definite.
  //! The eigenproblem is solved by the SAM library subroutine \a SPRLAN.
  //! \sa SAM library documentation.
  //! \param B Symmetric and positive definite mass matrix.
  //! \param[out] eigVal Computed eigenvalues
  //! \param[out] eigVec Computed eigenvectors stored column by column
  //! \param[in] nev The number of eigenvalues and eigenvectors to compute
  //! \param[in] shift Eigenvalue shift
  //! \param[in] iop Option telling whether to factorize matrix \a A or \b B.
  bool solveEig(SPRMatrix& B, RealArray& eigVal, Matrix& eigVec, int nev,
		Real shift = 0.0, int iop = 1);

  //! \brief Returns the L-infinity norm of the matrix.
  virtual Real Linfnorm() const;

private:
  int mpar[NS] = {};      //!< Matrix of sparse PARameters
  int* msica = nullptr;   //!< Matrix of Storage Information for CA
  int* msifa = nullptr;   //!< Matrix of Storage Information for FA
  int* mtrees = nullptr;  //!< Matrix of elimination assembly TREES
  int* mvarnc = nullptr;  //!< Matrix of VARiable to Node Correspondence
  Real* values = nullptr; //!< The actual matrix VALUES

  std::vector<int>  iWork; //!< Integer work array
  std::vector<Real> rWork; //!< Real work array
};

#endif
