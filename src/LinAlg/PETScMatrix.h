// $Id: PETScMatrix.h,v 1.6 2010-12-06 09:17:56 rho Exp $
//==============================================================================
//!
//! \file PETScMatrix.h
//!
//! \date Jan 15 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Representation of the system matrix in PETSc format with interface
//! to PETSc routines for assembling and solving linear equation systems.
//!
//==============================================================================

#ifndef _PETSC_MATRIX_H
#define _PETSC_MATRIX_H

#include "SystemMatrix.h"
#ifdef HAS_PETSC
#include "LinSolParams.h"
#include "petscksp.h"
#include "slepceps.h"
#endif


/*!
  \brief Class for representing the system vector in PETSc format.
  \details It is an interface to PETSc modules for assembling and solving
  linear systems of equations.
*/

class PETScVector : public SystemVector
{
public:
#ifdef HAS_PETSC
  //! \brief Constructor creating an empty vector.
  PETScVector();
  //! \brief Constructor creating a vector of length \a n.
  PETScVector(size_t n);
  //! \brief Constructor creating a vector from an array.
  PETScVector(const real* values, size_t n);
  //! \brief Copy constructor.
  PETScVector(const PETScVector& vec);
  //! \brief Destructor.
  virtual ~PETScVector();
#endif

  //! \brief Returns the vector type.
  virtual Type getType() const { return PETSC; }

#ifdef HAS_PETSC
  //! \brief Returns the dimension of the system vector.
  virtual size_t dim() const;

  //! \brief Sets the dimension of the system vector.
  virtual void redim(size_t n);

  //! \brief Creates a copy of the system vector and returns a pointer to it.
  virtual SystemVector* copy() const { return new PETScVector(*this); }

  //! \brief Access through pointer.
  virtual real* getPtr();
  //! \brief Reference through pointer.
  virtual const real* getRef() const;

  //! \brief Restores the vector contents from an array.
  virtual void restore(const real* ptr);

  //! \brief Initializes the vector to a given scalar value.
  virtual void init(real value = real(0));

  //! \brief Begins communication step needed in parallel vector assembly.
  //! \details Must be called together with endAssembly after vector assembly 
  //! is completed on each processor and before the linear system is solved.
  bool beginAssembly();

  //! \brief Ends communication step needed in parallel vector assembly.
  //! \details Must be called together with beginAssembly after vector assembly 
  //! is completed on each processor and before the linear system is solved.
  bool endAssembly();

  //! \brief Multiplication with a scalar.
  virtual void mult(real alpha);

  //! \brief L1-norm of vector.
  virtual real L1norm() const;

  //! \brief L2-norm of vector.
  virtual real L2norm() const;

  //! \brief Linfinity-norm of vector.
  virtual real Linfnorm() const;

  //! \brief Return PETSc vector.
  Vec& getVector() { return x; }

  //! \brief Return PETSc vector.
  const Vec& getVector() const { return x; }

private:
  Vec x; //!< The actual PETSc vector.

#else // dummy implementation when PETSc is not included
  PETScVector() {}
  virtual SystemVector* copy() const { return 0; }
  virtual size_t dim() const { return 0; }
  virtual void redim(size_t) {}
  virtual real* getPtr() { return 0; }
  virtual const real* getRef() const { return 0; }
  virtual void init(real = real(0)) {}
  virtual void mult(real) {}
  virtual real L1norm() const { return real(0); }
  virtual real L2norm() const { return real(0); }
  virtual real Linfnorm() const { return real(0); }
#endif
};


/*!
  \brief Class for representing the system matrix in PETSc format.
  \details It is an interface to PETSc modules for assembling and solving
  linear systems of equations.
*/

class PETScMatrix : public SystemMatrix
{
public:
#ifdef HAS_PETSC
  //! \brief Constructor.
  PETScMatrix(const LinSolParams& spar);
  //! \brief Copy constructor.
  PETScMatrix(const PETScMatrix& A);
  //! \brief The destructor frees the dynamically allocated arrays.
  virtual ~PETScMatrix();
#endif

  //! \brief Returns the matrix type.
  virtual Type getType() const { return PETSC; }

  //! \brief Returns the dimension of the system matrix.
  virtual size_t dim(int = 1) const { return 0; }

#ifdef HAS_PETSC
  //! \brief Creates a copy of the system matrix and returns a pointer to it.
  virtual SystemMatrix* copy() const { return new PETScMatrix(*this); }

  //! \brief Initializes the matrix to zero assuming it is properly dimensioned.
  virtual void init();

  //! \brief Begins communication step needed in parallel matrix assembly.
  //! \details Must be called together with endAssembly after matrix assembly 
  //! is completed on each processor and before the linear system is solved.
  bool beginAssembly();

  //! \brief Ends communication step needed in parallel matrix assembly.
  //! \details Must be called together with beginAssembly after matrix assembly 
  //! is completed on each processor and before the linear system is solved.
  bool endAssembly();

  //! \brief Initializes the element assembly process.
  //! \details Must be called once before the element assembly loop.
  //! The PETSc data structures are initialized and the all symbolic operations
  //! that are need before the actual assembly can start are performed.
  //! \param[in] sam Auxilliary data describing the FE model topology, etc.
  virtual void initAssembly(const SAM& sam);

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

  //! \brief Performs the matrix-vector multiplication \b C = \a *this * \b B.
  virtual bool multiply(const SystemVector& B, SystemVector& C);

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param B Right-hand-side vector on input, solution vector on output
  virtual bool solve(SystemVector& B, bool newLHSmatrix = true);

  //! \brief Solves a generalized symmetric-definite eigenproblem.
  //! \details The eigenproblem is assumed to be on the form
  //! \b A \b x = \f$\lambda\f$ \b B \b x where \b A ( = \a *this ) and \b B
  //! both are assumed to be symmetric and \b B also to be positive definite.
  //! The eigenproblem is solved by the SLEPc library subroutine \a EPSSolve.
  //! \sa SLEPc library documentation.
  //! \param B Symmetric and positive definite mass matrix.
  //! \param[out] eigVal Computed eigenvalues
  //! \param[out] eigVec Computed eigenvectors stored column by column
  //! \param[in] nev The number of eigenvalues and eigenvectors to compute
  //! \param[in] shift Eigenvalue shift
  //! \param[in] iop Option telling whether to factorize matrix \a A or \b B.
  bool solveEig(PETScMatrix& B, RealArray& eigVal, Matrix& eigVec, int nev,
		real shift = real(0), int iop = 1);

private:
  Mat                 A;         //!< Linear system matrix
  KSP                 ksp;       //!< Linear solver
  const LinSolParams& solParams; //!< Linear solver parameters
  //EPS                 eps;       //!< Eigenvalue solver

#else // dummy implementation when PETSc is not included
  PETScMatrix(const LinSolParams&) {}
  virtual SystemMatrix* copy() const { return 0; }
  virtual void init() {}
  virtual void initAssembly(const SAM&) {}
  virtual bool assemble(const Matrix&, const SAM&, int) { return false; }
  virtual bool assemble(const Matrix&, const SAM&,
			SystemVector&, int) { return false; }
#endif
};

#endif
