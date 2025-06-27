// $Id$
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

#include "SparseMatrix.h"
#include "PETScSupport.h"
#include "PETScSolParams.h"
#include <array>

using PetscIntVec = std::vector<PetscInt>;    //!< PETSc integer vector
using PetscIntMat = std::vector<PetscIntVec>; //!< PETSc integer matrix
using PetscRealVec = std::vector<PetscReal>;  //!< PETSc real vector
using ISVec = std::vector<IS>;                //!< Index set vector
using ISMat = std::vector<ISVec>;             //!< Index set matrix


/*!
  \brief Class for representing the system vector in PETSc format.
  \details It is an interface to PETSc modules for assembling and solving
  linear systems of equations.
*/

class PETScVector : public StdVector
{
public:
  //! \brief Constructor creating an empty vector.
  explicit PETScVector(const ProcessAdm& padm);
  //! \brief Constructor creating a vector of length \a n.
  PETScVector(const ProcessAdm& padm, size_t n);
  //! \brief Constructor creating a vector from an array.
  PETScVector(const ProcessAdm& padm, const Real* values, size_t n);
  //! \brief Copy constructor.
  PETScVector(const PETScVector& vec);
  //! \brief Destructor.
  ~PETScVector() override;

  //! \brief Returns the vector type.
  LinAlg::MatrixType getType() const override { return LinAlg::PETSC; }

  //! \brief Creates a copy of the system vector and returns a pointer to it.
  SystemVector* copy() const override { return new PETScVector(*this); }

  //! \brief Initializes the vector to a given scalar value.
  void init(Real value) override;

  //! \brief Sets the dimension of the system vector.
  void redim(size_t n) override;

  //! \brief Finalizes the system vector assembly.
  bool endAssembly() override;

  //! \brief L1-norm of vector.
  Real L1norm() const override;

  //! \brief L2-norm of vector.
  Real L2norm() const override;

  //! \brief Linfinity-norm of vector.
  Real Linfnorm() const override;

  //! \brief Returns the PETSc vector (for assignment).
  Vec& getVector() { return x; }
  //! \brief Returns the PETSc vector (for read access).
  const Vec& getVector() const { return x; }

  //! \brief Return associated process administrator.
  const ProcessAdm& getAdm() const { return adm; }

protected:
  Vec x;                  //!< The actual PETSc vector
  const ProcessAdm& adm;  //!< Process administrator
};


/*!
  \brief Class for representing the system matrix in PETSc format.
  \details It is an interface to PETSc modules for assembling and solving
  linear systems of equations.
*/

class PETScMatrix : public SparseMatrix
{
public:
  //! \brief Constructor.
  PETScMatrix(const ProcessAdm& padm, const LinSolParams& spar);
  //! \brief The destructor frees the dynamically allocated arrays.
  ~PETScMatrix() override;

  //! \brief Returns the matrix type.
  LinAlg::MatrixType getType() const override { return LinAlg::PETSC; }

  //! \brief Creates a copy of the system matrix and returns a pointer to it.
  SystemMatrix* copy() const override;

  //! \brief Initializes the element assembly process.
  //! \param[in] sam Auxiliary data describing the FE model topology, etc.
  //!
  //! \details Must be called once before the element assembly loop.
  //! The PETSc data structures are initialized and the all symbolic operations
  //! that are needed before the actual assembly can start are performed.
  void initAssembly(const SAM& sam, char) override;

  //! \brief Initializes the matrix to zero assuming it is properly dimensioned.
  void init() override;

  //! \brief Finalizes the system matrix assembly.
  bool endAssembly() override;

  //! \brief Performs the matrix-vector multiplication \b C = \a *this * \b B.
  bool multiply(const SystemVector& B, SystemVector& C) const override;

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param B Right-hand-side vector on input, solution vector on output
  bool solve(SystemVector& B, Real*) override;

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param[in] B Right-hand-side vector
  //! \param[out] x Solution vector
  bool solve(const SystemVector& B, SystemVector& x) override;

  //! \brief Solves a generalized symmetric-definite eigenproblem.
  //! \details The eigenproblem is assumed to be on the form
  //! \b A \b x = &lambda; \b B \b x where \b A ( = \a *this ) and \b B
  //! both are assumed to be symmetric and \b B also to be positive definite.
  //! The eigenproblem is solved by the SLEPc library subroutine \a EPSSolve.
  //! \sa SLEPc library documentation.
  //! \param B Symmetric and positive definite mass matrix.
  //! \param[out] eigVal Computed eigenvalues
  //! \param[out] eigVec Computed eigenvectors stored column by column
  //! \param[in] nev The number of eigenvalues and eigenvectors to compute
  //! \param[in] shift Eigenvalue shift
  //! \param[in] iop Option telling whether to factorize matrix \a A or \b B.
  bool solveEig(PETScMatrix& B, RealArray& eigVal, Matrix& eigVec,
                int nev, Real shift = Real(0), int iop = 1);

  //! \brief Returns the L-infinity norm of the matrix.
  Real Linfnorm() const override;

  //! \brief Returns the PETSc matrix (for assignment).
  Mat& getMatrix() { return pA; }
  //! \brief Returns the PETSc matrix (for read access).
  const Mat& getMatrix() const { return pA; }

  //! \brief Get vector of block matrices. Used for tests only.
  const std::vector<Mat>& getBlockMatrices() const { return matvec; }

  //! \brief Get vector of index sets.
  const std::vector<IS>& getIS() const { return isvec; }

  //! \brief Set the linear solver parameters (solver type, preconditioner, tolerances).
  //! \param[in] setup True to setup KSP/PC
  //! \return True on success
  bool setParameters(bool setup);

  //! \brief Returns a const-ref to process administrator.
  const ProcessAdm& getAdm() const { return adm; }

protected:
  //! \brief Solve a linear system.
  bool solve(const Vec& b, Vec& x, bool knoll);

  //! \brief Solve system stored in the elem map.
  //! \details Create matrix from elem table, and solve for (possibly) multiple
  //!          right-hand-side vectors in B.
  //! \param B Vector with right-hand-sides to solve for
  bool solveDirect(PETScVector& B);

  //! \brief Assemble matrix directly from sparse matrix values.
  //! \details Assumes no DD, ie, do not use in parallel
  bool assembleDirect();

  //! \brief Disabled copy constructor.
  PETScMatrix(const PETScMatrix& A) = delete;

  //! \brief Clone sparse matrix data.
  PETScMatrix(const ProcessAdm& padm,
              const PETScSolParams& spar,
              const SparseMatrix& A);

  //! \brief Setup sparsity pattern for a DD partitioned model.
  //! \param[in] sam Auxiliary data describing the FE model topology, etc.
  void setupSparsityDD(const SAM& sam);
  //! \brief Setup sparsity pattern for a graph partitioned model.
  //! \param[in] sam Auxiliary data describing the FE model topology, etc.
  void setupSparsityPartitioned(const SAM& sam);
  //! \brief Setup sparsity pattern for a serial model.
  //! \param[in] sam Auxiliary data describing the FE model topology, etc.
  void setupSparsitySerial(const SAM& sam);

  //! \brief Setup sparsity pattern for block-matrices for a DD partitioned model.
  //! \param[in] sam Auxiliary data describing the FE model topology, etc.
  void setupBlockSparsityDD(const SAM& sam);
  //! \brief Setup sparsity pattern for block-matrices for a graph partitioned model.
  void setupBlockSparsityPartitioned(const SAM& sam);
  //! \brief Setup sparsity pattern for block-matrices for a serial model.
  void setupBlockSparsitySerial(const SAM& sam);

  //! \brief Calculates the global-to-block mapping for equations.
  std::vector<std::array<int,2>> setupGlb2Blk (const SAM& sam);
  //! \brief Calculates the global-to-block mapping for equations for a graph partitioned model.
  void setupGlb2BlkPart (const SAM& sam);

  Mat                 pA;              //!< The actual PETSc matrix
  KSP                 ksp;             //!< Linear equation solver
  MatNullSpace*       nsp;             //!< Null-space of linear operator
  const ProcessAdm&   adm;             //!< Process administrator
  PETScSolParams      solParams;       //!< Linear solver parameters
  bool                setParams;       //!< If linear solver parameters are set
  std::string         forcedKSPType;   //!< Force a KSP type ignoring the parameters
  PetscInt            ISsize;          //!< Number of index sets/elements
  PetscRealVec        coords;          //!< Coordinates of local nodes (x0,y0,z0,x1,y1,...)
  ISMat               dirIndexSet;     //!< Direction ordering
  int                 nLinSolves;      //!< Number of linear solves
  bool                assembled;       //!< True if PETSc matrix has been assembled

  IS glob2LocEq = nullptr; //!< Index set for global-to-local equations.
  std::vector<Mat> matvec; //!< Blocks for block matrices.

  std::vector<IS> isvec; //!< Index sets for blocks.
  std::vector<std::array<int,3>> glb2Blk; //!< Maps matrix entries in CSC order to block matrix entries.
};


//! \brief Matrix-vector product.
PETScVector operator*(const SystemMatrix& A, const PETScVector& b);

//! \brief Solve linear system.
PETScVector operator/(SystemMatrix& A, const PETScVector& b);

#endif
