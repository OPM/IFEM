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

#include "DomainDecomposition.h"
#include "SystemMatrix.h"
#include "PETScSupport.h"
#include "PETScSolParams.h"

#include <array>
#include <memory>

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


class PETScMatrix;

/*!
  \brief Class for representing a set of system vectors in PETSc format.
  \details Used for solving systems with multiple RHS vectors.
*/
class PETScVectors : public SystemVector
{
public:
  //! \brief Constructor creating a set of system vectors.
  //! \param A Matrix with vector layout to use
  //! \param nvec Number of RHS vectors
  PETScVectors(const PETScMatrix& A, int nvec);

  //! \brief Destructor frees up the dynamically allocated vectors.
  ~PETScVectors();

  //! \brief Returns the vector type.
  LinAlg::MatrixType getType() const override { return LinAlg::PETSC; }

  //! \brief Creates a copy of the system vector and returns a pointer to it.
  SystemVector* copy() const override { return nullptr; }

  //! \brief Returns the dimension of the system vector.
  size_t dim() const override { return myDim; }

  //! \brief Sets the dimension of the system vector.
  void redim(size_t n) override {}

  //! \brief Access through pointer.
  Real* getPtr() override { return nullptr; }
  //! \brief Reference through pointer.
  const Real* getRef() const override { return nullptr; }

  //! \brief Initializes the vector to a given scalar value.
  void init(Real value) override {}

  //! \brief Multiplication with a scalar.
  void mult(Real alpha) override {}

  //! \brief Addition of another system vector to this one.
  void add(const SystemVector& vec, Real scale) override {}

  //! \brief L1-norm of the vector.
  Real L1norm() const override { return 0.0; }

  //! \brief L2-norm of the vector.
  Real L2norm() const override { return 0.0; }

  //! \brief Linfinity-norm of the vector.
  Real Linfnorm() const override { return 0.0; }

  //! \brief Adds element vectors into the system vector.
  //! \param[in] vecs The element vectors
  //! \param[in] meqn Matrix of element equation numbers (0 based)
  void assemble(const Vectors& vecs, const IntVec& meqn, int = 0) override;

  //! \brief Returns a particular vector.
  //! \param idx Index of vector to return
  Vec get(size_t idx) { return vectors[idx]; }

  //! \brief Returns number of vectors.
  size_t size() const { return vectors.size(); }

private:
  size_t myDim; //!< Global dimension of vectors
  std::vector<Vec> vectors; //!< Array of vectors
  const PETScMatrix& myA; //!< Reference to matrix vectors are associated with
};


/*!
  \brief Class for representing the system matrix in PETSc format.
  \details It is an interface to PETSc modules for assembling and solving
  linear systems of equations.
*/

class PETScMatrix : public SystemMatrix
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

  //! \brief Returns the dimension of the system matrix.
  //! \param[in] idim Which direction to return the dimension in.
  size_t dim(int idim = 1) const override;

  //! \brief Initializes a scalar matrix based on element connections.
  //! \param[in] maxEq Maximum equation number
  //! \param[in] elms Element nodal connection
  //! \param[in] neighs Element graph
  //! \param[in] part Partitioning to use
  bool init(int maxEq,
            const IntMat* elms = nullptr,
            const IntMat* neighs = nullptr,
            const IntVec* part = nullptr);

  //! \brief Initializes the equation sparsity pattern based on element connections.
  //! \param[in] MMNPC Matrix of matrices of nodal point correspondances
  //! \param[in] nel Number of elements
  void preAssemble(const std::vector<IntVec>& MMNPC, size_t nel) override;

  //! \brief Initializes the element assembly process.
  //! \param[in] sam Auxiliary data describing the FE model topology, etc.
  //!
  //! \details Must be called once before the element assembly loop.
  //! The PETSc data structures are initialized and the all symbolic operations
  //! that are needed before the actual assembly can start are performed.
  void initAssembly(const SAM& sam, char) override;

  //! \brief Initializes the matrix to zero assuming it is properly dimensioned.
  void init() override;

  //! \brief Adds an element matrix into the associated system matrix.
  //! \param[in] eM  The element matrix
  //! \param[in] sam Auxiliary data for FE assembly management
  //! \param[in] e   Identifier for the element that \a eM belongs to
  //! \return \e true on successful assembly, otherwise \e false
  bool assemble(const Matrix& eM, const SAM& sam, int e) override;
  //! \brief Adds an element matrix into the associated system matrix.
  //! \param[in] eM  The element matrix
  //! \param[in] sam Auxiliary data for FE assembly management
  //! \param     B   The system right-hand-side vector
  //! \param[in] e   Identifier for the element that \a eM belongs to
  //! \return \e true on successful assembly, otherwise \e false
  //!
  //! \details When multi-point constraints are present, contributions from
  //! these are also added into the system right-hand-side vector, \a B.
  bool assemble(const Matrix& eM, const SAM& sam, SystemVector& B, int e) override;

  //! \brief Adds an element matrix into the associated system matrix.
  //! \param[in] eM  The element matrix
  //! \param[in] sam Auxiliary data for FE assembly management
  //! \param     B   The system right-hand-side vector
  //! \param[in] meq Matrix of element equation numbers
  //! \return \e true on successful assembly, otherwise \e false
  //!
  //! \details When multi-point constraints are present, contributions from
  //! these are also added into the system right-hand-side vector, \a B.
  bool assemble(const Matrix& eM, const SAM& sam,
                SystemVector& B, const IntVec& meq) override;

  //! \brief Adds an element matrix into the associated system matrix.
  //! \param[in] eM  The element matrix
  //! \param[in] meq Matrix of element equation numbers
  //! \return \e true on successful assembly, otherwise \e false
  bool assemble(const Matrix& eM, const IntVec& meq) override;

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

  //! \brief Multiplication with a scalar.
  void mult(Real alpha) override;

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

  //! \brief Solve for multiple right-hand-side vectors.
  //! \param B Vectors with right-hand-sides to solve for
  //! \param sField Resulting vectors stored as a matrix
  bool solveMultipleRhs(PETScVectors& B, Matrix& sField);

  //! \brief Returns a const-ref to domain decompositioning.
  const DomainDecomposition& getDD() const;

protected:
  //! \brief Solve a linear system.
  bool solve(const Vec& b, Vec& x, bool knoll);

  //! \brief Disabled copy constructor.
  PETScMatrix(const PETScMatrix& A) = delete;

  //! \brief Clone matrix data.
  PETScMatrix(const ProcessAdm& padm,
              const PETScSolParams& spar);

  //! \brief Setup sparsity pattern for model.
  //! \param[in] elms Elements on this process
  //! \param[in] sam Auxiliary data describing the FE model topology, etc.
  void setupSparsity(const IntVec& elms, const SAM& sam);

  //! \brief Setup sparsity pattern for block-matrices for a model.
  void setupBlockSparsity(const IntVec& elms, const SAM& sam);

  //! \brief Calculates blocks for global eqs.
  void setupGlb2Blk(const SAM& sam);

  //! \brief Sets up preallocator matrix.
  Mat preAllocator(const int nrows, const int ncols = 0) const;

  //! \brief Sets up preallocator matrices for blocks.
  std::vector<Mat> preAllocators() const;

  Mat                 pA;              //!< The actual PETSc matrix
  PetscInt            nrow;            //!< Number of matrix rows
  PetscInt            ncol;            //!< Number of matrix columns
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
  bool                factored;        //!< True if PETsc matrix has been factored

  IS glob2LocEq = nullptr; //!< Index set for global-to-local equations
  std::vector<Mat> matvec; //!< Blocks for block matrices

  std::vector<IS> isvec; //!< Index sets for blocks

  std::vector<std::array<int,2>> glb2Blk; //!< Maps equations to block and block eq
  std::unique_ptr<DomainDecomposition> m_dd{}; //!< Internal partitioning information
};


//! \brief Matrix-vector product.
PETScVector operator*(const SystemMatrix& A, const PETScVector& b);

//! \brief Solve linear system.
PETScVector operator/(SystemMatrix& A, const PETScVector& b);

#endif
