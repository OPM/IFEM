//==============================================================================
//!
//! \file PETScBlockMatrix.h
//!
//! \date Oct 20 2012
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Representation of the system matrix in PETSc block format with 
//! interface
//! to PETSc routines for assembling and solving linear equation systems.
//!
//==============================================================================

#ifndef _PETSC_BLOCK_MATRIX_H
#define _PETSC_BLOCK_MATRIX_H

#include "PETScMatrix.h"

typedef std::vector<int>                        IntVec;            //!< Vector of integers
typedef std::vector<std::vector<PetscIntVec> >  PetscIntVecVecVec; //!< Vector of vector of vector of PetscInt
typedef std::vector<IS>                         ISVec;             //!< Vector of PETSc index sets

/*!
  \brief Class for representing the system vector in PETSc format using
  block format.
  \details It is an interface to PETSc modules for assembling and solving
  linear systems of equations.
*/


class PETScBlockVector : public PETScVector
{
public:
#ifdef HAS_PETSC
  //! \brief Constructor creating an empty vector.
  PETScBlockVector();
  //! \brief Constructur setting number of blocks
  //! \param[in] n      Size of vector
  //! \param[in] ncomp  Number of components in each block
  PETScBlockVector(size_t n, IntVec ncomp);
  //! \brief Copy constructor.
  PETScBlockVector(const PETScBlockVector& vec);
  //! \brief Destructor.
  virtual ~PETScBlockVector();
#endif

  //! \brief Returns the vector type.
  virtual Type getType() const { return PETSCBLOCK; }

#ifdef HAS_PETSC
  //! \brief Returns the dimension of the system vector.
  virtual size_t dim() const;

//! \brief Returns the dimension of one block in the vector.
  virtual size_t dim(size_t n) const;

  //! \brief Sets the dimension of the system vector.
  virtual void redim(size_t n);

  //! \brief Sets the dimension of the system vector.
  virtual void redim(size_t n, IntVec ncomps);

  //! \brief Creates a copy of the system vector and returns a pointer to it.
  virtual SystemVector* copy() const { return new PETScBlockVector(*this); }

  //! \brief Access through pointer.
  virtual Real* getPtr();
  //! \brief Reference through pointer.
  virtual const Real* getRef() const;
  
  //! \brief Access to block through pointer.
  //! \param[in] i Block number 
  virtual Real* getPtr(size_t i);
  //! \brief Access to block through pointer.
  //! \param[in] i Block number 
  virtual const Real* getRef(size_t i) const;

  //! \brief Restores the vector contents from an array.
  virtual void restore(const Real* ptr);

   //! \brief Restores the vector contents from an array.
  //! \param[in] i Block number 
  virtual void restore(const Real* ptr, size_t i);

  //! \brief Initializes the vector to a given scalar value.
  virtual void init(Real value = Real(0));

   //! \brief Initializes the element assembly process.
  //! \details Must be called once before the element assembly loop.
  //! \param[in] sam Auxiliary data describing the FE model topology, etc.
  virtual void initAssembly(const SAM& sam);

  //! \brief Begins communication step needed in parallel vector assembly.
  //! \details Must be called together with endAssembly after vector assembly
  //! is completed on each processor and before the linear system is solved.
  virtual bool beginAssembly();

  //! \brief Ends communication step needed in parallel vector assembly.
  //! \details Must be called together with beginAssembly after vector assembly
  //! is completed on each processor and before the linear system is solved.
  virtual bool endAssembly();

  //! \brief Multiplication with a scalar.
  virtual void mult(Real alpha);

  //! \brief L1-norm of vector.
  virtual Real L1norm() const;

  //! \brief L2-norm of vector.
  virtual Real L2norm() const;

  //! \brief Linfinity-norm of vector.
  virtual Real Linfnorm() const;

  //! \brief Returns the PETSc vector (for assignment).
  virtual Vec& getVector() { return x; }
  //! \brief Returns the PETSc vector (for read access).
  virtual const Vec& getVector() const { return x; }

  //! \brief Returns the PETSc vector block (for assignment).
  virtual Vec& getVector(size_t i) { return bvecs[i]; }
  //! \brief Returns the PETSc vector block (for read access).
  virtual const Vec& getVector(size_t i) const
  { return bvecs[i]; }

 protected:
  size_t nblocks;       //!< Number of blocks
  IntVec ncomps;        //!< Number of components
  IS*    is;            //!< Global indices for each block
  Vec*   bvecs;         //!< Vectors for each block
#endif
};

/*!
  \brief Class for representing the system matrix in PETSc format.
  \details It is an interface to PETSc modules for assembling and solving
  linear systems of equations.
*/

class PETScBlockMatrix : public PETScMatrix
{
public:
#ifdef HAS_PETSC
  //! \brief Constructor.
  PETScBlockMatrix(const LinSolParams& spar);
  //! \brief Copy constructor.
  //PETScBlockMatrix(const PETScBlockMatrix& A);
  //! \brief Constructor defining the blocks
  //! \param[in] ncomp Number of components in each block
  //! \param[in] spar Linear solver parameters
  PETScBlockMatrix(IntVec ncomp, const LinSolParams& spar);
  //! \brief Destructor
  virtual ~PETScBlockMatrix();
#else
  //! \brief Constructor.
  PETScBlockMatrix(const LinSolParams&) {}
#endif

  //! \brief Returns the matrix type.
  virtual Type getType() const { return PETSCBLOCK; }

#ifdef HAS_PETSC
  //! \brief Returns the dimension of the system matrix.
  //virtual size_t dim(int i = 1) const;

  //! \brief Creates a copy of the system matrix and returns a pointer to it.
  virtual SystemMatrix* copy() const { return new PETScBlockMatrix(*this); }

  //! \brief Initialize blocks
  //! \param[in] nc Number of components for blocks
  virtual void setBlocks(IntVec nc) { ncomps = nc; }

  //! \brief Initializes the element assembly process.
  //! \details Must be called once before the element assembly loop.
  //! The PETSc data structures are initialized and the all symbolic operations
  //! that are needed before the actual assembly can start are performed.
  //! \param[in] sam Auxiliary data describing the FE model topology, etc.
  virtual void initAssembly(const SAM& sam, bool);

  //! \brief Initializes the matrix to zero assuming it is properly dimensioned.
  virtual void init();

  //! \brief Begins communication step needed in parallel matrix assembly.
  //! \details Must be called together with endAssembly after matrix assembly
  //! is completed on each processor and before the linear system is solved.
  virtual bool beginAssembly();
  //! \brief Ends communication step needed in parallel matrix assembly.
  //! \details Must be called together with beginAssembly after matrix assembly
  //! is completed on each processor and before the linear system is solved.
  virtual bool endAssembly();

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

  //! \brief Performs the matrix-vector multiplication \b C = \a *this * \b B.
  virtual bool multiply(const SystemVector& B, SystemVector& C);

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param B Right-hand-side vector on input, solution vector on output
  //! \param[in] newLHS \e true if the left-hand-side matrix has been updated
  virtual bool solve(SystemVector& B, bool newLHS = true);

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param[in] B Right-hand-side vector
  //! \param[out] x Solution vector
  //! \param[in] newLHS \e true if the left-hand-side matrix has been updated
  virtual bool solve(const SystemVector& B, SystemVector& x, bool newLHS = true);

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param b Right-hand-side vector, solution vector on output
  //! \param P Preconditioning matrix (if different than system matrix)
  //! \param[in] newLHS \e true if the left-hand-side matrix has been updated
  virtual bool solve(SystemVector& b, SystemMatrix& P, bool newLHS = true);

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
  virtual bool solveEig(PETScBlockMatrix& B, RealArray& eigVal, Matrix& eigVec, int nev,
			Real shift = Real(0), int iop = 1);

  //! \brief Returns the L-infinity norm of the matrix.
  virtual Real Linfnorm() const;

  //! \brief Returns the PETSc matrix (for assignment).
  virtual Mat& getMatrix() { return A; }
  //! \brief Returns the PETSc matrix (for read access).
  virtual const Mat& getMatrix() const { return A; }

  //! \brief Returns matrix block
  virtual Mat& getMatrixBlock(size_t i, size_t j) { return matvec[i*nblocks+j]; }
  //! \brief Returns matrix block (for read access).
  virtual const Mat& getMatrixBlock(size_t i, size_t j) const
  { return matvec[i*nblocks+j]; }

 protected:
  size_t nblocks;                     //!< Number of blocks
  IntVec ncomps;                      //!< Number of components
  IS* isvec;                          //!< Index set for blocks
  Mat* matvec;                        //!< Vector of matrix blocks 
  PetscIntVecVecVec locSubdDofsBlock; //!< Dofs for subdomains for block matrix
  PetscIntVecVecVec subdDofsBlock;    //!< Dofs for subdomains for block matrix


  void assemPETScBlock(const Matrix& eM, Mat SM, PETScVector& SV,
		       const std::vector<int>& meen, const int* meqn,
		       const int* mpmceq, const int* mmceq,
		       const Real* ttcc);
  
  //! \brief Function to decompose element matrix into different blocks
  //! \param[in] eM     Element system for coupled problem
  //! \param[in] l2g   Local to global mapping for element degrees of freedom
  //! \param[out] eMb  Element block matrices
  //! \param[out] l2gb Local to global mapping for block element degrees for freedom
  void getBlockElmMatData(const Matrix& eM, const PetscIntVec& l2g, 
			  std::vector<std::vector<Matrix> >& eMb, 
			  std::vector<PetscIntVec>& l2gb) const;

  void renumberRHS(const Vec& b, Vec& bnew, bool renum2block = true);

  void setParameters();
#endif
};

#endif
