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

#include "SystemMatrix.h"
#include "PETScSupport.h"
#include "LinAlgenums.h"

class LinSolParams;

typedef std::vector<PetscInt>    PetscIntVec;  //!< PETSc integer vector
typedef std::vector<PetscIntVec> PetscIntMat;  //!< PETSc integer matrix
typedef std::vector<PetscReal>   PetscRealVec; //!< PETSc real vector
typedef std::vector<IS>          ISVec;        //!< Index set vector
typedef std::vector<ISVec>       ISMat;        //!< Index set matrix


/*!
  \brief Class for representing the system vector in PETSc format.
  \details It is an interface to PETSc modules for assembling and solving
  linear systems of equations.
*/

class PETScVector : public SystemVector
{
public:
  //! \brief Constructor creating an empty vector.
  PETScVector(const ProcessAdm& padm);
  //! \brief Constructor creating a vector of length \a n.
  PETScVector(const ProcessAdm& padm, size_t n);
  //! \brief Constructor creating a vector from an array.
  PETScVector(const ProcessAdm& padm, const Real* values, size_t n);
  //! \brief Copy constructor.
  PETScVector(const PETScVector& vec);
  //! \brief Destructor.
  virtual ~PETScVector();

  //! \brief Returns the vector type.
  virtual Type getType() const { return PETSC; }

  //! \brief Returns the dimension of the system vector.
  virtual size_t dim() const;

  //! \brief Sets the dimension of the system vector.
  virtual void redim(size_t n);

  //! \brief Creates a copy of the system vector and returns a pointer to it.
  virtual SystemVector* copy() const { return new PETScVector(*this); }

  //! \brief Access through pointer.
  virtual Real* getPtr();
  //! \brief Reference through pointer.
  virtual const Real* getRef() const;

  //! \brief Restores the vector contents from an array.
  virtual void restore(const Real* ptr);

  //! \brief Initializes the vector to a given scalar value.
  virtual void init(Real value = Real(0));

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

  //! \brief Addition of another system vector.
  virtual void add(const SystemVector& vec, double s) {} // add this later...

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

  //! \brief Return associated process administrator
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

class PETScMatrix : public SystemMatrix
{
public:
  //! \brief Constructor.
  PETScMatrix(const ProcessAdm& padm, const LinSolParams& spar,
              LinAlg::LinearSystemType ltype);
  //! \brief Copy constructor.
  PETScMatrix(const PETScMatrix& A);
  //! \brief The destructor frees the dynamically allocated arrays.
  virtual ~PETScMatrix();

  //! \brief Returns the matrix type.
  virtual Type getType() const { return PETSC; }

  //! \brief Returns the dimension of the system matrix.
  virtual size_t dim(int = 1) const { return 0; }

  //! \brief Creates a copy of the system matrix and returns a pointer to it.
  virtual SystemMatrix* copy() const { return new PETScMatrix(*this); }

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
  virtual bool multiply(const SystemVector& B, SystemVector& C) const;

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param B Right-hand-side vector on input, solution vector on output
  //! \param[in] newLHS \e true if the left-hand-side matrix has been updated
  virtual bool solve(SystemVector& B, bool newLHS, Real*);

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param[in] B Right-hand-side vector
  //! \param[out] x Solution vector
  //! \param[in] newLHS \e true if the left-hand-side matrix has been updated
  virtual bool solve(const SystemVector& B, SystemVector& x, bool newLHS);

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
  virtual bool solveEig(PETScMatrix& B, RealArray& eigVal, Matrix& eigVec,
			int nev, Real shift = Real(0), int iop = 1);

  //! \brief Returns the L-infinity norm of the matrix.
  virtual Real Linfnorm() const;

  //! \brief Returns the PETSc matrix (for assignment).
  virtual Mat& getMatrix() { return A; }
  //! \brief Returns the PETSc matrix (for read access).
  virtual const Mat& getMatrix() const { return A; }

  //! \brief Set the linear solver parameters (solver type, preconditioner, tolerances).
  //! \param[in] P Preconditioner  matrix (ignored here)
  //! \param[in] Pb Preconditioner vector (ignored here)
  //! \return True on success
  virtual bool setParameters(PETScMatrix* P = nullptr, PETScVector* Pb = nullptr);
protected:
  //! \brief Constructs index set needed for element-by-element preconditioner.
  bool makeElementIS(const SAM& sam);

  //! \brief Constructs the EBE preconditioner of the given matrix.
  bool makeEBEpreconditioner(const Mat A, Mat* AeI);

  //! \brief Solve a linear system
  bool solve(const Vec& b, Vec& x, bool newLHS, bool knoll);

  Mat                 A;               //!< The actual PETSc matrix
  KSP                 ksp;             //!< Linear equation solver
  MatNullSpace*       nsp;             //!< Null-space of linear operator
  const ProcessAdm&   adm;             //!< Process administrator
  const LinSolParams& solParams;       //!< Linear solver parameters
  bool                setParams;       //!< If linear solver parameters are set
  IS*                 elmIS;           //!< Element index sets
  PetscInt            ISsize;          //!< Number of index sets/elements
  PetscIntMat         locSubdDofs;     //!< Degrees of freedom for unique subdomains
  PetscIntMat         subdDofs;        //!< Degrees of freedom for subdomains
  PetscRealVec        coords;          //!< Coordinates of local nodes (x0,y0,z0,x1,y1,...)
  ISMat               dirIndexSet;     //!< Direction ordering
  int                 nLinSolves;      //!< Number of linear solves
  LinAlg::LinearSystemType linsysType; //!< Linear system type
};


//! \brief Matrix-vector product
PETScVector operator*(const SystemMatrix& A, const PETScVector& b);

//! \brief Solve linear system
PETScVector operator/(const SystemMatrix& A, const PETScVector& b);

#endif
