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
//! \cond DO_NOT_DOCUMENT
#ifdef USE_INT64
#include <cstdint>
#define Int_ std::int64_t
#else
#define Int_ int
#endif
#define NS 60
//! \endcond


/*!
  \brief Class for representing the system matrix on the SPR-format.
  \details It is an interface to a Fortran module for assembling and solving
  linear systems of equations.
*/

class SPRMatrix : public SystemMatrix
{
public:
  //! \brief Default constructor.
  SPRMatrix();
  //! \brief Copy constructor.
  SPRMatrix(const SPRMatrix& A);
  //! \brief The destructor frees the dynamically allocated arrays.
  virtual ~SPRMatrix();

  //! \brief Returns the matrix type.
  virtual LinAlg::MatrixType getType() const { return LinAlg::SPR; }

  //! \brief Creates a copy of the system matrix and returns a pointer to it.
  virtual SystemMatrix* copy() const { return new SPRMatrix(*this); }

  //! \brief This class is not copyable.
  SPRMatrix& operator=(const SPRMatrix&) = delete;

  //! \brief Returns the dimension of the system matrix.
  virtual size_t dim(int idim) const;

  //! \brief Initializes the element assembly process.
  //! \param[in] sam Auxiliary data describing the FE model topology, etc.
  virtual void initAssembly(const SAM& sam, char);

  //! \brief Initializes the matrix to zero assuming it is properly dimensioned.
  virtual void init();

  using SystemMatrix::assemble;
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
                        SystemVector& B, const IntVec& meq);

  //! \brief Adds an element matrix into the associated system matrix.
  //! \param[in] eM  The element matrix
  //! \param[in] meq Matrix of element equation numbers (0 based)
  //! \return \e true on successful assembly, otherwise \e false
  //!
  //! \details To be used when there is no underlying SAM
  virtual bool assemble(const Matrix& eM, const IntVec& meq);

  //! \brief Multiplication with a scalar.
  virtual void mult(Real alpha);

  //! \brief Adds a matrix with similar sparsity pattern to the current matrix.
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
  virtual bool solve(SystemVector& B, Real*);

  //! \brief Solves a generalized symmetric-definite eigenproblem.
  //! \param B Symmetric and positive definite mass matrix.
  //! \param[out] val Computed eigenvalues
  //! \param[out] vec Computed eigenvectors stored column by column
  //! \param[in] nev The number of eigenvalues and eigenvectors to compute
  //! \param[in] shift Eigenvalue shift
  //! \param[in] iop Option telling whether to factorize matrix \a A or \b B.
  bool solveEig(SPRMatrix& B, RealArray& val, Matrix& vec, int nev,
		Real shift = 0.0, int iop = 1);

  //! \brief Returns the L-infinity norm of the matrix.
  virtual Real Linfnorm() const;

  //! \brief Converts to a dense matrix.
  bool convert(Matrix& fullMat) const;

  //! \brief Dumps the system matrix on a specified format.
  virtual void dump(std::ostream& os, LinAlg::StorageFormat fmt, const char*);
  //! \brief Loads the system matrix from specified file.
  virtual bool load(const char* fileName, bool binary);

protected:
  //! \brief Adds an element matrix into the associated system matrix.
  //! \param[in] eM  The element matrix
  //! \param[in] sam Auxiliary data for FE assembly management
  //! \param     B   Pointer to the system right-hand-side vector
  //! \param[in] e   Identifier for the element that \a eM belongs to
  //! \return \e true on successful assembly, otherwise \e false
  bool assemble(int e, const Matrix& eM, const SAM& sam, Real* B = nullptr);
  //! \brief Performs the matrix-vector multiplication \b c = \a *this * \b b.
  bool multiply(size_t n, const Real* b, Real* c);

  //! \brief Writes the system matrix to the given output stream.
  virtual std::ostream& write(std::ostream& os) const;

private:
  Int_ ierr;     //!< Internal error flag
  Int_ mpar[NS]; //!< Matrix of sparse PARameters
  Int_* msica;   //!< Matrix of Storage Information for CA
  Int_* msifa;   //!< Matrix of Storage Information for FA
  Int_* mtrees;  //!< Matrix of elimination assembly TREES
  Int_* mvarnc;  //!< Matrix of VARiable to Node Correspondence
  Real* values;  //!< The actual matrix VALUES

  std::vector<Int_>  iWork; //!< Integer work array
  std::vector<Int_>* jWork; //!< Integer work arrays for multi-threaded assembly
  std::vector<Real>  rWork; //!< Real work array

#ifdef USE_INT64
  //! \brief Struct with SAM arrays needed in the matrix assembly process.
  struct SAM64
  {
    Int_* mpar;   //!< Matrix of parameters
    Int_* mpmnpc; //!< Matrix of pointers to MNPCs in MMNPC
    Int_* mmnpc;  //!< Matrix of matrices of nodal point correspondances
    Int_* madof;  //!< Matrix of accumulated DOFs
    Int_* mpmceq; //!< Matrix of pointers to MCEQs in MMCEQ
    Int_* mmceq;  //!< Matrix of matrices of constraint equation definitions
    Int_* meqn;   //!< Matrix of equation numbers

    //! \brief The constructor copies the 32-bit arrays into 64-bit versions.
    SAM64(const SAM& sam);
    //! \brief The destructor frees the dynamically allocated arrays.
    ~SAM64();
  };

  bool sharedSam = false; //!< True if SAM object is shared with another matrix
#else
  using SAM64 = SAM; //!< Convenience alias when using 32-bit int version
#endif
  SAM64* mySam; //!< Possibly 64-bit integer version of SAM arrays
};

#endif
