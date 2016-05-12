// $Id$
//==============================================================================
//!
//! \file ISTLMatrix.h
//!
//! \date Mar 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Representation of the system matrix in ISTL format with interface
//! to ISTL routines for solving linear equation systems.
//!
//==============================================================================

#ifndef _ISTL_MATRIX_H
#define _ISTL_MATRIX_H

#include "SystemMatrix.h"
#include "SparseMatrix.h"
#include "ISTLSolParams.h"
#include "LinAlgenums.h"
#include <memory>

class LinSolParams;
class SAMpatch;


/*!
  \brief Class for representing the system vector in ISTL format.
*/

class ISTLVector : public StdVector
{
public:
  //! \brief Constructor creating an empty vector.
  ISTLVector(const ProcessAdm& padm);
  //! \brief Constructor creating a vector of length \a n.
  ISTLVector(const ProcessAdm& padm, size_t n);
  //! \brief Constructor creating a vector from an array.
  ISTLVector(const ProcessAdm& padm, const Real* values, size_t n);
  //! \brief Copy constructor.
  ISTLVector(const ISTLVector& vec);
  //! \brief Destructor.
  virtual ~ISTLVector();
#endif

  //! \brief Returns the vector type.
  virtual Type getType() const { return ISTL; }

  //! \brief Initializes the vector to a given scalar value.
  virtual void init(Real value = Real(0));

  //! \brief Returns the dimension of the system vector.
  virtual size_t dim() const;

  //! \brief Sets the dimension of the system vector.
  virtual void redim(size_t n);

  //! \brief Begins communication step needed in parallel vector assembly.
  //! \details Must be called together with endAssembly after vector assembly
  //! is completed on each processor and before the linear system is solved.
  virtual bool beginAssembly();

  //! \brief Ends communication step needed in parallel vector assembly.
  //! \details Must be called together with beginAssembly after vector assembly
  //! is completed on each processor and before the linear system is solved.
  virtual bool endAssembly();

  //! \brief L1-norm of vector.
  virtual Real L1norm() const;

  //! \brief L2-norm of vector.
  virtual Real L2norm() const;

  //! \brief Linfinity-norm of vector.
  virtual Real Linfnorm() const;

  //! \brief Return associated process administrator
  const ProcessAdm& getAdm() const { return adm; }
  typedef Dune::BlockVector<Dune::FieldVector<double,1>> Vec;

  Vec& getVector() { return x; }
  const Vec& getVector() const { return x; }

protected:
  ISTL::Vec x; //!< ISTL vector
  const ProcessAdm& adm;  //!< Process administrator
};


/*!
  \brief Class for representing the system matrix in ISTL format.
*/

class ISTLMatrix : public SparseMatrix
{
public:
  //! \brief Constructor.
  ISTLMatrix(const ProcessAdm& padm, const LinSolParams& spar,
             LinAlg::LinearSystemType ltype);
  //! \brief Copy constructor.
  ISTLMatrix(const ISTLMatrix& A);
  //! \brief The destructor frees the dynamically allocated arrays.
  virtual ~ISTLMatrix();

  //! \brief Returns the matrix type.
  virtual Type getType() const { return ISTL; }

  //! \brief Returns the dimension of the system matrix.
  virtual size_t dim(int = 1) const { return 0; }

  //! \brief Creates a copy of the system matrix and returns a pointer to it.
  virtual SystemMatrix* copy() const { return new ISTLMatrix(*this); }

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

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param B Right-hand-side vector on input, solution vector on output
  //! \param[in] newLHS \e true if the left-hand-side matrix has been updated
  virtual bool solve(SystemVector& B, bool newLHS, Real*);

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param[in] B Right-hand-side vector
  //! \param[out] x Solution vector
  //! \param[in] newLHS \e true if the left-hand-side matrix has been updated
  virtual bool solve(const SystemVector& B, SystemVector& x, bool newLHS);

  //! \brief Returns the L-infinity norm of the matrix.
  virtual Real Linfnorm() const;

  //! \brief Returns the ISTL matrix (for assignment).
  virtual ISTL::Mat& getMatrix() { return A; }
  //! \brief Returns the ISTL matrix (for read access).
  virtual const ISTL::Mat& getMatrix() const { return A; }

protected:
  ISTL::Mat A; //!< The actual ISTL matrix
  std::unique_ptr<ISTL::Operator> op; //!< The matrix adapter
  std::unique_ptr<ISTL::InverseOperator> solver; //!< Solver to use
  std::unique_ptr<ISTL::Preconditioner> pre; //!< Preconditioner to use
  const ProcessAdm&   adm;             //!< Process administrator
  ISTLSolParams       solParams;       //!< Linear solver parameters
  bool                setParams;       //!< If linear solver parameters are set
  int                 nLinSolves;      //!< Number of linear solves
  LinAlg::LinearSystemType linsysType; //!< Linear system type
};
