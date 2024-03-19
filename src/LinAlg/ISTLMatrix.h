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

#include "SparseMatrix.h"
#include "ISTLSolParams.h"
#include <memory>

class LinSolParams;


/*!
  \brief Class for representing the system vector in ISTL format.
*/

class ISTLVector : public StdVector
{
public:
  //! \brief Constructor creating an empty vector.
  explicit ISTLVector(const ProcessAdm& padm);
  //! \brief Constructor creating a vector of length \a n.
  ISTLVector(const ProcessAdm& padm, size_t n);
  //! \brief Constructor creating a vector from an array.
  ISTLVector(const ProcessAdm& padm, const Real* values, size_t n);
  //! \brief Copy constructor.
  ISTLVector(const ISTLVector& vec);
  //! \brief Destructor.
  ~ISTLVector() override;

  //! \brief Returns the vector type.
  LinAlg::MatrixType getType() const override { return LinAlg::ISTL; }

  //! \brief Initializes the vector to a given scalar value.
  void init(Real value = Real(0)) override;

  //! \brief Returns the dimension of the system vector.
  size_t dim() const override { return x.size(); }

  //! \brief Sets the dimension of the system vector.
  void redim(size_t n) override;

  //! \brief Copies the assembled vector into \ref x.
  bool endAssembly() override;

  //! \brief L1-norm of vector.
  Real L1norm() const override;

  //! \brief L2-norm of vector.
  Real L2norm() const override;

  //! \brief Linfinity-norm of vector.
  Real Linfnorm() const override;

  //! \brief Return associated process administrator
  const ProcessAdm& getAdm() const { return adm; }

  //! \brief Returns the ISTL vector (for assignment).
  ISTL::Vec& getVector() { return x; }
  //! \brief Returns the ISTL vector (for read access).
  const ISTL::Vec& getVector() const { return x; }

protected:
  ISTL::Vec         x;   //!< The actual ISTL vector
  const ProcessAdm& adm; //!< Process administrator
};


/*!
  \brief Class for representing the system matrix in ISTL format.
*/

class ISTLMatrix : public SparseMatrix
{
public:
  //! \brief Constructor.
  ISTLMatrix(const ProcessAdm& padm, const LinSolParams& spar);
  //! \brief Copy constructor.
  ISTLMatrix(const ISTLMatrix& A);
  //! \brief The destructor frees the dynamically allocated arrays.
  ~ISTLMatrix() override;

  //! \brief Returns the matrix type.
  LinAlg::MatrixType getType() const override { return LinAlg::ISTL; }

  //! \brief Returns the dimension of the system matrix.
  size_t dim(int = 1) const override { return 0; }

  //! \brief Creates a copy of the system matrix and returns a pointer to it.
  SystemMatrix* copy() const override { return new ISTLMatrix(*this); }

  //! \brief Initializes the element assembly process.
  //! \details Must be called once before the element assembly loop.
  //! The PETSc data structures are initialized and the all symbolic operations
  //! that are needed before the actual assembly can start are performed.
  //! \param[in] sam Auxiliary data describing the FE model topology, etc.
  //! \param[in] delayLocking If \e true, do not lock the sparsity pattern yet
  void initAssembly(const SAM& sam, bool delayLocking) override;

  //! \brief Initializes the matrix to zero assuming it is properly dimensioned.
  void init() override;

  //! \brief Copies the assembled matrix into \ref iA.
  bool endAssembly() override;

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param B Right-hand-side vector on input, solution vector on output
  bool solve(SystemVector& B, Real*) override;

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param[in] B Right-hand-side vector
  //! \param[out] x Solution vector
  bool solve(const SystemVector& B, SystemVector& x) override;

  //! \brief Returns the L-infinity norm of the matrix.
  Real Linfnorm() const override;

  //! \brief Returns the ISTL matrix (for assignment).
  ISTL::Mat& getMatrix() { return iA; }
  //! \brief Returns the ISTL matrix (for read access).
  const ISTL::Mat& getMatrix() const { return iA; }

protected:
  //! \brief Handles solver resets.
  void handleSolverReset();

  std::unique_ptr<ISTL::Operator> op; //!< The matrix adapter
  std::unique_ptr<ISTL::InverseOperator> solver; //!< Solver to use
  std::unique_ptr<ISTL::Preconditioner> pre; //!< Preconditioner to use

  ISTL::Mat         iA;        //!< The actual ISTL matrix
  const ProcessAdm& adm;       //!< Process administrator
  ISTLSolParams     solParams; //!< Linear solver parameters
  int               nLinSolves = 0; //!< Number of linear solves
};

#endif
