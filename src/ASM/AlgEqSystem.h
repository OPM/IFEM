// $Id$
//==============================================================================
//!
//! \file AlgEqSystem.h
//!
//! \date Nov 11 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Storage of an algebraic equation system for a FEM problem.
//!
//==============================================================================

#ifndef _ALG_EQ_SYSTEM_H
#define _ALG_EQ_SYSTEM_H

#include "GlobalIntegral.h"
#include "SystemMatrix.h"

class SAM;
class LinSolParams;
class LocalIntegral;


/*!
  \brief Class for storage of general algebraic system of equations.
*/

class AlgEqSystem : public GlobalIntegral
{
public:
  //! \brief The constructor only sets its reference to the SAM object.
  AlgEqSystem(const SAM& _sam) : sam(_sam) {}

  //! \brief The destructor frees the dynamically allocated objects.
  virtual ~AlgEqSystem() { this->clear(); }

  //! \brief Allocates the system matrices of the specified format.
  //! \param[in] mtype The matrix format to use for all matrices
  //! \param[in] spar Input parameters for the linear equation solver
  //! \param[in] nmat Number of system matrices to allocate
  //! \param[in] nvec Number of system vectors to allocate
  //! \param[in] num_threads_SLU Number of threads for SuperLU_MT
  void init(SystemMatrix::Type mtype, const LinSolParams* spar,
	    size_t nmat, size_t nvec, int num_threads_SLU = 1);

  //! \brief Initializes the system matrices to zero.
  //! \param[in] initLHS If \e false, only initialize the right-hand-side vector
  void init(bool initLHS = true);

  //! \brief Erases the system matrices and frees dynamically allocated storage.
  void clear();

  //! \brief Associates a system vector to a system matrix.
  //! \param[in] imat Index of a coefficient matrix
  //! \param[in] ivec Index of the system vector to associate with the matrix
  //!
  //! \details The purpose of this method is to define which right-hand-side
  //! vector (if any) should receive contributions when assembling a
  //! coefficient matrix, when the system to be assembled has explicit
  //! constraint equations for which the dependent DOFs have been eliminated.
  bool setAssociatedVector(size_t imat, size_t ivec);

  //! Initializes the matrices to proper size for element assembly.
  void initAssembly();

  //! \brief Adds a set of element matrices into the algebraic equation system.
  //! \param[in] elmObj Pointer to the element matrices to add into \a *this
  //! \param[in] elmId Global number of the element associated with \a *elmObj
  virtual bool assemble(const LocalIntegral* elmObj, int elmId);

  //! \brief Returns the \a i'th matrix of the equation system.
  SystemMatrix* getMatrix(size_t i = 0) { return i < A.size() ? A[i]._A : 0; }

  //! \brief Returns the \a i'th right-hand-side vector of the equation system.
  SystemVector* getVector(size_t i = 0) { return i < b.size() ? b[i] : 0; }

  //! \brief Returns a pointer to the nodal reaction forces, if any.
  const Vector* getReactions() const { return R.empty() ? 0 : &R; }

private:
  //! \brief Struct defining a coefficient matrix and an associated RHS-vector.
  struct SysMatrixPair
  {
    SystemMatrix* _A; //!< The coefficient matrix
    SystemVector* _b; //!< Pointer to the associated right-hand-side vector
    //! \brief Constructor initializing the pointers to zero.
    SysMatrixPair() { _A = 0; _b = 0; }
  };

  std::vector<SysMatrixPair> A; //!< The actual coefficient matrices
  std::vector<SystemVector*> b; //!< The actual right-hand-side vectors
  Vector                     R; //!< Nodal reaction forces

  const SAM& sam; //!< Data for FE assembly management
};

#endif
