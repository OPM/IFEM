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


/*!
  \brief Class for storage of general algebraic system of equations.
*/

class AlgEqSystem : public GlobalIntegral
{
public:
  //! \brief The constructor sets its reference to SAM and ProcessAdm objects.
  AlgEqSystem(const SAM& s, const ProcessAdm* a = nullptr);

  //! \brief The destructor frees the dynamically allocated objects.
  virtual ~AlgEqSystem() { this->clear(); }

  //! \brief Allocates the system matrices of the specified format.
  //! \param[in] mtype The matrix format to use for all matrices
  //! \param[in] spar Input parameters for the linear equation solver
  //! \param[in] nmat Number of system matrices to allocate
  //! \param[in] nvec Number of system vectors to allocate
  //! \param[in] nscl Number of scalar quantities to allocate
  //! \param[in] withReactions If \e false, no reaction forces will be computed
  //! \param[in] num_threads_SLU Number of threads for SuperLU_MT
  bool init(LinAlg::MatrixType mtype, const LinSolParams* spar = nullptr,
            size_t nmat = 1, size_t nvec = 1, size_t nscl = 0,
            bool withReactions = false, int num_threads_SLU = 1);

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

  //! \brief Initializes the system matrices to zero.
  //! \param[in] initLHS If \e false, only initialize right-hand-side vectors
  virtual void initialize(bool initLHS);
  //! \brief Finalizes the system matrices after element assembly.
  //! \param[in] newLHS If \e false, only right-hand-side vectors was assembled
  virtual bool finalize(bool newLHS);

  //! \brief Adds a set of element matrices into the algebraic equation system.
  //! \param[in] elmObj Pointer to the element matrices to add into \a *this
  //! \param[in] elmId Global number of the element associated with \a *elmObj
  virtual bool assemble(const LocalIntegral* elmObj, int elmId);

  //! \brief Returns the number of right-hand-side vectors allocated.
  size_t getNoRHS() const { return b.size(); }

  //! \brief Returns the \a i'th matrix of the equation system.
  SystemMatrix* getMatrix(size_t i = 0) { return i < A.size() ? A[i]._A : 0; }
  //! \brief Returns the \a i'th right-hand-side vector of the equation system.
  SystemVector* getVector(size_t i = 0) { return i < b.size() ? b[i] : 0; }
  //! \brief Returns the \a i'th scalar quantity.
  double getScalar(size_t i = 0) { return i < c.size() ? c[i] : 0.0; }

  //! \brief Returns a pointer to the nodal reaction forces, if any.
  const Vector* getReactions() const { return R.empty() ? 0 : &R; }

private:
  //! \brief Struct defining a coefficient matrix and an associated RHS-vector.
  struct SysMatrixPair
  {
    SystemMatrix* _A; //!< The coefficient matrix
    SystemVector* _b; //!< Pointer to the associated right-hand-side vector

    //! \brief Constructor initializing the pointers to zero.
    SysMatrixPair() : _A(nullptr), _b(nullptr) {}
  };

  std::vector<SysMatrixPair> A; //!< The actual coefficient matrices
  std::vector<SystemVector*> b; //!< The actual right-hand-side vectors
  std::vector<double>        c; //!< Global scalar quantities
  std::vector<double>*       d; //!< Multithreading buffer for the scalar values
  Vector                     R; //!< Nodal reaction forces

  const SAM&        sam; //!< Data for FE assembly management
  const ProcessAdm* adm; //!< Parallel process administrator
};

#endif
