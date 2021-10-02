// $Id$
//==============================================================================
//!
//! \file ElmMats.h
//!
//! \date Nov 12 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of the element matrices for a FEM problem.
//!
//==============================================================================

#ifndef _ELM_MATS_H
#define _ELM_MATS_H

#include "LocalIntegral.h"
#include "MatVec.h"


/*!
  \brief Class collecting the element matrices associated with a FEM problem.
  \details The class is derived from LocalIntegral such that it may be passed
  as argument to the Integrand::evalInt and Integrand::evalBou methods.

  The class has two virtual methods returning the Newton matrix and the
  associated right-hand-side vector, which can be reimplemented by sub-classes
  to implement different time-integration schemes, mixed formulations, etc.
  The default versions use the first element matrix as the Newton matrix, and
  the first element vector as the associated right-hand-side vector.
*/

class ElmMats : public LocalIntegral
{
public:
  //! \brief Default constructor.
  explicit ElmMats(bool lhs = true) : rhsOnly(false), withLHS(lhs) {}
  //! \brief Empty destructor.
  virtual ~ElmMats() {}

  //! \brief Defines the number of element matrices and vectors.
  //! \param[in] nA Number of element matrices
  //! \param[in] nB Number of element vectors
  //! \param[in] nC Number of scalar quantities
  void resize(size_t nA, size_t nB, size_t nC = 0);

  //! \brief Sets the dimension of the element matrices and vectors.
  //! \param[in] ndim Number of rows and columns in the matrices/vectors
  void redim(size_t ndim);

  //! \brief Checks if the element matrices are empty.
  virtual bool empty() const { return A.empty() && b.empty(); }

  //! \brief Updates the time step size.
  virtual void setStepSize(double, int) {}

  //! \brief Returns the element-level Newton matrix.
  virtual const Matrix& getNewtonMatrix() const;
  //! \brief Returns the element-level right-hand-side vector.
  virtual const Vector& getRHSVector() const;

  std::vector<Matrix> A; //!< The element coefficient matrices
  std::vector<Vector> b; //!< The element right-hand-side vectors
  std::vector<double> c; //!< The scalar quantities

  bool rhsOnly; //!< If \e true, only the right-hand-sides are assembled
  bool withLHS; //!< If \e true, left-hand-side element matrices are present
};

#endif
