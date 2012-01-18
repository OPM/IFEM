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
  ElmMats() { rhsOnly = withLHS = false; }
  //! \brief Empty destructor.
  virtual ~ElmMats() {}

  //! \brief Defines the number of element matrices and vectors.
  //! \param[in] nA Number of element matrices
  //! \param[in] nB Number of element vectors
  void resize(size_t nA, size_t nB) { A.resize(nA); b.resize(nB); }

  //! \brief Returns the element-level Newton matrix.
  virtual const Matrix& getNewtonMatrix() const
  {
    if (A.empty())
      std::cerr <<" *** ElMats::getNewtonMatrix: No element matrices!"
		<< std::endl;
#if SP_DEBUG > 2 || INT_DEBUG > 0
    else
      std::cout <<"\nElement coefficient matrix"<< A.front();
#endif
    return A.front();
  }

  //! \brief Returns the element-level right-hand-side vector
  //! associated with the Newton matrix.
  virtual const Vector& getRHSVector() const
  {
    if (b.empty())
      std::cerr <<" *** ElMats::getNewtonMatrix: No element vectors!"
		<< std::endl;
#if SP_DEBUG > 2 || INT_DEBUG > 0
    else
      std::cout <<"\nElement right-hand-side vector"<< b.front();
#endif
    return b.front();
  }

  //! \brief Virtual destruction method to clean up after numerical integration.
  virtual void destruct() { delete this; }

  std::vector<Matrix> A; //!< The element coefficient matrices
  std::vector<Vector> b; //!< The element right-hand-side vectors

  bool rhsOnly; //!< If \e true, only the right-hand-sides are assembled
  bool withLHS; //!< If \e true, left-hand-side element matrices are present
};

#endif
