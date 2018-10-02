//==============================================================================
//!
//! \file SIMMxV.h
//!
//! \date Nov 28 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Use a simulator as a matrix-vector product for iterative solvers.
//!
//==============================================================================
#ifndef SIM_MXV_H_
#define SIM_MXV_H_

#ifdef HAS_PETSC

#include "PETScMatrix.h"


/*!
 \brief Template for applying the matrix assembled in a simulator to a vector
        for use in iterative solvers.
*/

template<class Sim>
class SIMMxV : public PETScMxV {
public:
  //! \brief Default constructor.
  //! \param sim The simulator to wrap
  SIMMxV(Sim& sim) : S1(sim)
  {
    S1.getMatrix()->setMxV(this,true);
  }

  //! \brief Evaluate the matrix-vector product y = A*y
  bool evalMxV(Vec& x, Vec& y) override
  {
    MatMult(S1.getMatrix()->getBlockMatrices()[0], x, y);
    return true;
  }

  //! \brief Set a matrix-free preconditioner.
  void setPC(PETScPC* pc)
  {
    S1.getMatrix()->setPC(pc);
  }

  //! \brief Solve at time level.
  bool solveStep()
  {
    TimeStep tp;
    return S1.solveStep(tp);
  }

  //! \brief The matrix used for building the preconditioner.
  //! \details Defaults to the system matrix
  bool evalPC(Mat& P) override
  {
    MatDuplicate(S1.getMatrix()->getMatrix(), MAT_COPY_VALUES, &P);
    return true;
  }

protected:
  Sim& S1; //!< Reference to wrapped simulator
};

#endif

#endif
