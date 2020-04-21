// $Id$
//==============================================================================
//!
//! \file SIMCoupledSI.h
//!
//! \date Mar 19 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Template for a semi-implicit coupling of two nonlinear solvers.
//!
//==============================================================================

#ifndef SIM_COUPLED_SI_H_
#define SIM_COUPLED_SI_H_

#include "SIMCoupled.h"
#include "SIMenums.h"


/*!
  \brief Template class for semi-implicitly coupled simulators.
*/

template<class T1, class T2>
class SIMCoupledSI : public SIMCoupled<T1,T2>
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  SIMCoupledSI(T1& s1, T2& s2) : SIMCoupled<T1,T2>(s1,s2), maxIter(-1) {}
  //! \brief Empty destructor.
  virtual ~SIMCoupledSI() {}

  //! \brief Enable/disable the staggering iteration cycles.
  virtual void enableStaggering(bool enable = true)
  {
    maxIter = enable ? std::min(this->S1.getMaxit(),this->S2.getMaxit()) : 0;
  }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp, bool firstS1 = true)
  {
    if (maxIter < 0)
      maxIter = std::min(this->S1.getMaxit(),this->S2.getMaxit());

    if (tp.multiSteps())
      this->S1.getProcessAdm().cout <<"\n  step="<< tp.step
                                    <<"  time="<< tp.time.t << std::endl;

    SIM::ConvStatus conv = SIM::OK;
    for (tp.iter = 0; tp.iter <= maxIter && conv != SIM::CONVERGED; tp.iter++)
    {
      SIM::ConvStatus status1 = SIM::OK, status2 = SIM::OK;
      if (firstS1 && (status1 = this->S1.solveIteration(tp)) <= SIM::DIVERGED)
        return false;

      if ((status2 = this->S2.solveIteration(tp)) <= SIM::DIVERGED)
        return false;

      if (!firstS1 && (status1 = this->S1.solveIteration(tp)) <= SIM::DIVERGED)
        return false;

      if ((conv = this->checkConvergence(tp,status1,status2)) <= SIM::DIVERGED)
        return false;
    }

    this->S1.postSolve(tp);
    this->S2.postSolve(tp);
    tp.time.first = false;

    return true;
  }

  //! \brief Override this method to add additional convergence criteria.
  virtual SIM::ConvStatus checkConvergence(const TimeStep&,
                                           SIM::ConvStatus status1,
                                           SIM::ConvStatus status2)
  {
    if (status1 == status2)
      return status1;

    if (status1 == SIM::FAILURE || status2 == SIM::FAILURE)
      return SIM::FAILURE;

    if (status1 == SIM::DIVERGED || status2 == SIM::DIVERGED)
      return SIM::DIVERGED;

    return SIM::OK;
  }

protected:
  int maxIter; //!< Maximum number of iterations
};

#endif
