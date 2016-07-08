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


/*!
  \brief Template class for semi-implicitly coupled simulators.
*/

template<class T1, class T2>
class SIMCoupledSI : public SIMCoupled<T1, T2>
{
public:
  //! \brief The constructor initializes the references to the two solvers.
  SIMCoupledSI(T1& t1_, T2& t2_) : SIMCoupled<T1,T2>(t1_,t2_) {}
  //! \brief Empty destructor.
  virtual ~SIMCoupledSI() {}

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp)
  {
    int maxit = std::min(this->S1.getMaxit(),this->S2.getMaxit());

    this->S1.updateDirichlet(tp.time.t,&this->S1.getSolution(0));
    this->S2.updateDirichlet(tp.time.t,&this->S2.getSolution(0));
    this->S1.getProcessAdm().cout <<"\n  step="<< tp.step
                                  <<"  time="<< tp.time.t << std::endl;

    SIM::ConvStatus status1 = SIM::OK, status2 = SIM::OK;
    for (tp.iter = 0; tp.iter < maxit; tp.iter++)
    {
      if ((status1 = this->S1.solveIteration(tp)) == SIM::DIVERGED)
        return false;

      if ((status2 = this->S2.solveIteration(tp)) == SIM::DIVERGED)
        return false;

      SIM::ConvStatus cstatus = this->checkConvergence(tp, status1, status2);
      if (cstatus == SIM::CONVERGED)
	break;
      else if (cstatus == SIM::DIVERGED)
        return false;

      this->S1.updateDirichlet();
      this->S2.updateDirichlet();
    }

    tp.time.first = false;
    this->S1.postSolve(tp);
    this->S2.postSolve(tp);

    return true;
  }
};

#endif
