// $Id$
//==============================================================================
//!
//! \file SIMCoupledSemiImplicit.h
//!
//! \date Mar 19 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Template for a semi-implicit coupling of two solvers
//!
//==============================================================================

#ifndef SIM_COUPLED_SEMIIMPLICIT_H_
#define SIM_COUPLED_SEMIIMPLICIT_H_

#include "SIMCoupled.h"

/*!
  \brief Template class for semi-implicitly coupled simulators.
*/

template<class T1, class T2>
class SIMCoupledSemiImplicit : public SIMCoupled<T1, T2> {
public:
  //! \brief The constructor initializes the references to the two solvers.
  SIMCoupledSemiImplicit(T1& t1_, T2& t2_) :
    SIMCoupled<T1,T2>(t1_, t2_)
  {
  }

  //! \brief Empty destructor.
  virtual ~SIMCoupledSemiImplicit() {}

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp)
  {
    double rtol, atol, dtol;
    int its1, its2;

    this->S2.getTolerances(atol, rtol, dtol, its2);
    this->S1.getTolerances(atol, rtol, dtol, its1);

    double r1=rtol+1.0, r2=rtol+1.0, r1r=1.0, r2r=1.0;

    this->S1.updateDirichlet(tp.time.t,&this->S1.getSolution(0));
    this->S2.updateDirichlet(tp.time.t,&this->S2.getSolution(0));
    int old1 = this->S1.msgLevel;
    this->S1.msgLevel = 0;
    if (this->S1.getGlobalProcessID() == 0)
      std::cout <<"\n  step="<< tp.step <<"  time="<< tp.time.t << std::endl;
    int it=0;
    bool abs1=true, abs2 = true, rel1=true, rel2=true;
    while (true) {
      ++it;
      for (int i=0;i<its1;++i && r1 > atol) {
        if (!this->S1.solveLinearizedSystem(tp,r1))
          return false;
        this->S1.updateDirichlet();
        if (it == 1 && i == 0)
          r1r = r1>0.0?r1:1.0;
      }
      for (int i=0;i<its2;++i && r2 > atol) {
        if (!this->S2.solveLinearizedSystem(tp,r2))
          return false;
        this->S2.updateDirichlet();
        if (it == 1 && i == 0)
          r2r = r2>0.0?r2:1.0;
      }

      if (r1 < atol)
        abs1 = false;
      if (r2 < atol)
        abs2 = false;
      if (r1/r1r < rtol)
        rel1 = false;
      if (r2/r2r < rtol)
        rel2 = false;

      if (r1/r1r > dtol || r2/r2r > dtol)
        return false;

      if (this->S1.getGlobalProcessID() == 0) {
        std::streamsize oldPrec = std::cout.precision(3);
        std::ios::fmtflags oldFlags = std::cout.flags(std::ios::scientific);
        std::cout << "  iter=" << it << " " << this->S1.getName() << ": conv=" << r1/r1r << "  incn=" << r1 << " " << this->S2.getName() << ": conv=" << r2/r2r << "  incn=" << r2 << std::endl;
        std::cout.precision(oldPrec);
        std::cout.flags(oldFlags);
      }
      if ((!abs1 || !rel1) && (!abs2 || !rel2))
        break;
    }
    this->S1.msgLevel = old1;

    this->S1.printSolutionSummary(this->S1.getSolution(0),0,nullptr,6);
    this->S2.printSolutionSummary(this->S2.getSolution(0),0,nullptr,6);

    return true;
  }
};

#endif
