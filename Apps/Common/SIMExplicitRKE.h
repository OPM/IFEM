//==============================================================================
//!
//! \file SIMExplicitRKE.h
//!
//! \date Nov 28 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Explicit embedded Runge-Kutta based time stepping for SIM classes
//!
//==============================================================================

#ifndef SIM_EXPLICIT_RKE_H_
#define SIM_EXPLICIT_RKE_H_

#include "TimeIntUtils.h"

namespace TimeIntegration {

  //! \brief Explicit embedded Runge-Kutta based time stepping for SIM classes.
  //! \details Template can be instanced over any SIM implementing ISolver,
  //            and which derive from SIMbase.
  template<class Solver>
class SIMExplicitRKE : public SIMExplicitRK<Solver>
{
public:
  //! \brief Constructor
  //! \param solv The simulator to do time stepping for
  //! \param type The Runge-Kutta scheme to use
  //! \param tol Tolerance for ERK truncation error control
  SIMExplicitRKE(Solver& solv, Method type, double tol) :
    SIMExplicitRK<Solver>(solv, NONE), errTol(tol)
  {
    if (type == HEUNEULER) {
      this->RK.order = 1;
      this->RK.b.push_back(1.0);
      this->RK.b.push_back(0.0);
      bs.push_back(0.5);
      bs.push_back(0.5);
      this->RK.c.push_back(0.0);
      this->RK.A.redim(2,2);
      this->RK.A(2,1) = 1.0;
    }
    if (type == BOGACKISHAMPINE) {
      this->RK.order = 2;
      this->RK.b.push_back(7.0/24.0);
      this->RK.b.push_back(1.0/4.0);
      this->RK.b.push_back(1.0/3.0);
      this->RK.b.push_back(1.0/8.0);
      bs.push_back(2.0/9.0);
      bs.push_back(1.0/3.0);
      bs.push_back(4.0/9.0);
      bs.push_back(0.0);
      this->RK.c.push_back(0.0);
      this->RK.c.push_back(0.5);
      this->RK.c.push_back(0.75);
      this->RK.c.push_back(1.0);
      this->RK.A.redim(4,4);
      this->RK.A(2,1) = 0.5;
      this->RK.A(3,2) = 0.75;
      this->RK.A(4,1) = 2.0/9.0;
      this->RK.A(4,2) = 1.0/3.0;
      this->RK.A(4,3) = 2.0/9.0;
    }
    if (type == FEHLBERG) {
      this->RK.order = 4;
      this->RK.b.push_back(25.0/216.0);
      this->RK.b.push_back(0.0);
      this->RK.b.push_back(1408.0/2565.0);
      this->RK.b.push_back(2197.0/4104.0);
      this->RK.b.push_back(-1.0/5.0);
      this->RK.b.push_back(0.0);
      bs.push_back(16.0/135);
      bs.push_back(0.0);
      bs.push_back(6656.0/12825.0);
      bs.push_back(28561.0/56430.0);
      bs.push_back(-9.0/50.0);
      bs.push_back(2.0/55.0);
      this->RK.c.push_back(0.0);
      this->RK.c.push_back(1.0/4.0);
      this->RK.c.push_back(3.0/8.0);
      this->RK.c.push_back(12.0/13.0);
      this->RK.c.push_back(1.0);
      this->RK.c.push_back(1.0/2.0);
      this->RK.A.redim(6,6);
      this->RK.A(2,1) = 0.25;
      this->RK.A(3,1) = 3.0/32.0;
      this->RK.A(3,2) = 9.0/32.0;
      this->RK.A(4,1) =  1932.0/2197.0;
      this->RK.A(4,2) = -7200.0/2197.0;
      this->RK.A(4,3) =  7296.0/2197.0;
      this->RK.A(5,1) = 439.0/216.0;
      this->RK.A(5,2) = -8.0;
      this->RK.A(5,3) = 3680.0/513.0;
      this->RK.A(5,4) = -845.0/4104.0;
      this->RK.A(6,1) = -8.0/27;
      this->RK.A(6,2) = 2.0;
      this->RK.A(6,3) = -3544.0/2565.0;
      this->RK.A(6,4) = 1859.0/4104.0;
      this->RK.A(6,5) = -11.0/40.0;
    }
  }

  //! \copydoc ISolver::solveStep(TimeStep&)
  bool solveStep(TimeStep& tp)
  {
    this->solver.getProcessAdm().cout <<"\n  step = "<< tp.step <<"  time = "<< tp.time.t << std::endl;

    std::vector<Vector> stages;
    Vector prevSol = this->solver.getSolution();
    bool ok = this->solveRK(stages, tp);
    double prevEst = 1.0;
    while (ok && prevEst > errTol) {
      Vector error(prevSol);
      // construct the error estimate
      for (size_t i=0;i<stages.size();++i)
        error.add(stages[i], tp.time.dt*bs[i]);
      this->solver.applyDirichlet(error);
      error -= this->solver.getSolution();

      const size_t nf = this->solver.getNoFields(1);
      size_t iMax[nf];
      double dMax[nf];
      prevEst = this->solver.solutionNorms(error,dMax,iMax,nf);
      this->solver.getProcessAdm().cout << "Error estimate: " << prevEst << std::endl;
      if (prevEst > errTol) {
        if (!tp.cutback())
          return false;
        this->solver.getSolution() = prevSol;
        ok = this->solveRK(stages, tp);
      }
    }

    if (ok && prevEst > 0.0) {
      tp.time.dt = tp.time.dt * pow(errTol/prevEst,1.0/(this->RK.order+1));
      this->solver.getProcessAdm().cout << "adjusting step size to " << tp.time.dt << std::endl;
    }

    return ok;
  }
private:
  std::vector<double> bs; //!< Runge-Kutta coefficients for embedded method
  double errTol; //!< Truncation error tolerance
};

}

#endif
