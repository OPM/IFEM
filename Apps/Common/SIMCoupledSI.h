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
#include "MatVec.h"
#include "TimeStep.h"
#include "Functions.h"
#include "Utilities.h"
#include "IFEM.h"
#include <tinyxml.h>


/*!
  \brief Template class for semi-implicitly coupled simulators.
*/

template<class T1, class T2>
class SIMCoupledSI : public SIMCoupled<T1,T2>
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  SIMCoupledSI(T1& s1, T2& s2) : SIMCoupled<T1,T2>(s1,s2)
  {
    maxIter = 50;
    maxItFunc = nullptr;
    noStg = aitken = false;
    omega = omega0 = 0.0;
  }

  //! \brief The destructor deletes the max iteration function.
  virtual ~SIMCoupledSI() { delete maxItFunc; }

  //! \brief Enable/disable the staggering iteration cycles.
  virtual void enableStaggering(bool enable = true) { noStg = !enable; }

  //! \brief Returns residual to use for aitken acceleration.
  virtual const Vector& getAitkenResidual() const
  {
    static Vector empty;
    return empty;
  }

  //! \brief Returns solution to use for relaxation.
  virtual const Vector& getRelaxationVector() const
  {
    static Vector empty;
    return empty;
  }

  //! \brief Sets the relaxed solution.
  virtual void setRelaxedSolution(const Vector&) {}

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp, bool firstS1 = true)
  {
    if (tp.multiSteps())
      this->S1.getProcessAdm().cout <<"\n  step="<< tp.step
                                    <<"  time="<< tp.time.t << std::endl;

    int maxit = this->getMaxit(tp.step);
    SIM::ConvStatus conv = SIM::OK;
    for (tp.iter = 0; tp.iter <= maxit && conv != SIM::CONVERGED; tp.iter++)
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

      if (omega0 != 0.0) {
        if (tp.iter > 0) {
          // Calculation of Aitken acceleration factor
          if (aitken) {
            Vector r1 = this->getAitkenResidual();
            r1 -= prevRes;
            omega *= -prevRes.dot(r1) / r1.dot(r1);
            if (fabs(omega) < 1e-6) {
              std::cerr <<"\n  ** Relaxation weight "<< omega <<" too small,"
                        <<" resetting to default "<< omega0 << std::endl;
              omega = omega0;
            }
          }

          // Perform relaxation
          IFEM::cout <<", omega="<< omega;
          prevSol *= 1.0 - omega;
          prevSol.add(this->getRelaxationVector(), omega);
          this->setRelaxedSolution(prevSol);
        } else {
          omega = omega0;
          prevSol = this->getRelaxationVector();
        }
        if (aitken)
          prevRes = this->getAitkenResidual();
      }
      IFEM::cout << std::endl;
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
  //! \brief Parses sub-iteration setup from an XML tag.
  void parseIterations(const TiXmlElement* elem)
  {
    IFEM::cout <<"\tUsing sub-iterations\n";

    std::string func;
    if (utl::getAttribute(elem,"maxFunc",func) ||
        utl::getAttribute(elem,"max",func))
      if (func.find_first_of('t') != std::string::npos)
      {
        IFEM::cout <<"\t\tmax = ";
        maxItFunc = utl::parseTimeFunc(func.c_str(),"expression");
      }

    if (!maxItFunc && utl::getAttribute(elem,"max",maxIter))
      IFEM::cout <<"\t\tmax = "<< maxIter << std::endl;

    if ((utl::getAttribute(elem,"relax",omega) ||
         utl::getAttribute(elem,"omega",omega)) && omega != 0.0) {
      IFEM::cout <<"\t\trelaxation = "<< omega;
      if (utl::getAttribute(elem,"aitken",aitken) && aitken)
        IFEM::cout <<" (aitken)";
      IFEM::cout << std::endl;
      omega0 = omega;
    }
  }

  //! \brief Returns the maximum number of sub-iteration cycles.
  int getMaxit(int iStep = 0) const
  {
    if (noStg)
      return 0; // staggering is disabled
    else if (maxItFunc)
      return static_cast<int>((*maxItFunc)(iStep));

    return maxIter;
  }

private:
  int maxIter; //!< Maximum number of iterations
  ScalarFunc* maxItFunc; //!< Maximum number of iterations as a function
  bool noStg;  //!< If \e true, sub-iterations is disabled

  double omega0; //!< Initial relaxation parameter
  double omega;  //!< Current relaxation parameter
  bool   aitken; //!< If \e true, Aitken-acceleration is enabled

  Vector prevSol; //!< Previous solution for relaxed field
  Vector prevRes; //!< Previous residual used for aitken-acceleration factor
};

#endif
