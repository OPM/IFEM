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

#include "MatVec.h"
#include "SIMCoupled.h"
#include "SIMenums.h"
#include "TimeStep.h"


/*!
  \brief Template class for semi-implicitly coupled simulators.
*/

template<class T1, class T2>
class SIMCoupledSI : public SIMCoupled<T1,T2>
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  SIMCoupledSI(T1& s1, T2& s2) : SIMCoupled<T1,T2>(s1,s2), maxIter(-1)
  {
    omega = omega0 = 0.0;
  }

  //! \brief Empty destructor.
  virtual ~SIMCoupledSI() {}

  //! \brief Enable/disable the staggering iteration cycles.
  virtual void enableStaggering(bool enable = true)
  {
    maxIter = enable ? std::min(this->S1.getMaxit(),this->S2.getMaxit()) : 0;
  }

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

  //! \brief Set the relaxed solution.
  virtual void setRelaxedSolution(const Vector&) {}

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

      // Calculate aitken acceleration factor
      if (aitken && omega0 != 0.0) {
        if (tp.iter > 0) {
            Vector r1 = this->getAitkenResidual();
            r1 -= prevRes;
            omega *= -prevRes.dot(r1) / r1.dot(r1);
        }
        prevRes = this->getAitkenResidual();
      }

      // Perform relaxation
      if (omega0 != 0.0) {
        if (tp.iter > 0) {
          IFEM::cout << "  relaxing field update, omega=" << omega << std::endl;
          prevSol *= 1.0 - omega;
          prevSol.add(this->getRelaxationVector(), omega);
          this->setRelaxedSolution(prevSol);
        } else {
          omega = omega0;
          prevSol = this->getRelaxationVector();
        }
      }
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
  double omega; //!< Relaxation parameter
  double omega0; //!< Initial relaxation parameter
  bool aitken; //!< True to enable aitken-acceleration
  Vector prevSol; //!< Previous solution for relaxed field
  Vector prevRes; //!< Previous residual used for aitken-acceleration factor
};

#endif
