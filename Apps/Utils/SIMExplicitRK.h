//==============================================================================
//!
//! \file SIMExplicitRK.h
//!
//! \date Nov 28 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Explicit Runge-Kutta based time stepping for SIM classes
//!
//==============================================================================

#ifndef SIM_EXPLICIT_RK_H_
#define SIM_EXPLICIT_RK_H_

#include "ProcessAdm.h"
#include "SIMenums.h"
#include "TimeIntUtils.h"
#include "TimeStep.h"

class DataExporter;


namespace TimeIntegration {

//! \brief Explicit Runge-Kutta based time stepping for SIM classes
//! \details Template can be instanced over any SIM implementing ISolver,
//!          and which derive from SIMbase.
template<class Solver>
class SIMExplicitRK
{
public:
  //! \brief Constructor
  //! \param solv The simulator to do time stepping for
  //! \param type The Runge-Kutta scheme to use
  //! \param standalone If true, this is a standalone solver
  SIMExplicitRK(Solver& solv, Method type, bool standalone = true) :
    solver(solv), alone(standalone)
  {
    if (type == EULER) {
      RK.order = 1;
      RK.b = {1.0};
      RK.c = {0.0};
      RK.A.resize(1,1);
    }
    else if (type == HEUN) {
      RK.order = 2;
      RK.b = {0.5, 0.5};
      RK.c = {0.0, 1.0};
      RK.A.resize(2,2);
      RK.A(2,1) = 1.0;
    }
    else if (type == RK3) {
      RK.order = 3;
      RK.b = {1.0/6.0, 2.0/3.0, 1.0/6.0};
      RK.c = {0.0, 0.5, 1.0};
      RK.A.resize(3,3);
      RK.A(2,1) =  0.5;
      RK.A(3,1) = -1.0;
      RK.A(3,2) =  2.0;
    }
    else if (type == RK4) {
      RK.order = 4;
      RK.b = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
      RK.c = {0.0, 0.5, 0.5, 1.0};
      RK.A.resize(4,4);
      RK.A(2,1) = 0.5;
      RK.A(3,2) = 0.5;
      RK.A(4,3) = 1.0;
    }
    else
      RK.order = 0;

    Solver::msgLevel = 1; // prints primary solution summary only
  }

  //! \copydoc ISolver::getProcessAdm()
  const ProcessAdm& getProcessAdm() const { return solver.getProcessAdm(); }

  //! \copydoc ISolver::solveStep(TimeStep&)
  virtual bool solveStep(TimeStep& tp)
  {
    if (alone)
      solver.getProcessAdm().cout <<"\n  step = "<< tp.step <<"  time = "<< tp.time.t << std::endl;

    solver.setMode(this->assemble ? SIM::DYNAMIC : SIM::RHS_ONLY);

    if (!solver.initDirichlet(tp.time.t))
      return false;

    Vectors stages;
    return this->solveRK(stages, tp);
  }

  //! \brief Applies the Runge-Kutta scheme.
  //! \param stages Vector of stage vectors
  //! \param[in] tp Time stepping information
  bool solveRK(Vectors& stages, const TimeStep& tp)
  {
    TimeDomain time(tp.time);
    Vector dum;

    stages.resize(RK.b.size());

    for (size_t i = 0; i < stages.size(); ++i) {
      Vector tmp(solver.getSolution());
      for (size_t j = 0; j < i; ++j)
        tmp.add(stages[j], tp.time.dt*RK.A(i+1,j+1));
      time.t = tp.time.t+tp.time.dt*(RK.c[i]-1.0);

      // 1) Update to stage time and impose stage-state boundary values.
      //    Passing an empty vector means absolute Dirichlet values.
      if (!solver.updateDirichlet(time.t, &dum))
        return false;

      if (!solver.applyDirichlet(tmp))
        return false;

      // 2) Constrain stage unknowns to tangent boundary values dg/dt.
      Vector zeroRef(solver.getSolution().size());
      if (!solver.updateDirichlet(time.t, &zeroRef, true))
        return false;

      if (!solver.assembleSystem(time, Vectors(1, tmp),
                                 !linear || (tp.step == 1 && i == 0)))
        return false;

      // solve Mk = Au + f
      if (!solver.solveSystem(stages[i]))
        return false;

      if (linear) {
        solver.setMode(SIM::RHS_ONLY);
        assemble = false;
      }
    }

    // finally construct solution as weighted stages
    for (size_t i = 0; i < RK.b.size(); ++i)
      solver.getSolution().add(stages[i], tp.time.dt*RK.b[i]);

    // Enforce absolute boundary values at the new time level.
    if (!solver.updateDirichlet(tp.time.t, &dum))
      return false;
    if (!solver.applyDirichlet(solver.getSolution()))
      return false;

    if (alone)
      solver.printSolutionSummary(solver.getSolution(), 0,
                                  solver.getProblem()->getField1Name(1).c_str());

    return true;
  }

  //! \copydoc ISolver::advanceStep(TimeStep&)
  bool advanceStep(TimeStep& tp)
  {
    return solver.advanceStep(tp);
  }

  //! \copydoc ISolver::saveModel(char*,int&,int&)
  bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    return solver.saveModel(fileName, geoBlk, nBlock);
  }

  //! \copydoc ISolver::saveStep(const TimeStep&,int&)
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    return solver.saveStep(tp, nBlock);
  }

  //! \copydoc ISolver::registerFields(DataExporter&)
  void registerFields(DataExporter& exporter)
  {
    solver.registerFields(exporter);
  }

  //! \brief Serialize internal state for restarting purposes.
  //! \param data Container for serialized data
  bool serialize(std::map<std::string,std::string>& data)
  {
    return solver.serialize(data);
  }

  //! \brief Set internal state from a serialized state.
  //! \param[in] data Container for serialized data
  bool deSerialize(const std::map<std::string,std::string>& data)
  {
    return solver.deSerialize(data);
  }

  //! \brief Mark operator as linear to avoid repeated assembly and factorization.
  void setLinear(bool enable) { linear = enable; }

protected:
  Solver& solver; //!< Reference to simulator
  RKTableaux RK;  //!< Tableaux of Runge-Kutta coefficients
  bool alone; //!< If true, this is a standalone solver
  bool linear = false; //!< If true operators are constant
  bool assemble = true; //!< If true, assemble operators
};

}

#endif
