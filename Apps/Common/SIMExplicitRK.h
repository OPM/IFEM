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

#include "SIMadmin.h"
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
      RK.b.push_back(1.0);
      RK.c.push_back(0.0);
      RK.A.resize(1,1);
    }
    else if (type == HEUN) {
      RK.order = 2;
      RK.b.push_back(0.5);
      RK.b.push_back(0.5);
      RK.c.push_back(0.0);
      RK.c.push_back(1.0);
      RK.A.resize(2,2);
      RK.A(2,1) = 1.0;
    }
    else if (type == RK3) {
      RK.order = 3;
      RK.b.push_back(1.0/6.0);
      RK.b.push_back(2.0/3.0);
      RK.b.push_back(1.0/6.0);
      RK.c.push_back(0.0);
      RK.c.push_back(0.5);
      RK.c.push_back(1.0);
      RK.A.resize(3,3);
      RK.A(2,1) =  0.5;
      RK.A(3,1) = -1.0;
      RK.A(3,2) =  2.0;
    }
    else if (type == RK4) {
      RK.order = 4;
      RK.b.push_back(1.0/6.0);
      RK.b.push_back(1.0/3.0);
      RK.b.push_back(1.0/3.0);
      RK.b.push_back(1.0/6.0);
      RK.c.push_back(0.0);
      RK.c.push_back(0.5);
      RK.c.push_back(0.5);
      RK.c.push_back(1.0);
      RK.A.resize(4,4);
      RK.A(2,1) = 0.5;
      RK.A(3,2) = 0.5;
      RK.A(4,3) = 1.0;
    }
    else
      RK.order = 0;
  }

  //! \brief Returns the parallel process administrator.
  //! \copydoc ISolver::getProcessAdm
  const ProcessAdm& getProcessAdm() const { return solver.getProcessAdm(); }

  //! \copydoc ISolver::solveStep(TimeStep&)
  virtual bool solveStep(TimeStep& tp)
  {
    int msgLevel = 0;
    if (alone)
      solver.getProcessAdm().cout <<"\n  step = "<< tp.step <<"  time = "<< tp.time.t << std::endl;
    else
      std::swap(SIMadmin::msgLevel, msgLevel);

    Vectors stages;
    bool result = this->solveRK(stages, tp);
    if (!alone)
      std::swap(SIMadmin::msgLevel, msgLevel);

    return result;
  }

  //! \brief Applies the Runge-Kutta scheme.
  //! \param stages Vector of stage vectors
  //! \param[in] tp Time stepping information
  bool solveRK(Vectors& stages, const TimeStep& tp)
  {
    TimeDomain time(tp.time);
    Vector dum;

    stages.resize(RK.b.size());

    for (size_t i=0;i<stages.size();++i) {
      Vector tmp(solver.getSolution());
      for (size_t j=0;j<i;++j)
        tmp.add(stages[j], tp.time.dt*RK.A(i+1,j+1));
      time.t = tp.time.t+tp.time.dt*(RK.c[i]-1.0);
      solver.updateDirichlet(time.t, &dum);
      solver.applyDirichlet(tmp);
      if (!solver.assembleSystem(time, Vectors(1, tmp),
                                 !linear || (tp.step == 1 && i == 0)))
        return false;

      // solve Mu = Au + f
      if (!solver.solveSystem(stages[i]))
        return false;

      if (linear)
        solver.setMode(SIM::RHS_ONLY);
    }

    // finally construct solution as weighted stages
    solver.updateDirichlet(tp.time.t, &dum);
    for (size_t i=0;i<RK.b.size();++i)
      solver.getSolution().add(stages[i], tp.time.dt*RK.b[i]);
    solver.applyDirichlet(solver.getSolution());

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
  bool linear = false; //!< If true mass matrix is constant
};

}

#endif
