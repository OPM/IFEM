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

#include "DataExporter.h"
#include "TimeIntUtils.h"
#include "TimeStep.h"

namespace TimeIntegration {

  //! \brief Explicit Runge-Kutta based time stepping for SIM classes
  //! \details Template can be instanced over any SIM implementing ISolver,
  //            and which derive from SIMbase.
  template<class Solver>
class SIMExplicitRK
{
public:
  //! \brief Constructor
  //! \param solv The simulator to do time stepping for
  //! \param type The Runge-Kutta scheme to use
  SIMExplicitRK(Solver& solv, Method type) : solver(solv)
  {
    if (type == EULER) {
      RK.order = 1;
      RK.b.push_back(1.0);
      RK.c.push_back(0.0);
      RK.A.redim(1,1);
    }
    if (type == HEUN) {
      RK.order = 2;
      RK.b.push_back(0.5);
      RK.b.push_back(0.5);
      RK.c.push_back(0.0);
      RK.c.push_back(1.0);
      RK.A.redim(2,2);
      RK.A(2,1) = 1.0;
    }
    if (type == RK3) {
      RK.order = 3;
      RK.b.push_back(1.0/6.0);
      RK.b.push_back(2.0/3.0);
      RK.b.push_back(1.0/6.0);
      RK.c.push_back(0.0);
      RK.c.push_back(0.5);
      RK.c.push_back(1.0);
      RK.A.redim(3,3);
      RK.A(2,1) =  0.5;
      RK.A(3,1) = -1.0;
      RK.A(3,2) =  2.0;
    }
    if (type == RK4) {
      RK.order = 4;
      RK.b.push_back(1.0/6.0);
      RK.b.push_back(1.0/3.0);
      RK.b.push_back(1.0/3.0);
      RK.b.push_back(1.0/6.0);
      RK.c.push_back(0.0);
      RK.c.push_back(0.5);
      RK.c.push_back(0.5);
      RK.c.push_back(1.0);
      RK.A.redim(4,4);
      RK.A(2,1) = 0.5;
      RK.A(3,2) = 0.5;
      RK.A(4,3) = 1.0;
    }
  }

  //! \copydoc ISolver::solveStep(TimeStep&)
  bool solveStep(TimeStep& tp)
  {
    std::cout <<"\n  step = "<< tp.step <<"  time = "<< tp.time.t << std::endl;

    std::vector<Vector> stages;
    return solveRK(stages, tp);

    return true;
  }

  //! \brief Apply the Runge-Kutta scheme
  //! \param stages Vector of stage vectors
  //! \param tp Time stepping information
  bool solveRK(std::vector<Vector>& stages, TimeStep& tp)
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
      if (!solver.assembleSystem(time, Vectors(1, tmp)))
        return false;

      // solve Mu = Au + f
      if (!solver.solveSystem(stages[i]))
        return false;
    }

    // finally construct solution as weighted stages
    solver.updateDirichlet(tp.time.t, &dum);
    for (size_t i=0;i<RK.b.size();++i)
      solver.getSolution().add(stages[i], tp.time.dt*RK.b[i]);
    solver.applyDirichlet(solver.getSolution());

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

  //! \brief Serialize internal state for restarting purposes.
  //! \param data Container for serialized data
  bool serialize(DataExporter::SerializeData& data)
  {
    return solver.serialize(data);
  }

  //! \brief Set internal state from a serialized state.
  //! \param[in] data Container for serialized data
  bool deSerialize(const DataExporter::SerializeData& data)
  {
    return solver.deSerialize(data);
  }

protected:
  Solver& solver; //!< Reference to simulator
  RKTableaux RK;  //!< Tableaux of Runge-Kutta coefficients
};

}

#endif
