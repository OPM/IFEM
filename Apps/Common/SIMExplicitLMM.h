//==============================================================================
//!
//! \file SIMExplicitLMM.h
//!
//! \date Aug 20 2019
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Explicit linear multistep time stepping for SIM classes.
//!
//==============================================================================

#ifndef SIM_EXPLICIT_LMM_H_
#define SIM_EXPLICIT_LMM_H_

#include "SystemMatrix.h"
#include "SIMenums.h"
#include "TimeIntUtils.h"
#include "TimeStep.h"

class DataExporter;


namespace TimeIntegration {

  //! \brief Explicit linear multistep time stepping for SIM classes.
  //! \details Template can be instanced over any SIM implementing ISolver,
  //            and which derive from SIMbase.
  template<class Solver>
class SIMExplicitLMM
{
public:
  //! \brief Constructor.
  //! \param solv The simulator to do time stepping for
  //! \param type The linear multistep scheme to use
  //! \param standalone If true, this is a standalone solver
  //! \param solField Name of primary solution fields (for ICs)
  SIMExplicitLMM(Solver& solv, Method type, bool standalone = true,
                 const std::string& solField = "") :
    solver(solv), alone(standalone), fieldName(solField)
  {
    if (type == AB2)
      order = 2;
    else if (type == AB3)
      order = 3;
    else if (type == AB4)
      order = 4;
    else if (type == AB5)
      order = 5;
    else
      order = 1;

    loads.resize(order, nullptr);
  }

  //! \brief Destructor frees up the load vectors.
  ~SIMExplicitLMM()
  {
    for (auto& v : loads)
      delete v;
  }

  //! \copydoc ISolver::solveStep(TimeStep&)
  bool solveStep(TimeStep& tp)
  {
    if (alone)
      solver.getProcessAdm().cout <<"\n  step = "<< tp.step <<"  time = "<< tp.time.t << std::endl;

    TimeDomain time(tp.time);
    time.t = tp.time.t - tp.time.dt;
    if (!solver.assembleSystem(time, Vectors(1, solver.getSolution(1)),
                               !linear || (tp.step == 1)))
      return false;

    loads[0] = solver.getRHSvector(0, true);

    const std::vector<std::vector<double>> AB_coefs = 
      {{1.0},
       {-0.5, 1.5},
       {5.0/12.0, -16.0/12.0, 23.0/12.0},
       {-9.0/24.0, 37.0/24, -59.0/24.0, 55.0/24.0},
       {251.0/720.0, -1274.0/720.0, 2616.0/720.0, -2774.0/720.0, 1901.0/720.0}};

    int c_order = hasICs? order-1 : std::min(order-1, tp.step-1);
    const auto& AB_coef = AB_coefs[c_order];
    solver.getRHSvector(0, false)->mult(AB_coef.back() * tp.time.dt);
    for (size_t j = 0; j < AB_coef.size()-1; ++j)
      solver.addToRHSvector(0, *loads[c_order-j], AB_coef[j]*tp.time.dt);

    if (!solver.solveSystem(solver.getSolution()))
      return false;

    if (linear)
      solver.setMode(SIM::RHS_ONLY);

    solver.getSolution() += solver.getSolution(1);

    Vector dum;
    solver.updateDirichlet(tp.time.t, &dum);
    solver.applyDirichlet(solver.getSolution());

    if (alone)
      solver.printSolutionSummary(solver.getSolution(), 0,
                                  solver.getProblem()->getField1Name(1).c_str());

    return true;
  }

  //! \copydoc ISolver::advanceStep(TimeStep&)
  bool advanceStep(TimeStep& tp)
  {
    // Evaluate fluxes for initial conditions.
    if (tp.step == 1) {
      hasICs = true;
      for (int j = 2; j <= order; ++j) {
        std::stringstream str;
        str << fieldName << j;
        if (solver.hasIC(str.str())) {
          TimeDomain time(tp.time);
          time.t = tp.time.t - j*tp.time.dt;
          if (!solver.assembleSystem(time, Vectors(1, solver.getSolution(j-1))))
            return false;

          loads[j-2] = solver.getRHSvector(0, true);
        } else {
          hasICs = false;
          break;
        }
      }
    }

    delete loads.back();
    for (int j = order-2; j >= 0; --j)
      loads[j+1] = loads[j];

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
  std::vector<SystemVector*> loads; //!< Unscaled load vectors
  int order; //!< Order of method
  bool alone; //!< If true, this is a standalone solver
  const std::string fieldName; //!< Name of primary solution fields (for ICs)
  bool hasICs = false; //!< If true, start with full order
  bool linear = false; //!< If true, mass matrix is constant
};

}

#endif
