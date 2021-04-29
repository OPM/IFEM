//==============================================================================
//!
//! \file SIMImplicitLMM.h
//!
//! \date Aug 22 2019
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Implicit linear multistep time stepping for SIM classes.
//!
//==============================================================================

#ifndef SIM_IMPLICIT_LMM_H_
#define SIM_IMPLICIT_LMM_H_

#include "NonLinSIM.h"
#include "SystemMatrix.h"
#include "SIMenums.h"
#include "TimeIntUtils.h"
#include "TimeStep.h"

class DataExporter;


namespace TimeIntegration {

  //! \brief Implicit linear multistep time stepping for SIM classes.
  //! \details Template can be instanced over any SIM implementing ISolver,
  //            and which derive from SIMbase.
  template<class Solver>
class SIMImplicitLMM
{
public:
  //! \brief Constructor.
  //! \param solv The simulator to do time stepping for
  //! \param type The linear multistep scheme to use
  //! \param standalone If true, this is a standalone solver
  //! \param solField Name of primary solution fields (for ICs)
  SIMImplicitLMM(Solver& solv, Method type, bool standalone = true,
                 const std::string& solField = "") :
    solver(solv), nSim(solver, loads),
    alone(standalone), fieldName(solField)
  {
    if (type == AM2)
      order = 2;
    else if (type == AM3)
      order = 3;
    else if (type == AM4)
      order = 4;
    else
      order = 1;

    loads.resize(order, nullptr);
  }

  //! \brief Destructor frees up the load vectors.
  ~SIMImplicitLMM()
  {
    for (SystemVector* v : loads)
      delete v;
  }

  //! \brief Returns the parallel process administrator.
  //! \copydoc ISolver::getProcessAdm
  const ProcessAdm& getProcessAdm() const { return solver.getProcessAdm(); }

  //! \copydoc ISolver::solveStep(TimeStep&)
  bool solveStep(TimeStep& tp)
  {
    if (alone)
      solver.getProcessAdm().cout <<"\n  step = "<< tp.step <<"  time = "<< tp.time.t << std::endl;

    const std::vector<std::vector<double>> AM_coefs =
      {{1.0},
       {0.5, 0.5},
       {-1.0/12.0, 2.0/3.0, 5.0/12.0},
       {1.0/24.0, -5.0/24, 19.0/24.0, 9.0/24.0},
       {-19.0/720.0, 106.0/720.0, -264.0/720.0, 646.0/720.0, 251.0/720.0}};

    int c_order = hasICs? order-1 : std::min(order-1, tp.step-1);
    nSim.setCoefs(AM_coefs[c_order]);

    nSim.initSol(2);
    nSim.theSolutions()[0] = solver.getSolution();
    nSim.theSolutions()[1] = solver.getSolution();
    solver.setTimeScale(AM_coefs[c_order].back()*tp.time.dt);
    if (nSim.solveStep(tp, SIM::DYNAMIC) != SIM::CONVERGED)
      return false;

    solver.getSolution() = nSim.getSolution(0);

    solver.setMode(SIM::RHS_ONLY);
    solver.setTimeScale(1.0);
    if (!solver.assembleSystem(tp.time, Vectors(1, solver.getSolution()), false))
      return false;

    loads[0] = solver.getRHSvector(0, true);

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
  //! \brief Specialized nonlinear solver for implicit LMM methods.
  class LMMNonLinSIM : public NonLinSIM
  {
  public:
    //! \brief The constructor initializes default solution parameters.
    //! \param sim Pointer to the spline FE model
    //! \param load Vector of system load vectors
    //! \param[in] n Which type of iteration norm to use in convergence checks
    LMMNonLinSIM(SIMbase& sim,
                 std::vector<SystemVector*>& load,
                 CNORM n = ENERGY) :
      NonLinSIM(sim, n), loads(load) {}

    //! \brief Set scaling coefficients
    void setCoefs(const std::vector<double>& coef) { coefs = coef; }

  protected:
    //! \brief Administers assembly of the linear equation system.
    //! \param[in] time Parameters for nonlinear/time-dependent simulations
    //! \param[in] pSol Previous primary solution vectors in DOF-order
    //! \param[in] newLHSmatrix If \e false, only integrate the RHS vector
    //! \param[in] poorConvg If \e true, the nonlinear driver is converging poorly
    bool assembleSystem(const TimeDomain& time, const Vectors& pSol,
                        bool newLHSmatrix = true, bool poorConvg = false)
    {
      if (!model.assembleSystem(time, pSol, newLHSmatrix, poorConvg))
        return false;

      for (size_t i = 0; i < coefs.size()-1; ++i)
        model.addToRHSvector(0, *loads[coefs.size()-i-1], coefs[i]*time.dt);

      return true;
    }

    std::vector<SystemVector*>& loads; //!< Reference to load vectors
    std::vector<double> coefs; //!< Coefficients
  };

  Solver& solver; //!< Reference to simulator
  LMMNonLinSIM nSim; //!< Nonlinear solver
  std::vector<SystemVector*> loads; //!< Unscaled load vectors
  int order; //!< Order of method
  bool alone; //!< If true, this is a standalone solver
  const std::string fieldName; //!< Name of primary solution fields (for ICs)
  bool hasICs = false; //!< If true, start with full order
  bool linear = false; //!< If true, mass matrix is constant
};

}

#endif
