// $Id$
//==============================================================================
//!
//! \file SIMCoupled.h
//!
//! \date Oct 12 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Coupled SIM solver class template.
//!
//==============================================================================

#ifndef _SIM_COUPLED_H_
#define _SIM_COUPLED_H_

#include "Function.h"
#include "Property.h"

class SIMdependency;
class ASMbase;
class DataExporter;
class TimeStep;
class VTF;


/*!
  \brief Template class for coupled simulators.
*/

template<class T1, class T2> class SIMCoupled
{
public:
  //! \brief The constructor initializes the references to the two simulators.
  SIMCoupled(T1& s1, T2& s2) : S1(s1), S2(s2) {}

  //! \brief The destructor nullifies the VTF pointer for the second simulator.
  virtual ~SIMCoupled() { S2.setVTF(nullptr); }

  //! \brief Sets up field dependencies.
  virtual void setupDependencies() {}

  //! \brief Performs some pre-processing tasks on the FE model.
  bool preprocess()
  {
    return S1.preprocess() && S2.preprocess();
  }

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& tp)
  {
    return S1.advanceStep(tp) && S2.advanceStep(tp);
  }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp, bool firstS1 = true)
  {
    if (firstS1)
      return S1.solveStep(tp) && S2.solveStep(tp);
    else
      return S2.solveStep(tp) && S1.solveStep(tp);
  }

  //! \brief Postprocesses the solution of current time step.
  bool postSolve(const TimeStep& tp, bool restart = false)
  {
    return S1.postSolve(tp,restart) && S2.postSolve(tp,restart);
  }

  //! \brief Saves the converged results to VTF-file of a given time step.
  virtual bool saveStep(const TimeStep& tp, int& nBlock)
  {
    return S2.saveStep(tp,nBlock) && S1.saveStep(tp,nBlock);
  }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  virtual bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    if (!S1.saveModel(fileName,geoBlk,nBlock))
      return false;

    S2.setVTF(S1.getVTF());
    return true;
  }

  //! \brief Returns the current VTF-file object.
  VTF* getVTF() const { return S1.getVTF(); }

  //! \brief Initializes for time-dependent simulation.
  virtual bool init(const TimeStep& tp)
  {
    return S1.init(tp) && S2.init(tp);
  }

  //! \brief Registers a dependency on a field from another SIM object.
  virtual void registerDependency(SIMdependency* sim, const std::string& name,
                                  short int nvc,
                                  const std::vector<ASMbase*>& patches,
                                  char diffBasis = 0, int component = 1)
  {
    S1.registerDependency(sim, name, nvc, patches, diffBasis, component);
    S2.registerDependency(sim, name, nvc, patches, diffBasis, component);
  }

  //! \brief Registers a dependency on a field from another SIM object.
  virtual void registerDependency(SIMdependency* sim, const std::string& name,
                                  short int nvc = 1)
  {
    S1.registerDependency(sim, name, nvc);
    S2.registerDependency(sim, name, nvc);
  }

  //! \brief Returns a unique integer code for a Property set.
  int getUniquePropertyCode(const std::string& setName, int comp = 0)
  {
    return S1.getUniquePropertyCode(setName, comp);
  }

  //! \brief Creates a set of Property objects.
  bool createPropertySet(const std::string& setName, int pc)
  {
    return S1.createPropertySet(setName, pc);
  }

  //! \brief Defines a vector field property.
  size_t setVecProperty(int code, Property::Type ptype,
                        VecFunc* field = nullptr, int pflag = -1)
  {
    return S1.setVecProperty(code, ptype, field, pflag);
  }

  //! \brief Registers the field vectors for storage on HDF5 output.
  void registerFields(DataExporter& exporter)
  {
    S1.registerFields(exporter);
    S2.registerFields(exporter);
  }

  //! \brief Sets the initial conditions for the simulators.
  bool setInitialConditions()
  {
    return S1.setInitialConditions() && S2.setInitialConditions();
  }

  //! \brief Checks whether a named initial condition is present.
  bool hasIC(const std::string& name) const
  {
    return S1.hasIC(name) || S2.hasIC(name);
  }

  //! \brief Returns the nodal vector of named field in this SIM.
  utl::vector<double>* getField(const std::string& name)
  {
    utl::vector<double>* result = S1.getField(name);
    if (!result)
      result = S2.getField(name);

    return result;
  }

protected:
  T1& S1; //!< First substep
  T2& S2; //!< Second substep
};

#endif
