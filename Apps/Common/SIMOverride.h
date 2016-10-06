// $Id$
//==============================================================================
//!
//! \file SIMOverride.h
//!
//! \date Oct 25 2013
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Template overriding methods for wrapping template purposes.
//!
//==============================================================================

#ifndef SIM_OVERRIDE_H_
#define SIM_OVERRIDE_H_

#include "SIMdependency.h"
#include "Property.h"
#include "Function.h"

class ProcessAdm;
class VTF;


/*!
  \brief Overrides methods commonly used in applications.

  \details This template overrides a subset of SIMbase and SIMdependency,
  and parts of the SIMSolver concept.
  Used when augmenting classes through a wrapping template.
*/

template<class T> class SIMOverride : public SIMdependency
{
  T& base; //!< Reference to wrapped simulator

public:
  //! \brief The constructor initializes the reference to the wrapped simulator.
  SIMOverride(T& t) : base(t) {}
  //! \brief Empty destructor.
  virtual ~SIMOverride() {}

  //! \copydoc SIMadmin::read(const char*)
  bool read(const char* fileName) { return base.read(fileName); }

  //! \copydoc SIMinput::setInitialConditions()
  bool setInitialConditions() { return base.setInitialConditions(); }
  //! \copydoc SIMinput::hasIC(const std::string&) const
  bool hasIC(const std::string& name) const { return base.hasIC(name); }

  //! \copydoc SIMdependency::getNoSpaceDim() const
  virtual size_t getNoSpaceDim() const { return base.getNoSpaceDim(); }
  //! \copydoc SIMdependency::getName() const
  virtual std::string getName() const { return base.getName(); }

  //! \copydoc SIMbase::getFEModel() const
  const PatchVec& getFEModel() const { return base.getFEModel(); }

  //! \copydoc ISolver::saveModel(char*,int&,int&)
  bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    return base.saveModel(fileName, geoBlk, nBlock);
  }

  //! \copydoc SIMoutput::getVTF() const
  VTF* getVTF() const { return base.getVTF(); }

  //! \copydoc SIMbase::getProcessAdm() const
  const ProcessAdm& getProcessAdm() const { return base.getProcessAdm(); }

  //! \copydoc SIMdependency::registerDependency(SIMdependency*,const std::string&,short int,const PatchVec&,char)
  void registerDependency(SIMdependency* sim, const std::string& name,
                          short int nvc, const PatchVec& patches,
                          char diffBasis = 0, int component = 1) override
  {
    base.registerDependency(sim, name, nvc, patches, diffBasis, component);
  }

  //! \copydoc SIMdependency::registerDependency(SIMdependency*,const std::string&,short int)
  virtual void registerDependency(SIMdependency* sim, const std::string& name,
                                  short int nvc = 1)
  {
    base.registerDependency(sim, name, nvc);
  }

  //! \copydoc SIMdependency::fillField(const std::string&,const std::vector<double>&)
  bool fillField(const std::string& name, const std::vector<double>& values)
  {
    return base.fillField(name, values);
  }

  //! \copydoc SIMdependency::getField(const std::string&)
  virtual utl::vector<double>* getField(const std::string& name)
  {
    return base.getField(name);
  }

  //! \copydoc SIMdependency::getField(const std::string&) const
  virtual const utl::vector<double>* getField(const std::string& name) const
  {
    return const_cast<const T&>(base).getField(name);
  }

  //! \copydoc SIMdependency::getDependentField(const std::string&) const
  const utl::vector<double>* getDependentField(const std::string& name) const
  {
    return base.getDependentField(name);
  }

  //! \copydoc ISolver::setupDependencies()
  void setupDependencies()
  {
    base.setupDependencies();
  }

  //! \copydoc ISolver::init(const TimeStep&)
  bool init(const TimeStep& tp)
  {
    return base.init(tp);
  }

  //! \copydoc ISolver::registerFields(DataExporter&)
  void registerFields(DataExporter& exporter)
  {
    base.registerFields(exporter);
  }

  //! \copydoc ISolver::saveStep(const TimeStep&,int&)
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    return base.saveStep(tp, nBlock);
  }

  //! \copydoc ISolver::advanceStep(TimeStep&)
  bool advanceStep(TimeStep& tp)
  {
    return base.advanceStep(tp);
  }

  T& get() { return base; }
  const T& get() const { return base; }

  int getGlobalNode(int node) const { return base.getGlobalNode(node); }
  int getLocalNode(int node) const { return base.getLocalNode(node); }
};

#endif
