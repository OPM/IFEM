//==============================================================================
//!
//! \file SIMOverride.h
//!
//! \date Oct 25 2013
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Template overriding methods for wrapping template purposes
//!
//==============================================================================

#ifndef SIM_OVERRIDE_H_
#define SIM_OVERRIDE_H_

#include "Vec3.h"
#include "VTF.h"

class ProcessAdm;

/*! \brief Overrides methods commonly used in applications.
 *  \details This template overrides a subset of SIMbase,
 *           a subset of SIMdependency and parts of the SIMSolver concept.
 *           Used when augmenting classes through a wrapping template.
 */
template<class T> class SIMOverride : public SIMdependency
{
public:
  SIMOverride(T& t) : base(t) {}
  virtual ~SIMOverride() {}

  // SIMbase
  bool read( const char* fileName)
  { 
    return base.read(fileName);
  }

  void setInitialConditions()
  { 
    base.setInitialConditions();
  }

  bool hasIC(const std::string& name) const
  { 
    return base.hasIC(name);
  }

  size_t getNoSpaceDim() const
  { 
    return base.getNoSpaceDim();
  }

  void getBoundaryNodes(int pcode, std::vector<int>& glbNodes,
                        Vec3Vec* XYZ = NULL) const
  { 
    return base.getBoundaryNodes(pcode, glbNodes, XYZ);
  }

  const PatchVec& getFEModel() const
  { 
    return base.getFEModel();
  }

  int getUniquePropertyCode(const std::string& setName, int comp = 0)
  {
    return base.getUniquePropertyCode(setName, comp);
  }

  size_t setVecProperty(int code, Property::Type ptype,
                        VecFunc* field = NULL,
                        int pflag = -1)
  { 
    return base.setVecProperty(code, ptype, field, pflag);
  }

  bool saveModel(char* fileName)
  { 
    return base.saveModel(fileName);
  }

  bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  { 
    return base.saveModel(fileName, geoBlk, nBlock);
  }

  VTF* getVTF() const
  { 
    return base.getVTF();
  }

  //! \brief Returns process administrator 
  const ProcessAdm& getProcessAdm() const { return base.getProcessAdm(); }

  // SIMdependency
  void registerDependency(SIMdependency* sim, const std::string& name,
                          short int nvc, const PatchVec& patches,
                          bool diffBasis = false)
  {
    base.registerDependency(sim, name, nvc, patches, diffBasis);
  }
  void registerDependency(SIMdependency* sim, const std::string& name,
                          short int nvc=1)
  {
    base.registerDependency(sim, name, nvc);
  }

  bool fillField(const std::string& name, const std::vector<double>& values)
  {
    return base.fillField(name, values);
  }
  utl::vector<double>* getField(const std::string& name)
  {
    return base.getField(name);
  }
  const utl::vector<double>* getField(const std::string& name) const
  {
    return base.getField(name);
  }
  const utl::vector<double>* getDependentField(const std::string& name) const
  { 
    return base.getDependentField(name);
  }
  ASMbase* getDependentPatch(const std::string& name, int pindx) const
  {
    return base.getDependentPatch(name, pindx);
  }

  // SIMSolver
  void setupDependencies()
  {
    base.setupDependencies();
  }
  bool init(const TimeStep& tp)
  { 
    return base.init(tp);
  }
  void registerFields(DataExporter& exporter)
  { 
    base.registerFields(exporter);
  }
  bool saveStep(const TimeStep& tp, int& nBlock)
  { 
    return base.saveStep(tp, nBlock);
  }
  bool advanceStep(TimeStep& tp)
  { 
    return base.advanceStep(tp);
  }
protected:
  T& base;
};

#endif
