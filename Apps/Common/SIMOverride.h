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
  //! \brief Constructor
  //! \param t The simulator to wrap
  SIMOverride(T& t) : base(t) {}

  //! \brief Empty destructor
  virtual ~SIMOverride() {}

  //! \copydoc SIMbase::read(const char*)
  bool read( const char* fileName)
  { 
    return base.read(fileName);
  }

  //! \brief Set initial conditions
  void setInitialConditions()
  { 
    base.setInitialConditions();
  }

  //! \copydoc SIMdependency::hasIC(const std::string&) const
  bool hasIC(const std::string& name) const
  { 
    return base.hasIC(name);
  }

  //! \copydoc SIMbase::getNoSpaceDim()
  size_t getNoSpaceDim() const
  { 
    return base.getNoSpaceDim();
  }

  //! \copydoc SIMbase::getBoundaryNodes(int,std::vector<int>&,Vec3Vec*) const
  void getBoundaryNodes(int pcode, std::vector<int>& glbNodes,
                        Vec3Vec* XYZ = NULL) const
  { 
    return base.getBoundaryNodes(pcode, glbNodes, XYZ);
  }

  //! \copydoc SIMbase::getFEModel() const
  const PatchVec& getFEModel() const
  { 
    return base.getFEModel();
  }

  //! \copydoc SIMbase::getUniquePropertyCode(const std::string&,int)
  int getUniquePropertyCode(const std::string& setName, int comp = 0)
  {
    return base.getUniquePropertyCode(setName, comp);
  }

  //! \copydoc SIMbase::setVecProperty(int,Property::Type,VecFunc*,int)
  size_t setVecProperty(int code, Property::Type ptype,
                        VecFunc* field = NULL,
                        int pflag = -1)
  { 
    return base.setVecProperty(code, ptype, field, pflag);
  }

  //! \brief Open VTF and save model
  //! \param fileName The filename of the VTF file
  bool saveModel(char* fileName)
  { 
    return base.saveModel(fileName);
  }

  //! \brief Save model to an possibly already opened VTF file
  //! \param fileName The filename of the VTF file
  //! \param[out] geoBlk The starting geometry part for this model
  //! \param[out] nBlock Running VTF block counter
  bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  { 
    return base.saveModel(fileName, geoBlk, nBlock);
  }

  //! \copydoc SIMoutput::getVTF() const
  VTF* getVTF() const
  { 
    return base.getVTF();
  }

  //! \copydoc SIMbase::getProcessAdm() const
  const ProcessAdm& getProcessAdm() const { return base.getProcessAdm(); }

  //! \copydoc SIMdependency::registerDependency(SIMdependency*,const std::string&,short int,const PatchVec&,bool)
  void registerDependency(SIMdependency* sim, const std::string& name,
                          short int nvc, const PatchVec& patches,
                          bool diffBasis = false)
  {
    base.registerDependency(sim, name, nvc, patches, diffBasis);
  }

  //! \copydoc SIMdependency::registerDependency(SIMdependency*,const std::string&,short int)
  void registerDependency(SIMdependency* sim, const std::string& name,
                          short int nvc=1)
  {
    base.registerDependency(sim, name, nvc);
  }

  //! \copydoc SIMdependency::fillField(const std::string&,const std::vector<double>&)
  bool fillField(const std::string& name, const std::vector<double>& values)
  {
    return base.fillField(name, values);
  }

  //! \copydoc SIMdependency::getField(const std::string&)
  utl::vector<double>* getField(const std::string& name)
  {
    return base.getField(name);
  }

  //! \copydoc SIMdependency::getField(const std::string&) const
  const utl::vector<double>* getField(const std::string& name) const
  {
    return base.getField(name);
  }

  //! \copydoc SIMdependency::getDependentField(const std::string&) const 
  const utl::vector<double>* getDependentField(const std::string& name) const
  { 
    return base.getDependentField(name);
  }

  //! \copydoc SIMdependency::getDependentPatch(const std::string&,int) const
  ASMbase* getDependentPatch(const std::string& name, int pindx) const
  {
    return base.getDependentPatch(name, pindx);
  }

  //! \brief Setup inter-SIM dependencies
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
protected:
  T& base; //!< Reference to wrapper simulator
};

#endif
