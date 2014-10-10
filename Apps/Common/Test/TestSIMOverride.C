//==============================================================================
//!
//! \file TestSIMOverride.C
//!
//! \date Oct 13 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for template overriding methods for wrapping template purposes.
//!
//==============================================================================

#include "DataExporter.h"
#include "ProcessAdm.h"
#include "TimeStep.h"
#include "SIMOverride.h"

#include "gtest/gtest.h"

class SIMMockOverride {
public:
  SIMMockOverride() :
    read_called(false),
    setinitialconditions_called(false),
    hasic_called(false),
    getnospacedim_called(false),
    getfemodel_called(false),
    savemodel_called(false),
    getvtf_called(false),
    getprocessadm_called(false),
    registerdependency_called(false),
    registerdependency2_called(false),
    fillfield_called(false),
    getfield_called(false),
    getfield2_called(false),
    getdependentfield_called(false),
    setupdependencies_called(false),
    init_called(false),
    registerfields_called(false),
    savestep_called(false),
    advancestep_called(false),
    postsolve_called(false),
    getglobalnode_called(false),
    getlocalnode_called(false)
  {}

  bool read(const char* fileName) { read_called = true; return true;}
  void setInitialConditions() { setinitialconditions_called = true; }
  bool hasIC(const std::string& name) const { hasic_called = true; return true; }
  size_t getNoSpaceDim() const { getnospacedim_called = true; return 2;}
  const std::vector<ASMbase*>& getFEModel() const
  {
    static std::vector<ASMbase*> result;
    getfemodel_called = true;
    return result;
  }
  bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    savemodel_called = true;
    return true;
  }
  VTF* getVTF() const { getvtf_called = true; return NULL; }
  const ProcessAdm& getProcessAdm() const
  {
    static ProcessAdm result;
    getprocessadm_called = true;
    return result;
  }
  void registerDependency(SIMdependency* sim, const std::string& name,
                          short int nvc, const std::vector<ASMbase*>& patches,
                          bool diffBasis = false)
  {
    registerdependency_called = true;
  }
  void registerDependency(SIMdependency* sim, const std::string& name,
                          short int nvc = 1)
  {
    registerdependency2_called = true;
  }
  bool fillField(const std::string& name, const std::vector<double>& values)
  {
    fillfield_called = true;
    return true;
  }
  utl::vector<double>* getField(const std::string& name)
  {
    getfield_called = true;
    return NULL;
  }
  const utl::vector<double>* getField(const std::string& name) const
  {
    getfield2_called = true;
    return NULL;
  }
  const utl::vector<double>* getDependentField(const std::string& name) const
  {
    getdependentfield_called = true;
    return NULL;
  }
  void setupDependencies() { setupdependencies_called = true; }
  bool init(const TimeStep& tp) { init_called = true; return true; }
  void registerFields(DataExporter& exporter) { registerfields_called = true; }
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    savestep_called = true;
    return true;
  }
  bool advanceStep(TimeStep& tp) { advancestep_called = true; return true; }
  void postSolve(const TimeStep& tp) { postsolve_called = true; }
  int getGlobalNode(int node) const { getglobalnode_called = true; return 0; }
  int getLocalNode(int node) const { getlocalnode_called = true; return 0; }

  bool read_called;
  bool setinitialconditions_called;
  mutable bool hasic_called;
  mutable bool getnospacedim_called;
  mutable bool getfemodel_called;
  bool savemodel_called;
  mutable bool getvtf_called;
  mutable bool getprocessadm_called;
  bool registerdependency_called;
  bool registerdependency2_called;
  bool fillfield_called;
  bool getfield_called;
  mutable bool getfield2_called;
  mutable bool getdependentfield_called;
  bool setupdependencies_called;
  bool init_called;
  bool registerfields_called;
  bool savestep_called;
  bool advancestep_called;
  bool postsolve_called;
  mutable bool getglobalnode_called;
  mutable bool getlocalnode_called;
};

TEST(TestSIMOverride, Override)
{
  SIMMockOverride ovr;
  SIMOverride<SIMMockOverride> sim(ovr);
  sim.read("");
  sim.setInitialConditions();
  sim.hasIC("");
  sim.getNoSpaceDim();
  sim.getFEModel();
  int geoBlk, nBlock;
  sim.saveModel((char*)"", geoBlk, nBlock);
  sim.getVTF();
  sim.getProcessAdm();
  sim.registerDependency(NULL, "", 1, std::vector<ASMbase*>());
  sim.registerDependency(NULL, "", 1);
  sim.fillField("", std::vector<double>());
  sim.getField("");
  const_cast<const SIMOverride<SIMMockOverride>&>(sim).getField("");
  sim.getDependentField("");
  sim.setupDependencies();
  TimeStep tp;
  sim.init(tp);
  DataExporter exp;
  sim.registerFields(exp);
  sim.saveStep(tp, nBlock);
  sim.advanceStep(tp);
  sim.postSolve(tp);
  sim.getGlobalNode(1);
  sim.getLocalNode(1);

  ASSERT_TRUE(ovr.read_called);
  ASSERT_TRUE(ovr.setinitialconditions_called);
  ASSERT_TRUE(ovr.hasic_called);
  ASSERT_TRUE(ovr.getnospacedim_called);
  ASSERT_TRUE(ovr.getfemodel_called);
  ASSERT_TRUE(ovr.savemodel_called);
  ASSERT_TRUE(ovr.getvtf_called);
  ASSERT_TRUE(ovr.getprocessadm_called);
  ASSERT_TRUE(ovr.registerdependency_called);
  ASSERT_TRUE(ovr.registerdependency2_called);
  ASSERT_TRUE(ovr.fillfield_called);
  ASSERT_TRUE(ovr.getfield_called);
  ASSERT_TRUE(ovr.getfield2_called);
  ASSERT_TRUE(ovr.getdependentfield_called);
  ASSERT_TRUE(ovr.setupdependencies_called);
  ASSERT_TRUE(ovr.init_called);
  ASSERT_TRUE(ovr.registerfields_called);
  ASSERT_TRUE(ovr.savestep_called);
  ASSERT_TRUE(ovr.advancestep_called);
  ASSERT_TRUE(ovr.postsolve_called);
  ASSERT_TRUE(ovr.getglobalnode_called);
  ASSERT_TRUE(ovr.getlocalnode_called);
}
