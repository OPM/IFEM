//==============================================================================
//!
//! \file TestSIMCoupled.C
//!
//! \date Oct 13 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for template overriding methods for wrapping template purposes.
//!
//==============================================================================

#include "DataExporter.h"
#include "Function.h"
#include "ProcessAdm.h"
#include "Property.h"
#include "TimeStep.h"
#include "SIMdependency.h"
#include "SIMCoupled.h"

#include "gtest/gtest.h"

class VTF;

template<class T1, class T2>
class SIMMockCoupling : public SIMCoupled<T1,T2>
{
  public:
    SIMMockCoupling<T1,T2>(T1& t1, T2& t2) : SIMCoupled<T1, T2>(t1, t2) {}
    void setupDependencies() {}
};

class SIMMockCoupled {
public:
  SIMMockCoupled() :
    preprocess_called(false),
    advancestep_called(false),
    solvestep_called(false),
    postsolve_called(false),
    savestep_called(false),
    savemodel_called(false),
    init_called(false),
    registerdependency_called(false),
    registerdependency2_called(false),
    getuniquepropertycode_called(false),
    createpropertyset_called(false),
    setvecproperty_called(false),
    registerfields_called(false),
    setinitialconditions_called(false),
    hasic_called(false),
    getfield_called(false),
    setvtf_called(false),
    getvtf_called(false)
  {}

  bool preprocess() { preprocess_called = true; return true; }
  bool advanceStep(TimeStep& tp) { advancestep_called = true; return true; }
  bool solveStep(TimeStep& tp) { solvestep_called = true; return true; }
  void postSolve(const TimeStep&tp, bool restart=false) { postsolve_called = true; }
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    savestep_called = true;
    return true;
  }
  bool saveModel(char* filename, int& geoBlk, int& nBlock)
  {
    savemodel_called = true;
    return true;
  }
  bool init(const TimeStep& tp) { init_called = true; return true; }
  void registerDependency(SIMdependency* sim, const std::string& name,
                          short int nvc,
                          const SIMdependency::PatchVec& patches,
                          bool diffBasis = false)
  {
    registerdependency_called = true;
  }
  void registerDependency(SIMdependency* sim, const std::string& name,
                          short int nvc = 1)
  {
    registerdependency2_called = true;
  }
  int getUniquePropertyCode(const std::string& setName, int comp=0)
  {
    getuniquepropertycode_called = true;
    return 1;
  }
  bool createPropertySet(const std::string& setName, int pc)
  {
    createpropertyset_called = true;
    return true;
  }
  size_t setVecProperty(int code, Property::Type ptype,
                        VecFunc* field=NULL, int pflag=-1)
  {
    setvecproperty_called = true;
    return 1;
  }
  void registerFields(DataExporter& exporter) { registerfields_called = true; }
  void setInitialConditions() { setinitialconditions_called = true; }
  bool hasIC(const std::string& name) const
  {
    hasic_called = true;
    return false;
  }
  utl::vector<double>* getField(const std::string& name)
  {
    getfield_called = true;
    return NULL;
  }
  void setVTF(VTF* vtf) { setvtf_called = true; }
  VTF* getVTF() { getvtf_called = true; return NULL; }

  bool preprocess_called;
  bool advancestep_called;
  bool solvestep_called;
  bool postsolve_called;
  bool savestep_called;
  bool savemodel_called;
  bool init_called;
  bool registerdependency_called;
  bool registerdependency2_called;
  bool getuniquepropertycode_called;
  bool createpropertyset_called;
  bool setvecproperty_called;
  bool registerfields_called;
  bool setinitialconditions_called;
  mutable bool hasic_called;
  bool getfield_called;
  bool setvtf_called;
  bool getvtf_called;
};

TEST(TestSIMCoupled, Override)
{
  SIMMockCoupled ovr1, ovr2;
  SIMMockCoupling<SIMMockCoupled, SIMMockCoupled> sim(ovr1, ovr2);
  TimeStep tp;
  int geoBlk, nBlock;
  DataExporter exporter;
  sim.preprocess();
  sim.advanceStep(tp);
  sim.solveStep(tp);
  sim.postSolve(tp);
  sim.saveStep(tp, nBlock);
  sim.saveModel((char*)"", geoBlk, nBlock);
  sim.init(tp);
  sim.registerDependency(NULL, "", 2, SIMdependency::PatchVec());
  sim.registerDependency(NULL, "", 2);
  sim.getUniquePropertyCode("", 0);
  sim.createPropertySet("", 0);
  sim.setVecProperty(0, Property::MATERIAL, NULL, -1);
  sim.registerFields(exporter);
  sim.setInitialConditions();
  sim.hasIC("");
  sim.getField("");

  ASSERT_TRUE(ovr1.preprocess_called);
  ASSERT_TRUE(ovr2.preprocess_called);
  ASSERT_TRUE(ovr1.advancestep_called);
  ASSERT_TRUE(ovr2.advancestep_called);
  ASSERT_TRUE(ovr1.solvestep_called);
  ASSERT_TRUE(ovr2.solvestep_called);
  ASSERT_TRUE(ovr1.postsolve_called);
  ASSERT_TRUE(ovr2.postsolve_called);
  ASSERT_TRUE(ovr1.savemodel_called);
  ASSERT_TRUE(ovr1.init_called);
  ASSERT_TRUE(ovr2.init_called);
  ASSERT_TRUE(ovr1.registerdependency_called);
  ASSERT_TRUE(ovr2.registerdependency_called);
  ASSERT_TRUE(ovr1.registerdependency2_called);
  ASSERT_TRUE(ovr2.registerdependency2_called);
  ASSERT_TRUE(ovr1.getuniquepropertycode_called);
  ASSERT_FALSE(ovr2.getuniquepropertycode_called);
  ASSERT_TRUE(ovr1.createpropertyset_called);
  ASSERT_FALSE(ovr2.createpropertyset_called);
  ASSERT_TRUE(ovr1.setvecproperty_called);
  ASSERT_FALSE(ovr2.setvecproperty_called);
  ASSERT_TRUE(ovr1.registerfields_called);
  ASSERT_TRUE(ovr2.registerfields_called);
  ASSERT_TRUE(ovr1.setinitialconditions_called);
  ASSERT_TRUE(ovr2.setinitialconditions_called);
  ASSERT_TRUE(ovr1.hasic_called);
  ASSERT_TRUE(ovr2.hasic_called);
  ASSERT_TRUE(ovr1.getfield_called);
  ASSERT_TRUE(ovr2.getfield_called);
  ASSERT_TRUE(ovr1.getvtf_called);
  ASSERT_FALSE(ovr2.getvtf_called);
  ASSERT_FALSE(ovr1.setvtf_called);
  ASSERT_TRUE(ovr2.setvtf_called);
}
