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
#include "Property.h"
#include "TimeStep.h"
#include "matrix.h"
#include "SIMCoupled.h"

#include "Catch2Support.h"


namespace {

class SIMMockCoupled {
public:
  SIMMockCoupled() :
    preprocess_called(false),
    advancestep_called(false),
    solvestep_called(false),
    savestep_called(false),
    savemodel_called(false),
    init_called(false),
    registerdependency1_called(false),
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

  bool preprocess() { return preprocess_called = true; }
  bool advanceStep(TimeStep& tp) { return advancestep_called = true; }
  bool solveStep(TimeStep& tp) { return solvestep_called = true; }
  bool saveStep(const TimeStep&, int&) { return savestep_called = true; }
  bool saveModel(char*, int&, int&) { return savemodel_called = true; }
  bool init(const TimeStep&) { return init_called = true; }
  void registerDependency(SIMdependency*, const std::string&,
                          short int, const std::vector<ASMbase*>&, char, int)
  {
    registerdependency1_called = true;
  }
  void registerDependency(SIMdependency*, const std::string&, short int)
  {
    registerdependency2_called = true;
  }
  int getUniquePropertyCode(const std::string&, int)
  {
    getuniquepropertycode_called = true;
    return 1;
  }
  bool createPropertySet(const std::string&, int)
  {
    return createpropertyset_called = true;
  }
  size_t setVecProperty(int, Property::Type, VecFunc*, int)
  {
    setvecproperty_called = true;
    return 1;
  }
  void registerFields(DataExporter&) { registerfields_called = true; }
  bool setInitialConditions() { return setinitialconditions_called = true; }
  bool hasIC(const std::string&) const { hasic_called = true; return false; }
  utl::vector<double>* getField(const std::string& name)
  {
    getfield_called = true;
    return nullptr;
  }
  void setVTF(VTF* vtf) { setvtf_called = true; }
  VTF* getVTF() { getvtf_called = true; return nullptr; }

  bool preprocess_called;
  bool advancestep_called;
  bool solvestep_called;
  bool savestep_called;
  bool savemodel_called;
  bool init_called;
  bool registerdependency1_called;
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

}

TEST_CASE("TestSIMCoupled.Override")
{
  SIMMockCoupled ovr1, ovr2;
  SIMCoupled<SIMMockCoupled,SIMMockCoupled> sim(ovr1,ovr2);
  TimeStep tp;
  int geoBlk, nBlock;
  DataExporter exporter;

  sim.preprocess();
  sim.advanceStep(tp);
  sim.solveStep(tp);
  sim.saveStep(tp, nBlock);
  sim.saveModel(const_cast<char*>(""), geoBlk, nBlock);
  sim.init(tp);
  sim.registerDependency(nullptr, "", 2, std::vector<ASMbase*>());
  sim.registerDependency(nullptr, "", 2);
  sim.getUniquePropertyCode("", 0);
  sim.createPropertySet("", 0);
  sim.setVecProperty(0, Property::MATERIAL, nullptr, -1);
  sim.registerFields(exporter);
  sim.setInitialConditions();
  sim.hasIC("");
  sim.getField("");
  sim.getVTF();

  REQUIRE(ovr1.preprocess_called);
  REQUIRE(ovr2.preprocess_called);
  REQUIRE(ovr1.advancestep_called);
  REQUIRE(ovr2.advancestep_called);
  REQUIRE(ovr1.solvestep_called);
  REQUIRE(ovr2.solvestep_called);
  REQUIRE(ovr1.savemodel_called);
  REQUIRE(ovr1.init_called);
  REQUIRE(ovr2.init_called);
  REQUIRE(ovr1.registerdependency1_called);
  REQUIRE(ovr2.registerdependency1_called);
  REQUIRE(ovr1.registerdependency2_called);
  REQUIRE(ovr2.registerdependency2_called);
  REQUIRE(ovr1.getuniquepropertycode_called);
  REQUIRE(!ovr2.getuniquepropertycode_called);
  REQUIRE(ovr1.createpropertyset_called);
  REQUIRE(!ovr2.createpropertyset_called);
  REQUIRE(ovr1.setvecproperty_called);
  REQUIRE(!ovr2.setvecproperty_called);
  REQUIRE(ovr1.registerfields_called);
  REQUIRE(ovr2.registerfields_called);
  REQUIRE(ovr1.setinitialconditions_called);
  REQUIRE(ovr2.setinitialconditions_called);
  REQUIRE(ovr1.hasic_called);
  REQUIRE(ovr2.hasic_called);
  REQUIRE(ovr1.getfield_called);
  REQUIRE(ovr2.getfield_called);
  REQUIRE(ovr1.getvtf_called);
  REQUIRE(!ovr2.getvtf_called);
  REQUIRE(!ovr1.setvtf_called);
  REQUIRE(ovr2.setvtf_called);
}
