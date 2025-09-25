//==============================================================================
//!
//! \file TestMultiPatchLRRefine.C
//!
//! \date Mar 31 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for multi-patch LR refinement.
//!
//==============================================================================

#include "ASMbase.h"
#include "ASMunstruct.h"
#include "IntegrandBase.h"
#include "MultiPatchModelGenerator.h"
#include "SIMMultiPatchModelGen.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "tinyxml2.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <sstream>

namespace {

// Dummy SIM class.
template<class Dim>
class RefineSim : public SIMMultiPatchModelGen<Dim>
{
public:
  RefineSim() : SIMMultiPatchModelGen<Dim>(Dim::dimension)
  {
    Dim::opt.discretization = ASM::LRSpline;
  }

  bool parse(const tinyxml2::XMLElement* elem) override
  { return this->SIMMultiPatchModelGen<Dim>::parse(elem); }

  virtual ~RefineSim() {}
};

}


TEST_CASE("TestMultiPatchLRRefine.2D")
{
  const int param = GENERATE(0,1,2);

  SECTION(std::to_string(param) + "refinements") {
    RefineSim<SIM2D> sim;
    tinyxml2::XMLDocument doc;
    std::stringstream str;
    str << R"(<geometry dim="2" nx="2" ny="2">)"
        << R"(  <raiseorder lowerpatch="1" upperpatch="4)"
        << R"(" u=")" << param
        << R"(" v=")" << param << '"' << "/>"
        << "</geometry>";
    doc.Parse(str.str().c_str());
    MultiPatchModelGenerator2D gen(doc.RootElement());
    REQUIRE(sim.parse(doc.RootElement()));
    REQUIRE(sim.preprocess());

    srand(0);

    for (size_t i = 0; i < 4; ++i) {
      LR::RefineData prm;
      sim.getPatch(1 + (rand() % 4))->getBoundaryNodes(1 + (rand() % 4), prm.elements);

      prm.options.resize(3);
      prm.options[0] = 1;
      prm.options[1] = 1;
      prm.options[2] = 2;
      sim.refine(prm);
      sim.clearProperties();

      REQUIRE(gen.createTopology(sim));
      REQUIRE(sim.preprocess());
    }
  }
}


TEST_CASE("TestMultiPatchLRRefine.3D")
{
  const int param = GENERATE(0,1,2);

  SECTION(std::to_string(param) + "refinements") {
    RefineSim<SIM3D> sim;
    tinyxml2::XMLDocument doc;
    std::stringstream str;
    str << R"(<geometry dim="3" nx="2" ny="2" nz="2">)"
        << R"(  <raiseorder lowerpatch="1" upperpatch="8")"
        << R"( u=")" << param
        << R"(" v=")" << param
        << R"(" w=")" << param
        << '"' << "/>"
        << "</geometry>";
    doc.Parse(str.str().c_str());
    MultiPatchModelGenerator3D gen(doc.RootElement());
    REQUIRE(sim.parse(doc.RootElement()));
    REQUIRE(sim.preprocess());

    srand(0);

    for (size_t i = 0; i < 3; ++i) {
      LR::RefineData prm;
      sim.getPatch(1 + (rand() % 8))->getBoundaryNodes(1 + (rand() % 6), prm.elements);

      prm.options.resize(3);
      prm.options[0] = 1;
      prm.options[1] = 1;
      prm.options[2] = 2;
      sim.refine(prm);
      sim.clearProperties();

      REQUIRE(gen.createTopology(sim));
      REQUIRE(sim.preprocess());
    }
  }
}
