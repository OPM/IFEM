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
#include "tinyxml.h"

#include "gtest/gtest.h"

#include <fstream>


// Dummy SIM class.
template<class Dim>
class RefineSim : public SIMMultiPatchModelGen<Dim>
{
public:
  RefineSim() : SIMMultiPatchModelGen<Dim>(Dim::dimension)
  {
    Dim::opt.discretization = ASM::LRSpline;
  }

  bool parse(const TiXmlElement* elem) override
  { return this->SIMMultiPatchModelGen<Dim>::parse(elem); }

  virtual ~RefineSim() {}
};


class TestMultiPatchLRRefine2D :
  public testing::Test,
  public testing::WithParamInterface<int>
{
};


TEST_P(TestMultiPatchLRRefine2D, Refine)
{
  RefineSim<SIM2D> sim;
  TiXmlDocument doc;
  std::stringstream str;
  str << R"(<geometry dim="2" nx="2" ny="2">)"
      << R"(  <raiseorder lowerpatch="1" upperpatch="4)"
      << R"(" u=")" << GetParam()
      << R"(" v=")" << GetParam() << '"' << "/>"
      << "</geometry>";
  doc.Parse(str.str().c_str());
  MultiPatchModelGenerator2D gen(doc.RootElement());
  EXPECT_TRUE(sim.parse(doc.RootElement()));
  EXPECT_TRUE(sim.preprocess());

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

    ASSERT_TRUE(gen.createTopology(sim) && sim.preprocess());
  }
}

class TestMultiPatchLRRefine3D :
  public testing::Test,
  public testing::WithParamInterface<int>
{
};


TEST_P(TestMultiPatchLRRefine3D, Refine)
{
  RefineSim<SIM3D> sim;
  TiXmlDocument doc;
  std::stringstream str;
  str << R"(<geometry dim="3" nx="2" ny="2" nz="2">)"
      << R"(  <raiseorder lowerpatch="1" upperpatch="8")"
      << R"( u=")" << GetParam()
      << R"(" v=")" << GetParam()
      << R"(" w=")" << GetParam()
      << '"' << "/>"
      << "</geometry>";
  doc.Parse(str.str().c_str());
  MultiPatchModelGenerator3D gen(doc.RootElement());
  EXPECT_TRUE(sim.parse(doc.RootElement()));
  EXPECT_TRUE(sim.preprocess());

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

    ASSERT_TRUE(gen.createTopology(sim) && sim.preprocess());
  }
}


const std::vector<int> refValues = {0,1,2};
INSTANTIATE_TEST_CASE_P(TestMultiPatchLRRefine2D,
                        TestMultiPatchLRRefine2D,
                        testing::ValuesIn(refValues));
INSTANTIATE_TEST_CASE_P(TestMultiPatchLRRefine3D,
                        TestMultiPatchLRRefine3D,
                        testing::ValuesIn(refValues));
