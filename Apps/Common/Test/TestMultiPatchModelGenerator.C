//==============================================================================
//!
//! \file TestMultiPatchModelGenerator.C
//!
//! \date Sep 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for multi-patch model generators.
//!
//==============================================================================

#include "MultiPatchModelGenerator.h"
#include "SIMMultiPatchModelGen.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "Functions.h"
#include "Utilities.h"

#include "tinyxml2.h"

#include <GoTools/geometry/SplineCurve.h>
#include <GoTools/geometry/SplineSurface.h>
#include <GoTools/utils/Point.h>

#include "Catch2Support.h"


namespace {

template<class Generator>
class TestModelGeneratorWrapper : public Generator {
public:
  TestModelGeneratorWrapper(const tinyxml2::XMLElement* geo) : Generator(geo) {}
  virtual std::string createG2(int nsd, bool = false) const
  {
    bool rational = false;
    utl::getAttribute(Generator::geo,"rational",rational);
    return this->Generator::createG2(nsd,rational);
  }
};

struct GeomTest {
  const char* name;
  const char* xml;
  int dim;
  std::string g2;
  std::string sets;
};


void DoTest(const GeomTest& ref, const std::string& gen,
            const TopologySet& sets)
{
  REQUIRE(gen == ref.g2);

  if (ref.sets.empty())
    return;

  std::string gsets;
  for (const TopologySet::value_type& tset : sets) {
    gsets += tset.first + ":";
    for (const TopItem& item : tset.second) {
      std::stringstream str;
      str <<" "<< item.patch <<" "<< item.item <<" "<< item.idim;
      gsets += str.str();
    }
    gsets += " \n";
  }
  REQUIRE(gsets == ref.sets);
};

}


TEST_CASE("TestMultiPatchModelGenerator.Generate2D")
{
  const GeomTest param = GENERATE(
    GeomTest{"Sets",
            "<geometry sets=\"true\"/>", 2,
             "200 1 0 0\n"
             "2 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0\n"
             "1 0\n"
             "0 1\n"
             "1 1\n",
             "Boundary: 1 1 1 1 2 1 1 3 1 1 4 1 \n"
             "BoundaryX: 1 1 1 1 2 1 \n"
             "BoundaryY: 1 3 1 1 4 1 \n"
             "Corners: 1 1 0 1 2 0 1 3 0 1 4 0 \n"
             "Edge1: 1 1 1 \n"
             "Edge1Patches: 1 0 2 \n"
             "Edge2: 1 2 1 \n"
             "Edge2Patches: 1 0 2 \n"
             "Edge3: 1 3 1 \n"
             "Edge3Patches: 1 0 2 \n"
             "Edge4: 1 4 1 \n"
             "Edge4Patches: 1 0 2 \n"
             "InnerPatches: \n"
             "Vertex1: 1 1 0 \n"
             "Vertex2: 1 2 0 \n"
             "Vertex3: 1 3 0 \n"
             "Vertex4: 1 4 0 \n"},

    GeomTest{"Rational",
             "<geometry rational=\"1\"/>", 2,
             "200 1 0 0\n"
             "2 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0 1.0\n"
             "1 0 1.0\n"
             "0 1 1.0\n"
             "1 1 1.0\n", ""},

    GeomTest{"Scale",
             "<geometry scale=\"2\"/>", 2,
             "200 1 0 0\n"
             "2 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0\n"
             "2 0\n"
             "0 2\n"
             "2 2\n", ""},

    GeomTest{"Corner X",
             "<geometry X0=\"2 0\"/>", 2,
             "200 1 0 0\n"
             "2 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 0\n"
             "3 0\n"
             "2 1\n"
             "3 1\n"},

    GeomTest{"Corner Y",
             "<geometry X0=\"0 2\"/>", 2,
             "200 1 0 0\n"
             "2 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 2\n"
             "1 2\n"
             "0 3\n"
             "1 3\n", ""},

    GeomTest{"Rectangle X",
             "<geometry Lx=\"2\"/>", 2,
             "200 1 0 0\n"
             "2 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0\n"
             "2 0\n"
             "0 1\n"
             "2 1\n", ""},

    GeomTest{"Rectangle Y",
             "<geometry Ly=\"2\"/>", 2,
             "200 1 0 0\n"
             "2 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0\n"
             "1 0\n"
             "0 2\n"
             "1 2\n", ""},

    GeomTest{"Multiple X",
             "<geometry sets=\"true\" nx=\"2\"/>", 2,
             "200 1 0 0\n"
             "2 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0\n"
             "0.5 0\n"
             "0 1\n"
             "0.5 1\n"
             "200 1 0 0\n"
             "2 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0.5 0\n"
             "1 0\n"
             "0.5 1\n"
             "1 1\n",
             "Boundary: 1 1 1 1 3 1 1 4 1 2 2 1 2 3 1 2 4 1 \n"
             "BoundaryX: 1 1 1 2 2 1 \n"
             "BoundaryY: 1 3 1 1 4 1 2 3 1 2 4 1 \n"
             "Corners: 1 1 0 1 3 0 2 2 0 2 4 0 \n"
             "Edge1: 1 1 1 \n"
             "Edge1Patches: 1 0 2 \n"
             "Edge2: 2 2 1 \n"
             "Edge2Patches: 2 0 2 \n"
             "Edge3: 1 3 1 2 3 1 \n"
             "Edge3Patches: 1 0 2 2 0 2 \n"
             "Edge4: 1 4 1 2 4 1 \n"
             "Edge4Patches: 1 0 2 2 0 2 \n"
             "InnerPatches: \n"
             "Vertex1: 1 1 0 \n"
             "Vertex2: 2 2 0 \n"
             "Vertex3: 1 3 0 \n"
             "Vertex4: 2 4 0 \n"},

    GeomTest{"Multiple Y",
             "<geometry sets=\"true\" ny=\"2\"/>", 2,
             "200 1 0 0\n"
             "2 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0\n"
             "1 0\n"
             "0 0.5\n"
             "1 0.5\n"
             "200 1 0 0\n"
             "2 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0.5\n"
             "1 0.5\n"
             "0 1\n"
             "1 1\n",
             "Boundary: 1 1 1 1 2 1 1 3 1 2 1 1 2 2 1 2 4 1 \n"
             "BoundaryX: 1 1 1 1 2 1 2 1 1 2 2 1 \n"
             "BoundaryY: 1 3 1 2 4 1 \n"
             "Corners: 1 1 0 1 2 0 2 3 0 2 4 0 \n"
             "Edge1: 1 1 1 2 1 1 \n"
             "Edge1Patches: 1 0 2 2 0 2 \n"
             "Edge2: 1 2 1 2 2 1 \n"
             "Edge2Patches: 1 0 2 2 0 2 \n"
             "Edge3: 1 3 1 \n"
             "Edge3Patches: 1 0 2 \n"
             "Edge4: 2 4 1 \n"
             "Edge4Patches: 2 0 2 \n"
             "InnerPatches: \n"
             "Vertex1: 1 1 0 \n"
             "Vertex2: 1 2 0 \n"
             "Vertex3: 2 3 0 \n"
             "Vertex4: 2 4 0 \n"},

    GeomTest{"Multiple",
             "<geometry sets=\"true\" nx=\"2\" ny=\"2\"/>", 2,
             "200 1 0 0\n"
             "2 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0\n"
             "0.5 0\n"
             "0 0.5\n"
             "0.5 0.5\n"
             "200 1 0 0\n"
             "2 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0.5 0\n"
             "1 0\n"
             "0.5 0.5\n"
             "1 0.5\n"
             "200 1 0 0\n"
             "2 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0.5\n"
             "0.5 0.5\n"
             "0 1\n"
             "0.5 1\n"
             "200 1 0 0\n"
             "2 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0.5 0.5\n"
             "1 0.5\n"
             "0.5 1\n"
             "1 1\n",
             "Boundary: 1 1 1 1 3 1 2 2 1 2 3 1 3 1 1 3 4 1 4 2 1 4 4 1 \n"
             "BoundaryX: 1 1 1 2 2 1 3 1 1 4 2 1 \n"
             "BoundaryY: 1 3 1 2 3 1 3 4 1 4 4 1 \n"
             "Corners: 1 1 0 2 2 0 3 3 0 4 4 0 \n"
             "Edge1: 1 1 1 3 1 1 \n"
             "Edge1Patches: 1 0 2 3 0 2 \n"
             "Edge2: 2 2 1 4 2 1 \n"
             "Edge2Patches: 2 0 2 4 0 2 \n"
             "Edge3: 1 3 1 2 3 1 \n"
             "Edge3Patches: 1 0 2 2 0 2 \n"
             "Edge4: 3 4 1 4 4 1 \n"
             "Edge4Patches: 3 0 2 4 0 2 \n"
             "InnerPatches: \n"
             "Vertex1: 1 1 0 \n"
             "Vertex2: 2 2 0 \n"
             "Vertex3: 3 3 0 \n"
             "Vertex4: 4 4 0 \n"}
  );

  SECTION(param.name) {
    tinyxml2::XMLDocument doc;
    doc.Parse(param.xml);
    TestModelGeneratorWrapper<MultiPatchModelGenerator2D> gen(doc.RootElement());
    SIM2D sim;
    gen.createTopologySets(sim);
    DoTest(param, gen.createG2(param.dim), sim.getTopology());
  }
}


class SIMMultiPatch2D : public SIMMultiPatchModelGen<SIM2D>
{
public:
  explicit SIMMultiPatch2D(int n1) : SIMMultiPatchModelGen<SIM2D>(n1,false) {}
  SIMMultiPatch2D(unsigned char n1, unsigned char n2)
    : SIMMultiPatchModelGen<SIM2D>({n1,n2},false) {}
};


class SIMMultiPatch3D : public SIMMultiPatchModelGen<SIM3D>
{
public:
  explicit SIMMultiPatch3D(int n1) : SIMMultiPatchModelGen<SIM3D>(n1,false) {}
};


TEST_CASE("TestMultiPatchModelGenerator.Generate2DLR")
{
  SIMMultiPatch2D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  REQUIRE(sim.read("refdata/modelgen2d_lr.xinp"));
  sim.preprocess();
  REQUIRE(sim.getNoNodes() == 28);
}


TEST_CASE("TestMultiPatchModelGeneratorGenerate2DLRmx")
{
  SIMMultiPatch2D sim(2,1);
  sim.opt.discretization = ASM::LRSpline;
  REQUIRE(sim.read("refdata/modelgen2d_lr.xinp"));
  sim.preprocess();
  REQUIRE(sim.getNoNodes() == 73);
}


TEST_CASE("TestMultiPatchModelGenerator.InnerPatches2D")
{
  tinyxml2::XMLDocument doc;
  doc.Parse("<geometry nx=\"3\" ny=\"3\" sets=\"true\"/>");
  TestModelGeneratorWrapper<MultiPatchModelGenerator2D> gen(doc.RootElement());
  SIM2D sim;
  gen.createTopologySets(sim);
  REQUIRE(sim.topology("InnerPatches").size() == 1);
  REQUIRE(sim.topology("InnerPatches").begin()->patch == 5);
}


TEST_CASE("TestMultiPatchModelGenerator.Generate3D")
{
  const GeomTest param = GENERATE(
    GeomTest{"Sets",
             "<geometry sets=\"true\"/>", 3,
             "700 1 0 0\n"
             "3 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0 0\n"
             "1 0 0\n"
             "0 1 0\n"
             "1 1 0\n"
             "0 0 1\n"
             "1 0 1\n"
             "0 1 1\n"
             "1 1 1\n",
             "Boundary: 1 1 2 1 2 2 1 3 2 1 4 2 1 5 2 1 6 2 \n"
             "BoundaryX: 1 1 2 1 2 2 \n"
             "BoundaryY: 1 3 2 1 4 2 \n"
             "BoundaryZ: 1 5 2 1 6 2 \n"
             "Corners: 1 1 0 1 2 0 1 3 0 1 4 0 1 5 0 1 6 0 1 7 0 1 8 0 \n"
             "Edge1: 1 1 1 \n"
             "Edge10: 1 10 1 \n"
             "Edge11: 1 11 1 \n"
             "Edge12: 1 12 1 \n"
             "Edge2: 1 2 1 \n"
             "Edge3: 1 3 1 \n"
             "Edge4: 1 4 1 \n"
             "Edge5: 1 5 1 \n"
             "Edge6: 1 6 1 \n"
             "Edge7: 1 7 1 \n"
             "Edge8: 1 8 1 \n"
             "Edge9: 1 9 1 \n"
             "Face1: 1 1 2 \n"
             "Face1Patches: 1 0 3 \n"
             "Face2: 1 2 2 \n"
             "Face2Patches: 1 0 3 \n"
             "Face3: 1 3 2 \n"
             "Face3Patches: 1 0 3 \n"
             "Face4: 1 4 2 \n"
             "Face4Patches: 1 0 3 \n"
             "Face5: 1 5 2 \n"
             "Face5Patches: 1 0 3 \n"
             "Face6: 1 6 2 \n"
             "Face6Patches: 1 0 3 \n"
             "Frame: 1 1 1 1 2 1 1 3 1 1 4 1 1 5 1 1 6 1 1 7 1 1 8 1 1 9 1 1 10 1 1 11 1 1 12 1 \n"
             "InnerPatches: \n"
             "Vertex1: 1 1 0 \n"
             "Vertex2: 1 2 0 \n"
             "Vertex3: 1 3 0 \n"
             "Vertex4: 1 4 0 \n"
             "Vertex5: 1 5 0 \n"
             "Vertex6: 1 6 0 \n"
             "Vertex7: 1 7 0 \n"
             "Vertex8: 1 8 0 \n"},

    GeomTest{"Rational",
             "<geometry rational=\"1\"/>", 3,
             "700 1 0 0\n"
             "3 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0 0 1.0\n"
             "1 0 0 1.0\n"
             "0 1 0 1.0\n"
             "1 1 0 1.0\n"
             "0 0 1 1.0\n"
             "1 0 1 1.0\n"
             "0 1 1 1.0\n"
             "1 1 1 1.0\n", ""},

    GeomTest{"Scale",
             "<geometry scale=\"2\"/>", 3,
             "700 1 0 0\n"
             "3 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0 0\n"
             "2 0 0\n"
             "0 2 0\n"
             "2 2 0\n"
             "0 0 2\n"
             "2 0 2\n"
             "0 2 2\n"
             "2 2 2\n", ""},

    GeomTest{"Corner X",
             "<geometry X0=\"2 0 0\"/>", 3,
             "700 1 0 0\n"
             "3 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 0 0\n"
             "3 0 0\n"
             "2 1 0\n"
             "3 1 0\n"
             "2 0 1\n"
             "3 0 1\n"
             "2 1 1\n"
             "3 1 1\n", ""},

    GeomTest{"Corner Y",
             "<geometry X0=\"0 2 0\"/>", 3,
             "700 1 0 0\n"
             "3 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 2 0\n"
             "1 2 0\n"
             "0 3 0\n"
             "1 3 0\n"
             "0 2 1\n"
             "1 2 1\n"
             "0 3 1\n"
             "1 3 1\n", ""},

    GeomTest{"Corner Z",
             "<geometry X0=\"0 0 2\"/>", 3,
             "700 1 0 0\n"
             "3 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0 2\n"
             "1 0 2\n"
             "0 1 2\n"
             "1 1 2\n"
             "0 0 3\n"
             "1 0 3\n"
             "0 1 3\n"
             "1 1 3\n", ""},

    GeomTest{"Cuboid X",
             "<geometry Lx=\"2\"/>", 3,
             "700 1 0 0\n"
             "3 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0 0\n"
             "2 0 0\n"
             "0 1 0\n"
             "2 1 0\n"
             "0 0 1\n"
             "2 0 1\n"
             "0 1 1\n"
             "2 1 1\n", ""},

    GeomTest{"Cuboid Y",
             "<geometry Ly=\"2\"/>", 3,
             "700 1 0 0\n"
             "3 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0 0\n"
             "1 0 0\n"
             "0 2 0\n"
             "1 2 0\n"
             "0 0 1\n"
             "1 0 1\n"
             "0 2 1\n"
             "1 2 1\n", ""},

    GeomTest{"Cuboid Z",
             "<geometry Lz=\"2\"/>", 3,
             "700 1 0 0\n"
             "3 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0 0\n"
             "1 0 0\n"
             "0 1 0\n"
             "1 1 0\n"
             "0 0 2\n"
             "1 0 2\n"
             "0 1 2\n"
             "1 1 2\n", ""},

    GeomTest{"Multiple X",
             "<geometry sets=\"true\" nx=\"2\"/>", 3,
             "700 1 0 0\n"
             "3 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0 0\n"
             "0.5 0 0\n"
             "0 1 0\n"
             "0.5 1 0\n"
             "0 0 1\n"
             "0.5 0 1\n"
             "0 1 1\n"
             "0.5 1 1\n"
             "700 1 0 0\n"
             "3 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0.5 0 0\n"
             "1 0 0\n"
             "0.5 1 0\n"
             "1 1 0\n"
             "0.5 0 1\n"
             "1 0 1\n"
             "0.5 1 1\n"
             "1 1 1\n",
             "Boundary: 1 1 2 1 3 2 1 4 2 1 5 2 1 6 2 2 2 2 2 3 2 2 4 2 2 5 2 2 6 2 \n"
             "BoundaryX: 1 1 2 2 2 2 \n"
             "BoundaryY: 1 3 2 1 4 2 2 3 2 2 4 2 \n"
             "BoundaryZ: 1 5 2 1 6 2 2 5 2 2 6 2 \n"
             "Corners: 1 1 0 1 3 0 1 5 0 1 7 0 2 2 0 2 4 0 2 6 0 2 8 0 \n"
             "Edge1: 1 1 1 \n"
             "Edge10: 1 10 1 2 10 1 \n"
             "Edge11: 1 11 1 2 11 1 \n"
             "Edge12: 1 12 1 2 12 1 \n"
             "Edge2: 2 2 1 \n"
             "Edge3: 1 3 1 \n"
             "Edge4: 2 4 1 \n"
             "Edge5: 1 5 1 \n"
             "Edge6: 2 6 1 \n"
             "Edge7: 1 7 1 \n"
             "Edge8: 2 8 1 \n"
             "Edge9: 1 9 1 2 9 1 \n"
             "Face1: 1 1 2 \n"
             "Face1Patches: 1 0 3 \n"
             "Face2: 2 2 2 \n"
             "Face2Patches: 2 0 3 \n"
             "Face3: 1 3 2 2 3 2 \n"
             "Face3Patches: 1 0 3 2 0 3 \n"
             "Face4: 1 4 2 2 4 2 \n"
             "Face4Patches: 1 0 3 2 0 3 \n"
             "Face5: 1 5 2 2 5 2 \n"
             "Face5Patches: 1 0 3 2 0 3 \n"
             "Face6: 1 6 2 2 6 2 \n"
             "Face6Patches: 1 0 3 2 0 3 \n"
             "Frame: 1 1 1 1 3 1 1 5 1 1 7 1 1 9 1 1 10 1 1 11 1 1 12 1 2 2 1 2 4 1 2 6 1 2 8 1 2 9 1 2 10 1 2 11 1 2 12 1 \n"
             "InnerPatches: \n"
             "Vertex1: 1 1 0 \n"
             "Vertex2: 2 2 0 \n"
             "Vertex3: 1 3 0 \n"
             "Vertex4: 2 4 0 \n"
             "Vertex5: 1 5 0 \n"
             "Vertex6: 2 6 0 \n"
             "Vertex7: 1 7 0 \n"
             "Vertex8: 2 8 0 \n"},

    GeomTest{"Multiple Y",
             "<geometry sets=\"true\" ny=\"2\"/>", 3,
             "700 1 0 0\n"
             "3 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0 0\n"
             "1 0 0\n"
             "0 0.5 0\n"
             "1 0.5 0\n"
             "0 0 1\n"
             "1 0 1\n"
             "0 0.5 1\n"
             "1 0.5 1\n"
             "700 1 0 0\n"
             "3 0\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "2 2\n"
             "0 0 1 1\n"
             "0 0.5 0\n"
             "1 0.5 0\n"
             "0 1 0\n"
             "1 1 0\n"
             "0 0.5 1\n"
             "1 0.5 1\n"
             "0 1 1\n"
             "1 1 1\n",
             "Boundary: 1 1 2 1 2 2 1 3 2 1 5 2 1 6 2 2 1 2 2 2 2 2 4 2 2 5 2 2 6 2 \n"
             "BoundaryX: 1 1 2 1 2 2 2 1 2 2 2 2 \n"
             "BoundaryY: 1 3 2 2 4 2 \n"
             "BoundaryZ: 1 5 2 1 6 2 2 5 2 2 6 2 \n"
             "Corners: 1 1 0 1 2 0 1 5 0 1 6 0 2 3 0 2 4 0 2 7 0 2 8 0 \n"
             "Edge1: 1 1 1 2 1 1 \n"
             "Edge10: 2 10 1 \n"
             "Edge11: 1 11 1 \n"
             "Edge12: 2 12 1 \n"
             "Edge2: 1 2 1 2 2 1 \n"
             "Edge3: 1 3 1 2 3 1 \n"
             "Edge4: 1 4 1 2 4 1 \n"
             "Edge5: 1 5 1 \n"
             "Edge6: 1 6 1 \n"
             "Edge7: 2 7 1 \n"
             "Edge8: 2 8 1 \n"
             "Edge9: 1 9 1 \n"
             "Face1: 1 1 2 2 1 2 \n"
             "Face1Patches: 1 0 3 2 0 3 \n"
             "Face2: 1 2 2 2 2 2 \n"
             "Face2Patches: 1 0 3 2 0 3 \n"
             "Face3: 1 3 2 \n"
             "Face3Patches: 1 0 3 \n"
             "Face4: 2 4 2 \n"
             "Face4Patches: 2 0 3 \n"
             "Face5: 1 5 2 2 5 2 \n"
             "Face5Patches: 1 0 3 2 0 3 \n"
             "Face6: 1 6 2 2 6 2 \n"
             "Face6Patches: 1 0 3 2 0 3 \n"
             "Frame: 1 1 1 1 2 1 1 3 1 1 4 1 1 5 1 1 6 1 1 9 1 1 11 1 2 1 1 2 2 1 2 3 1 2 4 1 2 7 1 2 8 1 2 10 1 2 12 1 \n"
             "InnerPatches: \n"
             "Vertex1: 1 1 0 \n"
             "Vertex2: 1 2 0 \n"
             "Vertex3: 2 3 0 \n"
             "Vertex4: 2 4 0 \n"
             "Vertex5: 1 5 0 \n"
             "Vertex6: 1 6 0 \n"
             "Vertex7: 2 7 0 \n"
             "Vertex8: 2 8 0 \n"},

    GeomTest{"Multiple Z",
             "<geometry sets=\"true\" nz=\"2\"/>", 3,
              "700 1 0 0\n"
              "3 0\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "0 0 0\n"
              "1 0 0\n"
              "0 1 0\n"
              "1 1 0\n"
              "0 0 0.5\n"
              "1 0 0.5\n"
              "0 1 0.5\n"
              "1 1 0.5\n"
              "700 1 0 0\n"
              "3 0\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "0 0 0.5\n"
              "1 0 0.5\n"
              "0 1 0.5\n"
              "1 1 0.5\n"
              "0 0 1\n"
              "1 0 1\n"
              "0 1 1\n"
              "1 1 1\n",
              "Boundary: 1 1 2 1 2 2 1 3 2 1 4 2 1 5 2 2 1 2 2 2 2 2 3 2 2 4 2 2 6 2 \n"
              "BoundaryX: 1 1 2 1 2 2 2 1 2 2 2 2 \n"
              "BoundaryY: 1 3 2 1 4 2 2 3 2 2 4 2 \n"
              "BoundaryZ: 1 5 2 2 6 2 \n"
              "Corners: 1 1 0 1 2 0 1 3 0 1 4 0 2 5 0 2 6 0 2 7 0 2 8 0 \n"
              "Edge1: 1 1 1 \n"
              "Edge10: 1 10 1 \n"
              "Edge11: 2 11 1 \n"
              "Edge12: 2 12 1 \n"
              "Edge2: 1 2 1 \n"
              "Edge3: 2 3 1 \n"
              "Edge4: 2 4 1 \n"
              "Edge5: 1 5 1 2 5 1 \n"
              "Edge6: 1 6 1 2 6 1 \n"
              "Edge7: 1 7 1 2 7 1 \n"
              "Edge8: 1 8 1 2 8 1 \n"
              "Edge9: 1 9 1 \n"
              "Face1: 1 1 2 2 1 2 \n"
              "Face1Patches: 1 0 3 2 0 3 \n"
              "Face2: 1 2 2 2 2 2 \n"
              "Face2Patches: 1 0 3 2 0 3 \n"
              "Face3: 1 3 2 2 3 2 \n"
              "Face3Patches: 1 0 3 2 0 3 \n"
              "Face4: 1 4 2 2 4 2 \n"
              "Face4Patches: 1 0 3 2 0 3 \n"
              "Face5: 1 5 2 \n"
              "Face5Patches: 1 0 3 \n"
              "Face6: 2 6 2 \n"
              "Face6Patches: 2 0 3 \n"
              "Frame: 1 1 1 1 2 1 1 5 1 1 6 1 1 7 1 1 8 1 1 9 1 1 10 1 2 3 1 2 4 1 2 5 1 2 6 1 2 7 1 2 8 1 2 11 1 2 12 1 \n"
              "InnerPatches: \n"
              "Vertex1: 1 1 0 \n"
              "Vertex2: 1 2 0 \n"
              "Vertex3: 1 3 0 \n"
              "Vertex4: 1 4 0 \n"
              "Vertex5: 2 5 0 \n"
              "Vertex6: 2 6 0 \n"
              "Vertex7: 2 7 0 \n"
              "Vertex8: 2 8 0 \n"},

    GeomTest{"Multiple",
             "<geometry sets=\"true\" nx=\"2\" ny=\"2\" nz=\"2\"/>", 3,
              "700 1 0 0\n"
              "3 0\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "0 0 0\n"
              "0.5 0 0\n"
              "0 0.5 0\n"
              "0.5 0.5 0\n"
              "0 0 0.5\n"
              "0.5 0 0.5\n"
              "0 0.5 0.5\n"
              "0.5 0.5 0.5\n"
              "700 1 0 0\n"
              "3 0\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "0.5 0 0\n"
              "1 0 0\n"
              "0.5 0.5 0\n"
              "1 0.5 0\n"
              "0.5 0 0.5\n"
              "1 0 0.5\n"
              "0.5 0.5 0.5\n"
              "1 0.5 0.5\n"
              "700 1 0 0\n"
              "3 0\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "0 0.5 0\n"
              "0.5 0.5 0\n"
              "0 1 0\n"
              "0.5 1 0\n"
              "0 0.5 0.5\n"
              "0.5 0.5 0.5\n"
              "0 1 0.5\n"
              "0.5 1 0.5\n"
              "700 1 0 0\n"
              "3 0\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "0.5 0.5 0\n"
              "1 0.5 0\n"
              "0.5 1 0\n"
              "1 1 0\n"
              "0.5 0.5 0.5\n"
              "1 0.5 0.5\n"
              "0.5 1 0.5\n"
              "1 1 0.5\n"
              "700 1 0 0\n"
              "3 0\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "0 0 0.5\n"
              "0.5 0 0.5\n"
              "0 0.5 0.5\n"
              "0.5 0.5 0.5\n"
              "0 0 1\n"
              "0.5 0 1\n"
              "0 0.5 1\n"
              "0.5 0.5 1\n"
              "700 1 0 0\n"
              "3 0\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "0.5 0 0.5\n"
              "1 0 0.5\n"
              "0.5 0.5 0.5\n"
              "1 0.5 0.5\n"
              "0.5 0 1\n"
              "1 0 1\n"
              "0.5 0.5 1\n"
              "1 0.5 1\n"
              "700 1 0 0\n"
              "3 0\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "0 0.5 0.5\n"
              "0.5 0.5 0.5\n"
              "0 1 0.5\n"
              "0.5 1 0.5\n"
              "0 0.5 1\n"
              "0.5 0.5 1\n"
              "0 1 1\n"
              "0.5 1 1\n"
              "700 1 0 0\n"
              "3 0\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "2 2\n"
              "0 0 1 1\n"
              "0.5 0.5 0.5\n"
              "1 0.5 0.5\n"
              "0.5 1 0.5\n"
              "1 1 0.5\n"
              "0.5 0.5 1\n"
              "1 0.5 1\n"
              "0.5 1 1\n"
              "1 1 1\n",
              "Boundary: 1 1 2 1 3 2 1 5 2 "
              "2 2 2 2 3 2 2 5 2 "
              "3 1 2 3 4 2 3 5 2 "
              "4 2 2 4 4 2 4 5 2 "
              "5 1 2 5 3 2 5 6 2 "
              "6 2 2 6 3 2 6 6 2 "
              "7 1 2 7 4 2 7 6 2 "
              "8 2 2 8 4 2 8 6 2 \n"
              "BoundaryX: 1 1 2 2 2 2 3 1 2 4 2 2 5 1 2 6 2 2 7 1 2 8 2 2 \n"
              "BoundaryY: 1 3 2 2 3 2 3 4 2 4 4 2 5 3 2 6 3 2 7 4 2 8 4 2 \n"
              "BoundaryZ: 1 5 2 2 5 2 3 5 2 4 5 2 5 6 2 6 6 2 7 6 2 8 6 2 \n"
              "Corners: 1 1 0 2 2 0 3 3 0 4 4 0 5 5 0 6 6 0 7 7 0 8 8 0 \n"
              "Edge1: 1 1 1 3 1 1 \n"
              "Edge10: 3 10 1 4 10 1 \n"
              "Edge11: 5 11 1 6 11 1 \n"
              "Edge12: 7 12 1 8 12 1 \n"
              "Edge2: 2 2 1 4 2 1 \n"
              "Edge3: 5 3 1 7 3 1 \n"
              "Edge4: 6 4 1 8 4 1 \n"
              "Edge5: 1 5 1 5 5 1 \n"
              "Edge6: 2 6 1 6 6 1 \n"
              "Edge7: 3 7 1 7 7 1 \n"
              "Edge8: 4 8 1 8 8 1 \n"
              "Edge9: 1 9 1 2 9 1 \n"
              "Face1: 1 1 2 3 1 2 5 1 2 7 1 2 \n"
              "Face1Patches: 1 0 3 3 0 3 5 0 3 7 0 3 \n"
              "Face2: 2 2 2 4 2 2 6 2 2 8 2 2 \n"
              "Face2Patches: 2 0 3 4 0 3 6 0 3 8 0 3 \n"
              "Face3: 1 3 2 2 3 2 5 3 2 6 3 2 \n"
              "Face3Patches: 1 0 3 2 0 3 5 0 3 6 0 3 \n"
              "Face4: 3 4 2 4 4 2 7 4 2 8 4 2 \n"
              "Face4Patches: 3 0 3 4 0 3 7 0 3 8 0 3 \n"
              "Face5: 1 5 2 2 5 2 3 5 2 4 5 2 \n"
              "Face5Patches: 1 0 3 2 0 3 3 0 3 4 0 3 \n"
              "Face6: 5 6 2 6 6 2 7 6 2 8 6 2 \n"
              "Face6Patches: 5 0 3 6 0 3 7 0 3 8 0 3 \n"
              "Frame: 1 1 1 1 5 1 1 9 1 "
              "2 2 1 2 6 1 2 9 1 "
              "3 1 1 3 7 1 3 10 1 "
              "4 2 1 4 8 1 4 10 1 "
              "5 3 1 5 5 1 5 11 1 "
              "6 4 1 6 6 1 6 11 1 "
              "7 3 1 7 7 1 7 12 1 "
              "8 4 1 8 8 1 8 12 1 \n"
              "InnerPatches: \n"
              "Vertex1: 1 1 0 \n"
              "Vertex2: 2 2 0 \n"
              "Vertex3: 3 3 0 \n"
              "Vertex4: 4 4 0 \n"
              "Vertex5: 5 5 0 \n"
              "Vertex6: 6 6 0 \n"
              "Vertex7: 7 7 0 \n"
              "Vertex8: 8 8 0 \n"}
  );

  SECTION(param.name) {
    tinyxml2::XMLDocument doc;
    doc.Parse(param.xml);
    TestModelGeneratorWrapper<MultiPatchModelGenerator3D> gen(doc.RootElement());
    SIM3D sim;
    gen.createTopologySets(sim);
    DoTest(param, gen.createG2(param.dim), sim.getTopology());
  }
}


TEST_CASE("TestMultiPatchModelGenerator.Generate3DLR")
{
  SIMMultiPatch3D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  REQUIRE(sim.read("refdata/modelgen3d_lr.xinp"));
  sim.preprocess();
  REQUIRE(sim.getNoNodes() == 112);
}


TEST_CASE("TestMultiPatchModelGenerator.InnerPatches3D")
{
  tinyxml2::XMLDocument doc;
  doc.Parse("<geometry nx=\"3\" ny=\"3\" nz=\"3\"  sets=\"true\"/>");
  TestModelGeneratorWrapper<MultiPatchModelGenerator3D> gen(doc.RootElement());
  SIM3D sim;
  gen.createTopologySets(sim);
  REQUIRE(sim.topology("InnerPatches").size() == 1);
  REQUIRE(sim.topology("InnerPatches").begin()->patch == 14);
}
