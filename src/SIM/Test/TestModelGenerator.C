//==============================================================================
//!
//! \file TestModelGenerator.C
//!
//! \date Sep 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for default model generators.
//!
//==============================================================================

#include "IFEM.h"
#include "ModelGenerator.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "TopologySet.h"

#include "gtest/gtest.h"
#include "tinyxml.h"


template<class Generator>
class TestModelGeneratorWrapper : public Generator {
public:
  TestModelGeneratorWrapper(const TiXmlElement* geo) : Generator(geo) {}
  std::string createG2(int nsd)
  {
    return Generator::createG2(nsd);
  }
};

struct DefaultGeomTest {
  std::string xml;
  int dim;
  std::string g2;
  std::string sets;
};


class TestModelGenerator1D :
  public testing::Test,
  public testing::WithParamInterface<DefaultGeomTest>
{
};


class TestModelGenerator2D :
  public testing::Test,
  public testing::WithParamInterface<DefaultGeomTest>
{
};


class TestModelGenerator3D :
  public testing::Test,
  public testing::WithParamInterface<DefaultGeomTest>
{
};


auto&& DoTest = [](const DefaultGeomTest& ref, const std::string& gen,
                   const TopologySet& sets)
{
  ASSERT_STREQ(gen.c_str(), ref.g2.c_str());

  if (!ref.sets.empty()) {
    std::string gsets;
    for (auto& it : sets) {
      gsets += it.first + ": ";
      for (auto& it2 : it.second) {
        std::stringstream str;
        str << it2.patch << " " << it2.item << " " << it2.idim << " ";
        gsets += str.str();
      }
      gsets += "\n";
    }
    ASSERT_STREQ(gsets.c_str(), ref.sets.c_str());
  }
};


TEST_P(TestModelGenerator1D, Generate)
{
  TiXmlDocument doc;
  doc.Parse(GetParam().xml.c_str());
  TestModelGeneratorWrapper<DefaultGeometry1D> gen(doc.RootElement());
  std::string g2 = gen.createG2(GetParam().dim);
  SIM1D sim;
  TopologySet sets = gen.createTopologySets(sim);
  DoTest(GetParam(), g2, sets);
}


TEST_P(TestModelGenerator2D, Generate)
{
  TiXmlDocument doc;
  doc.Parse(GetParam().xml.c_str());
  TestModelGeneratorWrapper<DefaultGeometry2D> gen(doc.RootElement());
  std::string g2 = gen.createG2(GetParam().dim);
  SIM2D sim;
  TopologySet sets = gen.createTopologySets(sim);
  DoTest(GetParam(), g2, sets);
}


TEST_P(TestModelGenerator3D, Generate)
{
  TiXmlDocument doc;
  doc.Parse(GetParam().xml.c_str());
  TestModelGeneratorWrapper<DefaultGeometry3D> gen(doc.RootElement());
  std::string g2 = gen.createG2(GetParam().dim);
  SIM3D sim;
  TopologySet sets = gen.createTopologySets(sim);
  DoTest(GetParam(), g2, sets);
}


const std::vector<DefaultGeomTest> geometry1D =
  {{"<geometry sets=\"true\"/>", 1,
    "100 1 0 0\n"
    "1 0\n"
    "2 2\n"
    "0 0 1 1\n"
    "0.0\n"
    "1.0\n",
    "Boundary: 1 1 0 1 2 0 \n"
    "Corners: 1 1 0 1 2 0 \n"
    "Vertex1: 1 1 0 \n"
    "Vertex2: 1 2 0 \n"},

   {"<geometry/>", 3,
    "100 1 0 0\n"
    "3 0\n"
    "2 2\n"
    "0 0 1 1\n"
    "0.0 0.0 0.0\n"
    "1.0 0.0 0.0\n",
    ""},

   {"<geometry X0=\"1.0 1.0 0.0\" X1=\"1.0 2.0 0.0\"/>", 3,
    "100 1 0 0\n"
    "3 0\n"
    "2 2\n"
    "0 0 1 1\n"
    "1.0 1.0 0.0\n"
    "1.0 2.0 0.0\n",
    ""},

   {"<geometry L=\"2.0\"/>", 1,
    "100 1 0 0\n"
    "1 0\n"
    "2 2\n"
    "0 0 1 1\n"
    "0.0\n"
    "2.0\n",
    ""},

   {"<geometry rational=\"true\" sets=\"true\"/>", 1,
    "100 1 0 0\n"
    "1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "0.0 1.0\n"
    "1.0 1.0\n",
    "Boundary: 1 1 0 1 2 0 \n"
    "Corners: 1 1 0 1 2 0 \n"
    "Vertex1: 1 1 0 \n"
    "Vertex2: 1 2 0 \n"}};


INSTANTIATE_TEST_CASE_P(TestModelGenerator1D, TestModelGenerator1D, testing::ValuesIn(geometry1D));


const std::vector<DefaultGeomTest> geometry2D =
  {{"<geometry sets=\"true\"/>", 2,
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
    "Corners: 1 1 0 1 2 0 1 3 0 1 4 0 \n"
    "Edge1: 1 1 1 \n"
    "Edge2: 1 2 1 \n"
    "Edge3: 1 3 1 \n"
    "Edge4: 1 4 1 \n"
    "Vertex1: 1 1 0 \n"
    "Vertex2: 1 2 0 \n"
    "Vertex3: 1 3 0 \n"
    "Vertex4: 1 4 0 \n"},

   {"<geometry rational=\"1\"/>", 2,
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

   {"<geometry scale=\"2\"/>", 2,
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

   {"<geometry X0=\"2 0\"/>", 2,
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

    {"<geometry X0=\"0 2\"/>", 2,
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

    {"<geometry Lx=\"2\"/>", 2,
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

    {"<geometry Ly=\"2\"/>", 2,
     "200 1 0 0\n"
     "2 0\n"
     "2 2\n"
     "0 0 1 1\n"
     "2 2\n"
     "0 0 1 1\n"
     "0 0\n"
     "1 0\n"
     "0 2\n"
     "1 2\n", ""}};


INSTANTIATE_TEST_CASE_P(TestModelGenerator2D, TestModelGenerator2D, testing::ValuesIn(geometry2D));


const std::vector<DefaultGeomTest> geometry3D =
  {{"<geometry sets=\"true\"/>", 3,
    "700 1 0 0\n"
    "3 0\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "0 0 0 \n"
    "1 0 0 \n"
    "0 1 0 \n"
    "1 1 0 \n"
    "0 0 1 \n"
    "1 0 1 \n"
    "0 1 1 \n"
    "1 1 1 \n",
    "Boundary: 1 1 2 1 2 2 1 3 2 1 4 2 1 5 2 1 6 2 \n"
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
    "Face2: 1 2 2 \n"
    "Face3: 1 3 2 \n"
    "Face4: 1 4 2 \n"
    "Face5: 1 5 2 \n"
    "Face6: 1 6 2 \n"
    "Frame: 1 1 1 1 2 1 1 3 1 1 4 1 1 5 1 1 6 1 1 7 1 1 8 1 1 9 1 1 10 1 1 11 1 1 12 1 \n"
    "Vertex1: 1 1 0 \n"
    "Vertex2: 1 2 0 \n"
    "Vertex3: 1 3 0 \n"
    "Vertex4: 1 4 0 \n"
    "Vertex5: 1 5 0 \n"
    "Vertex6: 1 6 0 \n"
    "Vertex7: 1 7 0 \n"
    "Vertex8: 1 8 0 \n"},

  {"<geometry rational=\"1\"/>", 3,
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

  {"<geometry scale=\"2\"/>", 3,
   "700 1 0 0\n"
   "3 0\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "0 0 0 \n"
   "2 0 0 \n"
   "0 2 0 \n"
   "2 2 0 \n"
   "0 0 2 \n"
   "2 0 2 \n"
   "0 2 2 \n"
   "2 2 2 \n", ""},

  {"<geometry X0=\"2 0 0\"/>", 3,
   "700 1 0 0\n"
   "3 0\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 0 0 \n"
   "3 0 0 \n"
   "2 1 0 \n"
   "3 1 0 \n"
   "2 0 1 \n"
   "3 0 1 \n"
   "2 1 1 \n"
   "3 1 1 \n", ""},

  {"<geometry X0=\"0 2 0\"/>", 3,
   "700 1 0 0\n"
   "3 0\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "0 2 0 \n"
   "1 2 0 \n"
   "0 3 0 \n"
   "1 3 0 \n"
   "0 2 1 \n"
   "1 2 1 \n"
   "0 3 1 \n"
   "1 3 1 \n", ""},

  {"<geometry X0=\"0 0 2\"/>", 3,
   "700 1 0 0\n"
   "3 0\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "0 0 2 \n"
   "1 0 2 \n"
   "0 1 2 \n"
   "1 1 2 \n"
   "0 0 3 \n"
   "1 0 3 \n"
   "0 1 3 \n"
   "1 1 3 \n", ""},

  {"<geometry Lx=\"2\"/>", 3,
   "700 1 0 0\n"
   "3 0\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "0 0 0 \n"
   "2 0 0 \n"
   "0 1 0 \n"
   "2 1 0 \n"
   "0 0 1 \n"
   "2 0 1 \n"
   "0 1 1 \n"
   "2 1 1 \n", ""},

  {"<geometry Ly=\"2\"/>", 3,
   "700 1 0 0\n"
   "3 0\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "0 0 0 \n"
   "1 0 0 \n"
   "0 2 0 \n"
   "1 2 0 \n"
   "0 0 1 \n"
   "1 0 1 \n"
   "0 2 1 \n"
   "1 2 1 \n", ""},

  {"<geometry Lz=\"2\"/>", 3,
   "700 1 0 0\n"
     "3 0\n"
     "2 2\n"
     "0 0 1 1\n"
     "2 2\n"
     "0 0 1 1\n"
     "2 2\n"
     "0 0 1 1\n"
     "0 0 0 \n"
     "1 0 0 \n"
     "0 1 0 \n"
     "1 1 0 \n"
     "0 0 2 \n"
     "1 0 2 \n"
     "0 1 2 \n"
     "1 1 2 \n", ""}};


INSTANTIATE_TEST_CASE_P(TestModelGenerator3D, TestModelGenerator3D, testing::ValuesIn(geometry3D));
