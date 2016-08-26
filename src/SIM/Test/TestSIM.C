//==============================================================================
//!
//! \file TestSIM.C
//!
//! \date Apr 27 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for various SIM functionality.
//!
//==============================================================================

#include "ASMs3D.h"
#include "IntegrandBase.h"
#include "IFEM.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"

#include "gtest/gtest.h"
#include "tinyxml.h"


template<class Dim>
class CreateDefaultGeometrySIM : public Dim
{
public:
  CreateDefaultGeometrySIM(unsigned char n1) : Dim(n1)
  {}
  CreateDefaultGeometrySIM(const std::vector<unsigned char>& nf) : Dim(nf)
  {}

  virtual ~CreateDefaultGeometrySIM()
  {
    Dim::myProblem = nullptr;
  }

  std::string createDefaultGeom(const TiXmlElement* geo)
  {
    return Dim::createDefaultG2(geo);
  }

  bool createDefaultTop(const TiXmlElement* geo)
  {
    Dim::myModel = Dim::createDefaultGeometry(geo);
    Dim::preprocess();
    return Dim::createDefaultTopology(geo);
  }

  TopologySet createDefaultTopSets(const TiXmlElement* geo)
  {
    return Dim::createDefaultTopologySets(geo);
  }

  void setIntegrand(IntegrandBase* itg)
  {
    Dim::myProblem = itg;
  }
};


class TestProjectIntegrand : public IntegrandBase {
public:
  TestProjectIntegrand(int dim) : IntegrandBase(dim) {}

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  bool evalSol(Vector& s, const FiniteElement& fe, const Vec3& X,
                       const std::vector<int>& MNPC) const override
  {
    s.resize(1);
    s(1) = X[0] + X[1] + X[2];

    return true;
  }

  size_t getNoFields(int = 2) const override { return 1; }
};


class TestSIM1D : public testing::Test,
                  public testing::WithParamInterface<std::vector<std::string>>
{
};


class TestSIM2D : public testing::Test,
                  public testing::WithParamInterface<std::vector<std::string>>
{
};


class TestSIM3D : public testing::Test,
                  public testing::WithParamInterface<std::vector<std::string>>
{
};



TEST(TestSIM, UniqueBoundaryNodes)
{
  SIM2D sim(1);
  ASSERT_TRUE(sim.read("src/SIM/Test/refdata/boundary_nodes.xinp"));
  sim.preprocess();

  int bcode = sim.getUniquePropertyCode("dir",0);
  std::vector<int> vec;
  sim.getBoundaryNodes(bcode, vec);

  std::sort(vec.begin(), vec.end());
  ASSERT_TRUE(std::unique(vec.begin(), vec.end()) == vec.end());
}

auto&& DoTest = [](const std::vector<std::string>& param, const std::string gen, const TopologySet& sets)
{
  ASSERT_STREQ(gen.c_str(), param[1].c_str());
  if (param.size() > 2) {
    std::string gen;
    for (auto& it : sets) {
      gen += it.first + ": ";
      for (auto& it2 : it.second) {
        std::stringstream str;
        str << it2.patch << " " << it2.item << " " << it2.idim << " ";
        gen += str.str();
      }
      gen += "\n";
    }
    ASSERT_STREQ(gen.c_str(), param[2].c_str());
  }
};


TEST_P(TestSIM1D, CreateDefaultGeometry)
{
  CreateDefaultGeometrySIM<SIM1D> sim(1);
  TiXmlDocument doc;
  doc.Parse(GetParam()[0].c_str());
  std::string gen = sim.createDefaultGeom(doc.RootElement());
  TopologySet sets = sim.createDefaultTopSets(doc.RootElement());
  DoTest(GetParam(), gen, sets);
}


const std::vector<std::vector<std::string>> geometry1D =
  {{"<geometry/>", "100 1 0 0\n"
                   "1 0\n"
                   "2 2\n"
                   "0 0 1 1\n"
                   "0.0\n"
                   "1.0\n",
                   "Boundary: 1 1 0 1 2 0 \n"
                   "Corners: 1 1 0 1 2 0 \n"
                   "Vertex1: 1 1 0 \n"
                   "Vertex2: 1 2 0 \n"}};


INSTANTIATE_TEST_CASE_P(TestSIM, TestSIM1D, testing::ValuesIn(geometry1D));


TEST_P(TestSIM2D, CreateDefaultGeometry)
{
  CreateDefaultGeometrySIM<SIM2D> sim(1);
  TiXmlDocument doc;
  doc.Parse(GetParam()[0].c_str());
  std::string gen = sim.createDefaultGeom(doc.RootElement());
  TopologySet sets = sim.createDefaultTopSets(doc.RootElement());
  DoTest(GetParam(), gen, sets);
}


const std::vector<std::vector<std::string>> geometry2D =
  {{"<geometry/>", "200 1 0 0\n"
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

    {"<geometry rational=\"1\"/>",
                    "200 1 0 0\n"
                    "2 1\n"
                    "2 2\n"
                    "0 0 1 1\n"
                    "2 2\n"
                    "0 0 1 1\n"
                    "0 0 1.0\n"
                    "1 0 1.0\n"
                    "0 1 1.0\n"
                    "1 1 1.0\n"},

    {"<geometry scale=\"2\"/>",
                    "200 1 0 0\n"
                    "2 0\n"
                    "2 2\n"
                    "0 0 1 1\n"
                    "2 2\n"
                    "0 0 1 1\n"
                    "0 0\n"
                    "2 0\n"
                    "0 2\n"
                    "2 2\n"},

    {"<geometry X0=\"2 0\"/>",
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

    {"<geometry X0=\"0 2\"/>",
                    "200 1 0 0\n"
                    "2 0\n"
                    "2 2\n"
                    "0 0 1 1\n"
                    "2 2\n"
                    "0 0 1 1\n"
                    "0 2\n"
                    "1 2\n"
                    "0 3\n"
                    "1 3\n"},

    {"<geometry Lx=\"2\"/>",
                    "200 1 0 0\n"
                    "2 0\n"
                    "2 2\n"
                    "0 0 1 1\n"
                    "2 2\n"
                    "0 0 1 1\n"
                    "0 0\n"
                    "2 0\n"
                    "0 1\n"
                    "2 1\n"},

    {"<geometry Ly=\"2\"/>",
                    "200 1 0 0\n"
                    "2 0\n"
                    "2 2\n"
                    "0 0 1 1\n"
                    "2 2\n"
                    "0 0 1 1\n"
                    "0 0\n"
                    "1 0\n"
                    "0 2\n"
                    "1 2\n"},

    {"<geometry nx=\"2\"/>",
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
                    "Corners: 1 1 0 1 3 0 2 2 0 2 4 0 \n"
                    "Edge1: 1 1 1 \n"
                    "Edge2: 2 2 1 \n"
                    "Edge3: 1 3 1 2 3 1 \n"
                    "Edge4: 1 4 1 2 4 1 \n"
                    "Vertex1: 1 1 0 \n"
                    "Vertex2: 2 2 0 \n"
                    "Vertex3: 1 3 0 \n"
                    "Vertex4: 2 4 0 \n"},

    {"<geometry ny=\"2\"/>",
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
                    "Corners: 1 1 0 1 2 0 2 3 0 2 4 0 \n"
                    "Edge1: 1 1 1 2 1 1 \n"
                    "Edge2: 1 2 1 2 2 1 \n"
                    "Edge3: 1 3 1 \n"
                    "Edge4: 2 4 1 \n"
                    "Vertex1: 1 1 0 \n"
                    "Vertex2: 1 2 0 \n"
                    "Vertex3: 2 3 0 \n"
                    "Vertex4: 2 4 0 \n"},

    {"<geometry nx=\"2\" ny=\"2\"/>",
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
                    "Corners: 1 1 0 2 2 0 3 3 0 4 4 0 \n"
                    "Edge1: 1 1 1 3 1 1 \n"
                    "Edge2: 2 2 1 4 2 1 \n"
                    "Edge3: 1 3 1 2 3 1 \n"
                    "Edge4: 3 4 1 4 4 1 \n"
                    "Vertex1: 1 1 0 \n"
                    "Vertex2: 2 2 0 \n"
                    "Vertex3: 3 3 0 \n"
                    "Vertex4: 4 4 0 \n"}};

INSTANTIATE_TEST_CASE_P(TestSIM, TestSIM2D, testing::ValuesIn(geometry2D));


typedef std::vector<std::pair<std::array<int,2>,std::array<int,2>>> ConnectionPairVec;

const std::pair<std::string, ConnectionPairVec> periodic2D =
  {"<geometry nx=\"3\" ny=\"3\" periodic_x=\"1\" periodic_y=\"1\"/>",
                                       {{{{1,2},{2,1}},
                                         {{2,2},{3,1}},
                                         {{4,2},{5,1}},
                                         {{5,2},{6,1}},
                                         {{7,2},{8,1}},
                                         {{8,2},{9,1}},
                                         {{1,4},{4,3}},
                                         {{2,4},{5,3}},
                                         {{3,4},{6,3}},
                                         {{4,4},{7,3}},
                                         {{5,4},{8,3}},
                                         {{6,4},{9,3}},
                                         {{1,1},{3,2}},
                                         {{4,1},{6,2}},
                                         {{7,1},{9,2}},
                                         {{1,3},{7,4}},
                                         {{2,3},{8,4}},
                                         {{3,3},{9,4}}}}};


TEST(TestSIM2D, PeriodicConnections)
{
  CreateDefaultGeometrySIM<SIM2D> sim(1);
  TiXmlDocument doc;
  doc.Parse(periodic2D.first.c_str());
  std::stringstream str;
  IFEM::cout.setStream(str);
  sim.createDefaultTop(doc.RootElement());
  std::string text = str.str();
  for (auto& it : periodic2D.second) {
    std::stringstream str;
    str << "Connecting P" << it.second[0] << " E" << it.second[1] << " to P" << it.first[0] << " E" << it.first[1];
    auto pos = text.find(str.str());
    ASSERT_TRUE(pos != std::string::npos);
  }

  for (int p = 1; p <= 3; ++p)
    for (int n = 1; n <= 2; ++n)
      ASSERT_EQ(sim.getPatch(3*(p-1)+1)->getNodeID(2*(n-1)+1),
                sim.getPatch(3*p)->getNodeID(2*n));

  for (int p = 1; p <= 3; ++p)
    for (int n = 1; n <= 2; ++n)
      ASSERT_EQ(sim.getPatch(p)->getNodeID(n),
                sim.getPatch(6+p)->getNodeID(2+n));
}


TEST(TestSIM2D, PeriodicConnectionsMixed)
{
  CreateDefaultGeometrySIM<SIM2D> sim({1,1});
  TiXmlDocument doc;
  doc.Parse(periodic2D.first.c_str());
  IFEM::cout.setStream(std::cout);
  sim.createDefaultTop(doc.RootElement());

  for (int p = 1; p <= 3; ++p)
    for (int n = 1; n <= 2; ++n)
      ASSERT_NE(sim.getPatch(3*(p-1)+1)->getNodeID(9+2*(n-1)+1),
                sim.getPatch(3*p)->getNodeID(9+2*n));

  for (int p = 1; p <= 3; ++p)
    for (int n = 1; n <= 2; ++n)
      ASSERT_NE(sim.getPatch(p)->getNodeID(9+n),
                sim.getPatch(6+p)->getNodeID(9+2+n));
}


TEST_P(TestSIM3D, CreateDefaultGeometry)
{
  CreateDefaultGeometrySIM<SIM3D> sim(1);
  TiXmlDocument doc;
  doc.Parse(GetParam()[0].c_str());
  std::string gen = sim.createDefaultGeom(doc.RootElement());
  TopologySet sets = sim.createDefaultTopSets(doc.RootElement());
  DoTest(GetParam(), gen, sets);
}


const std::vector<std::vector<std::string>> geometry3D =
  {{"<geometry/>", "700 1 0 0\n"
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

  {"<geometry rational=\"1\"/>",
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
                   "1 1 1 1.0\n"},

  {"<geometry scale=\"2\"/>",
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
                   "2 2 2\n"},

  {"<geometry X0=\"2 0 0\"/>",
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
                   "3 1 1\n"},

  {"<geometry X0=\"0 2 0\"/>",
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
                   "1 3 1\n"},

  {"<geometry X0=\"0 0 2\"/>",
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
                   "1 1 3\n"},

  {"<geometry Lx=\"2\"/>",
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
                   "2 1 1\n"},

  {"<geometry Ly=\"2\"/>",
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
                   "1 2 1\n"},

  {"<geometry Lz=\"2\"/>",
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
                   "1 1 2\n"},

  {"<geometry nx=\"2\"/>",
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
                   "Face2: 2 2 2 \n"
                   "Face3: 1 3 2 2 3 2 \n"
                   "Face4: 1 4 2 2 4 2 \n"
                   "Face5: 1 5 2 2 5 2 \n"
                   "Face6: 1 6 2 2 6 2 \n"
                   "Frame: 1 1 1 1 3 1 1 5 1 1 7 1 1 9 1 1 10 1 1 11 1 1 12 1 2 2 1 2 4 1 2 6 1 2 8 1 2 9 1 2 10 1 2 11 1 2 12 1 \n"
                   "Vertex1: 1 1 0 \n"
                   "Vertex2: 2 2 0 \n"
                   "Vertex3: 1 3 0 \n"
                   "Vertex4: 2 4 0 \n"
                   "Vertex5: 1 5 0 \n"
                   "Vertex6: 2 6 0 \n"
                   "Vertex7: 1 7 0 \n"
                   "Vertex8: 2 8 0 \n"},

  {"<geometry ny=\"2\"/>",
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
                   "Face2: 1 2 2 2 2 2 \n"
                   "Face3: 1 3 2 \n"
                   "Face4: 2 4 2 \n"
                   "Face5: 1 5 2 2 5 2 \n"
                   "Face6: 1 6 2 2 6 2 \n"
                   "Frame: 1 1 1 1 2 1 1 3 1 1 4 1 1 5 1 1 6 1 1 9 1 1 11 1 2 1 1 2 2 1 2 3 1 2 4 1 2 7 1 2 8 1 2 10 1 2 12 1 \n"
                   "Vertex1: 1 1 0 \n"
                   "Vertex2: 1 2 0 \n"
                   "Vertex3: 2 3 0 \n"
                   "Vertex4: 2 4 0 \n"
                   "Vertex5: 1 5 0 \n"
                   "Vertex6: 1 6 0 \n"
                   "Vertex7: 2 7 0 \n"
                   "Vertex8: 2 8 0 \n"},

   {"<geometry nz=\"2\"/>",
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
                   "Face2: 1 2 2 2 2 2 \n"
                   "Face3: 1 3 2 2 3 2 \n"
                   "Face4: 1 4 2 2 4 2 \n"
                   "Face5: 1 5 2 \n"
                   "Face6: 2 6 2 \n"
                   "Frame: 1 1 1 1 2 1 1 5 1 1 6 1 1 7 1 1 8 1 1 9 1 1 10 1 2 3 1 2 4 1 2 5 1 2 6 1 2 7 1 2 8 1 2 11 1 2 12 1 \n"
                   "Vertex1: 1 1 0 \n"
                   "Vertex2: 1 2 0 \n"
                   "Vertex3: 1 3 0 \n"
                   "Vertex4: 1 4 0 \n"
                   "Vertex5: 2 5 0 \n"
                   "Vertex6: 2 6 0 \n"
                   "Vertex7: 2 7 0 \n"
                   "Vertex8: 2 8 0 \n"},

  {"<geometry nx=\"2\" ny=\"2\" nz=\"2\"/>",
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
                   "Face2: 2 2 2 4 2 2 6 2 2 8 2 2 \n"
                   "Face3: 1 3 2 2 3 2 5 3 2 6 3 2 \n"
                   "Face4: 3 4 2 4 4 2 7 4 2 8 4 2 \n"
                   "Face5: 1 5 2 2 5 2 3 5 2 4 5 2 \n"
                   "Face6: 5 6 2 6 6 2 7 6 2 8 6 2 \n"
                   "Frame: 1 1 1 1 5 1 1 9 1 "
                          "2 2 1 2 6 1 2 9 1 "
                          "3 1 1 3 7 1 3 10 1 "
                          "4 2 1 4 8 1 4 10 1 "
                          "5 3 1 5 5 1 5 11 1 "
                          "6 4 1 6 6 1 6 11 1 "
                          "7 3 1 7 7 1 7 12 1 "
                          "8 4 1 8 8 1 8 12 1 \n"
                   "Vertex1: 1 1 0 \n"
                   "Vertex2: 2 2 0 \n"
                   "Vertex3: 3 3 0 \n"
                   "Vertex4: 4 4 0 \n"
                   "Vertex5: 5 5 0 \n"
                   "Vertex6: 6 6 0 \n"
                   "Vertex7: 7 7 0 \n"
                   "Vertex8: 8 8 0 \n"}};


INSTANTIATE_TEST_CASE_P(TestSIM, TestSIM3D, testing::ValuesIn(geometry3D));


const std::pair<std::string, ConnectionPairVec> periodic3D =
  {"<geometry nx=\"3\" ny=\"3\" nz=\"3\" periodic_x=\"1\" periodic_y=\"1\" periodic_z=\"1\"/>",
                                        {{{{1,2}, {2,1}},
                                          {{2,2}, {3,1}},
                                          {{4,2}, {5,1}},
                                          {{5,2}, {6,1}},
                                          {{7,2}, {8,1}},
                                          {{8,2}, {9,1}},
                                         {{10,2},{11,1}},
                                         {{11,2},{12,1}},
                                         {{13,2},{14,1}},
                                         {{14,2},{15,1}},
                                         {{16,2},{17,1}},
                                         {{17,2},{18,1}},
                                         {{19,2},{20,1}},
                                         {{20,2},{21,1}},
                                         {{22,2},{23,1}},
                                         {{23,2},{24,1}},
                                         {{25,2},{26,1}},
                                         {{26,2},{27,1}},
                                         {{ 1,4}, { 4,3}},
                                         {{ 4,4}, { 7,3}},
                                         {{ 2,4}, { 5,3}},
                                         {{ 5,4}, { 8,3}},
                                         {{ 3,4}, { 6,3}},
                                         {{ 6,4}, { 9,3}},
                                         {{10,4}, {13,3}},
                                         {{13,4}, {16,3}},
                                         {{11,4}, {14,3}},
                                         {{14,4}, {17,3}},
                                         {{12,4}, {15,3}},
                                         {{15,4}, {18,3}},
                                         {{19,4}, {22,3}},
                                         {{22,4}, {25,3}},
                                         {{20,4}, {23,3}},
                                         {{23,4}, {26,3}},
                                         {{21,4}, {24,3}},
                                         {{24,4}, {27,3}},
                                         {{ 1,6}, {10,5}},
                                         {{ 2,6}, {11,5}},
                                         {{ 3,6}, {12,5}},
                                         {{ 4,6}, {13,5}},
                                         {{ 5,6}, {14,5}},
                                         {{ 6,6}, {15,5}},
                                         {{ 7,6}, {16,5}},
                                         {{ 8,6}, {17,5}},
                                         {{ 9,6}, {18,5}},
                                         {{10,6}, {19,5}},
                                         {{11,6}, {20,5}},
                                         {{12,6}, {21,5}},
                                         {{13,6}, {22,5}},
                                         {{14,6}, {23,5}},
                                         {{15,6}, {24,5}},
                                         {{16,6}, {25,5}},
                                         {{17,6}, {26,5}},
                                         {{18,6}, {27,5}}}}};


TEST(TestSIM3D, PeriodicConnections)
{
  CreateDefaultGeometrySIM<SIM3D> sim(1);
  TiXmlDocument doc;
  doc.Parse(periodic3D.first.c_str());
  std::stringstream str;
  IFEM::cout.setStream(str);
  sim.createDefaultTop(doc.RootElement());
  std::string text = str.str();
  for (auto& it : periodic3D.second) {
    std::stringstream str;
    str << "Connecting P" << it.second[0] << " F" << it.second[1] << " to P" << it.first[0] << " F" << it.first[1];
    auto pos = text.find(str.str());
    ASSERT_TRUE(pos != std::string::npos);
  }

  auto&& IJK = [](int i, int j, int k) { return 1 + (2*k + j)*2 + i; };

  for (int p = 1; p <= 9; ++p)
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        ASSERT_EQ(sim.getPatch(3*(p-1)+1)->getNodeID(IJK(0, j, k)),
                  sim.getPatch(3*p)->getNodeID(IJK(1, j, k)));

  int r = 1;
  for (int q = 1; q <= 3; ++q) {
    for (int p = 1; p <= 3; ++p, ++r)
      for (int i = 0; i < 2; ++i)
        for (int k = 0; k < 2; ++k) {
          ASSERT_EQ(sim.getPatch(r)->getNodeID(IJK(i, 0, k)),
                    sim.getPatch(r+6)->getNodeID(IJK(i, 1, k)));
        }
    r += 6;
  }

  for (int p = 1; p <= 9; ++p)
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        ASSERT_EQ(sim.getPatch(p)->getNodeID(IJK(i, j, 0)),
                  sim.getPatch(p+18)->getNodeID(IJK(i, j, 1)));
}


TEST(TestSIM3D, PeriodicConnectionsMixed)
{
  CreateDefaultGeometrySIM<SIM3D> sim({1,1});
  TiXmlDocument doc;
  doc.Parse(periodic3D.first.c_str());
  IFEM::cout.setStream(std::cout);
  sim.createDefaultTop(doc.RootElement());

  auto&& IJK = [](int i, int j, int k) { return 1 + (2*k + j)*2 + i; };

  for (int p = 1; p <= 9; ++p)
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        ASSERT_NE(sim.getPatch(3*(p-1)+1)->getNodeID(27+IJK(0, j, k)),
                  sim.getPatch(3*p)->getNodeID(27+IJK(1, j, k)));

  int r = 1;
  for (int q = 1; q <= 3; ++q) {
    for (int p = 1; p <= 3; ++p, ++r)
      for (int i = 0; i < 2; ++i)
        for (int k = 0; k < 2; ++k) {
          ASSERT_NE(sim.getPatch(r)->getNodeID(27+IJK(i, 0, k)),
                    sim.getPatch(r+6)->getNodeID(27+IJK(i, 1, k)));
        }
    r += 6;
  }

  for (int p = 1; p <= 9; ++p)
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        ASSERT_NE(sim.getPatch(p)->getNodeID(27+IJK(i, j, 0)),
                  sim.getPatch(p+18)->getNodeID(27+IJK(i, j, 1)));
}


TEST(TestSIM3D, ProjectSolution)
{
  CreateDefaultGeometrySIM<SIM3D> sim({1});
  TestProjectIntegrand itg(3);
  sim.setIntegrand(&itg);
  TiXmlDocument doc;
  doc.Parse("<geometry/>");
  ASSERT_TRUE(sim.createDefaultTop(doc.RootElement()));

  Matrix ssol;
  ASSERT_TRUE(sim.project(ssol, Vector(sim.getNoDOFs())));

  size_t n = 1;
  for (size_t k = 0; k < 2; ++k)
    for (size_t j = 0; j < 2; ++j)
      for (size_t i = 0; i < 2; ++i)
      ASSERT_FLOAT_EQ(ssol(1, n++), i + j + k);
}


TEST(TestSIM3D, ProjectSolutionMixed)
{
  CreateDefaultGeometrySIM<SIM3D> sim({1,1});
  TestProjectIntegrand itg(3);
  sim.setIntegrand(&itg);
  TiXmlDocument doc;
  doc.Parse("<geometry/>");
  ASSERT_TRUE(sim.createDefaultTop(doc.RootElement()));

  Matrix ssol;
  ASSERT_TRUE(sim.project(ssol, Vector(sim.getNoDOFs())));

  size_t n = 1;
  for (size_t k = 0; k < 3; ++k)
    for (size_t j = 0; j < 3; ++j)
      for (size_t i = 0; i < 3; ++i)
      ASSERT_FLOAT_EQ(ssol(1, n++), i/2.0 + j/2.0 + k/2.0);
}
