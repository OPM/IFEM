//==============================================================================
//!
//! \file TestSIM.C
//!
//! \date Apr 27 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for various expected SIM behavior.
//!
//==============================================================================

#include "ASMbase.h"
#include "SIM2D.h"

#include "gtest/gtest.h"
#include "tinyxml.h"


template<class Dim>
class CreateDefaultGeometrySIM : public Dim
{
public:
  CreateDefaultGeometrySIM(unsigned char n1) : Dim(n1)
  {}

  std::string createDefaultGeom(const TiXmlElement* geo) const
  {
    return Dim::createDefaultG2(geo);
  }
};


class TestSIM2D : public testing::Test,
                public testing::WithParamInterface<std::pair<std::string,
                                                             std::string>>
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


TEST_P(TestSIM2D, CreateDefaultGeometry)
{
  CreateDefaultGeometrySIM<SIM2D> sim(1);
  TiXmlDocument doc;
  doc.Parse(GetParam().first.c_str());
  std::string gen = sim.createDefaultGeom(doc.RootElement());
  ASSERT_STREQ(gen.c_str(), GetParam().second.c_str());
}

const std::vector<std::pair<std::string, std::string>> tests2D =
  {{"<geometry/>", "200 1 0 0\n"
                    "2 0\n"
                    "2 2\n"
                    "0 0 1 1\n"
                    "2 2\n"
                    "0 0 1 1\n"
                    "0 0\n"
                    "1 0\n"
                    "0 1\n"
                    "1 1\n"},
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
                    "1 2\n"}};

INSTANTIATE_TEST_CASE_P(TestSIM, TestSIM2D, testing::ValuesIn(tests2D));
