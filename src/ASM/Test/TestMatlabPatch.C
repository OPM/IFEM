//==============================================================================
//!
//! \file TestMatlabPatch.C
//!
//! \date Aug 12 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Unit tests for output & input of Matlab mesh files.
//!
//==============================================================================

#include "SIM2D.h"
#include "ASMs2D.h"
#include "IntegrandBase.h"
#include "ModelGenerator.h"
#include "tinyxml.h"
#include <fstream>

#include "gtest/gtest.h"


class Dummy : public IntegrandBase {};


class SIM2D_default : public SIM2D
{
public:
  SIM2D_default(int nu, int nv) : SIM2D(new Dummy())
  {
    // Create a bi-unit square grid with (nu+1)x(nv+1) linear elements
    this->createDefaultModel();
    ASMs2D* pch1 = static_cast<ASMs2D*>(myModel.front());
    pch1->uniformRefine(0,nu);
    pch1->uniformRefine(1,nv);
    EXPECT_TRUE(this->preprocess());

    // Create topological boundary entities (vertices and edges)
    TiXmlDocument doc;
    doc.Parse("<geometry sets=\"true\"/>",nullptr,TIXML_ENCODING_UTF8);
    myEntitys = DefaultGeometry2D(doc.RootElement()).createTopologySets(*this);
    EXPECT_EQ(myEntitys.size(),10U);
  }
  virtual ~SIM2D_default() {}
};


class SIM2D_matlab : public SIM2D
{
public:
  SIM2D_matlab(const char* inputFile) : SIM2D(new Dummy())
  {
    std::string xml("<geometry><patchfile type=\"matlab\">");
    xml.append(inputFile);
    xml.append("</geometry>");

    // Read the FE model from the provided Matlab file
    TiXmlDocument doc;
    doc.Parse(xml.c_str(),nullptr,TIXML_ENCODING_UTF8);
    const TiXmlElement* tag = doc.RootElement();
    EXPECT_TRUE(tag != nullptr);
    EXPECT_TRUE(this->parse(tag));
    EXPECT_TRUE(this->preprocess());
  }
  virtual ~SIM2D_matlab() {}
};


TEST(TestMatlabPatch, IO)
{
  // Create a 5x4 element mesh and write it to a matlab file
  SIM2D_default sim1(4,3);
  std::ofstream out("/tmp/testGrid.m");
  ASSERT_TRUE(sim1.dumpMatlabGrid(out,"TheMesh",{"Boundary","Edge2"}));
  out.close();

  // Read the matlab file into a new SIM and compare the models
  SIM2D_matlab sim2("/tmp/testGrid.m");
  EXPECT_EQ(sim1.getNoNodes(),sim2.getNoNodes());
  EXPECT_EQ(sim1.getNoElms(),sim2.getNoElms());
  ASMbase* pch1 = sim1.getPatch(1);
  ASMbase* pch2 = sim2.getPatch(1);
  ASSERT_TRUE(pch1 != nullptr);
  ASSERT_TRUE(pch2 != nullptr);
  IntVec b1, b2 = pch2->getNodeSet(pch2->getNodeSetIdx("Boundary"));
  for (int edge = 1; edge <= 4; edge++)
    pch1->getBoundaryNodes(edge,b1);
  std::set<int> b1set;
  for (int n : b1) b1set.insert(n);
  EXPECT_EQ(b1set.size(),b2.size());
  EXPECT_TRUE(IntVec(b1set.begin(),b1set.end()) == b2);
  IntVec e1, e2 = pch2->getNodeSet(pch2->getNodeSetIdx("Edge2"));
  pch1->getBoundaryNodes(2,e1);
  EXPECT_EQ(e1.size(),e2.size());
  EXPECT_TRUE(e1 == e2);
}
