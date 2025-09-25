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

#include "SIM1D.h"
#include "SIM2D.h"
#include "ASMs1D.h"
#include "ASMs2D.h"
#include "ModelGenerator.h"
#include "tinyxml2.h"
#include <fstream>

#include <catch2/catch_test_macros.hpp>


class SIM1D_default : public SIM1D
{
public:
  SIM1D_default(int nu)
  {
    // Create a unit grid with (nu+1) linear elements
    this->createDefaultModel();
    ASMs1D* pch1 = static_cast<ASMs1D*>(myModel.front());
    REQUIRE(pch1->uniformRefine(nu));
    REQUIRE(this->createFEMmodel());

    // Create topological boundary entities (vertices)
    tinyxml2::XMLDocument doc;
    doc.Parse("<geometry sets='true'/>");
    DefaultGeometry1D(doc.RootElement()).createTopologySets(*this);
    REQUIRE(myEntitys.size() == 4);
  }
  virtual ~SIM1D_default() {}
};


class SIM2D_default : public SIM2D
{
public:
  SIM2D_default(int nu, int nv)
  {
    // Create a bi-unit square grid with (nu+1)x(nv+1) linear elements
    this->createDefaultModel();
    ASMs2D* pch1 = static_cast<ASMs2D*>(myModel.front());
    pch1->uniformRefine(0,nu);
    pch1->uniformRefine(1,nv);
    REQUIRE(this->createFEMmodel());

    // Create topological boundary entities (vertices and edges)
    tinyxml2::XMLDocument doc;
    doc.Parse("<geometry sets='true'/>");
    DefaultGeometry2D(doc.RootElement()).createTopologySets(*this);
    REQUIRE(myEntitys.size() == 10);
  }
  virtual ~SIM2D_default() {}
};


class SIM1D_matlab : public SIM1D
{
public:
  SIM1D_matlab(const char* inputFile)
  {
    std::string xml("<geometry><patchfile type='matlab'>");
    xml.append(inputFile);
    xml.append("</patchfile></geometry>");

    // Read the FE model from the provided Matlab file
    REQUIRE(this->loadXML(xml.c_str()));
    REQUIRE(this->createFEMmodel());
  }
  virtual ~SIM1D_matlab() {}
};


class SIM2D_matlab : public SIM2D
{
public:
  SIM2D_matlab(const char* inputFile)
  {
    std::string xml("<geometry><patchfile type='matlab'>");
    xml.append(inputFile);
    xml.append("</patchfile></geometry>");

    // Read the FE model from the provided Matlab file
    REQUIRE(this->loadXML(xml.c_str()));
    REQUIRE(this->createFEMmodel());
  }
  virtual ~SIM2D_matlab() {}
};


TEST_CASE("TestMatlabPatch.IO")
{
  // Create a 5x4 element mesh and write it to a matlab file
  SIM2D_default sim1(4,3);
  std::ofstream out("/tmp/testGrid2.m");
  REQUIRE(sim1.dumpMatlabGrid(out,"TheMesh",{"Boundary","Edge2"}));
  out.close();

  // Read the matlab file into a new SIM and compare the models
  SIM2D_matlab sim2("/tmp/testGrid2.m");
  REQUIRE(sim1.getNoNodes() == sim2.getNoNodes());
  REQUIRE(sim1.getNoElms() == sim2.getNoElms());
  ASMbase* pch1 = sim1.getPatch(1);
  ASMbase* pch2 = sim2.getPatch(1);
  REQUIRE(pch1 != nullptr);
  REQUIRE(pch2 != nullptr);
  int idx1 = pch2->getNodeSetIdx("Boundary");
  int idx2 = pch2->getNodeSetIdx("Edge2");
  int idx3 = pch2->parseNodeSet("ACorner","1");
  REQUIRE(idx1 == 1);
  REQUIRE(idx2 == 2);
  REQUIRE(idx3 == 3);

  IntVec b1, b2 = pch2->getNodeSet(idx1);
  for (int edge = 1; edge <= 4; edge++)
    pch1->getBoundaryNodes(edge,b1);
  std::set<int> b1set;
  for (int n : b1) b1set.insert(n);
  REQUIRE(b1set.size() == b2.size());
  REQUIRE(IntVec(b1set.begin(),b1set.end()) == b2);
  IntVec e1, e2 = pch2->getNodeSet(idx2);
  pch1->getBoundaryNodes(2,e1);
  REQUIRE(e1.size() == e2.size());
  REQUIRE(e1 == e2);
  IntVec corner = pch2->getNodeSet(idx3);
  REQUIRE(corner.size() == 1);
  REQUIRE(corner.front() == 1);

  // Create a 6-element 1D mesh and write it to a matlab file
  SIM1D_default sim3(5);
  out.open("/tmp/testGrid1.m");
  REQUIRE(sim3.dumpMatlabGrid(out,"TheMesh",{"Boundary","Vertex2"}));
  out.close();

  // Read the matlab file into a new SIM and compare the models
  SIM1D_matlab sim4("/tmp/testGrid1.m");
  REQUIRE(sim3.getNoNodes() == sim3.getNoNodes());
  REQUIRE(sim3.getNoElms() == sim3.getNoElms());
  pch1 = sim3.getPatch(1);
  pch2 = sim4.getPatch(1);
  REQUIRE(pch1 != nullptr);
  REQUIRE(pch2 != nullptr);
  idx1 = pch2->getNodeSetIdx("Boundary");
  idx2 = pch2->getNodeSetIdx("Vertex2");
  REQUIRE(idx1 == 1);
  REQUIRE(idx2 == 2);

  b1.clear();
  b2 = pch2->getNodeSet(idx1);
  for (int vert = 1; vert <= 2; vert++)
    pch1->getBoundaryNodes(vert,b1);
  b1set.clear();
  for (int n : b1) b1set.insert(n);
  REQUIRE(b1set.size() == b2.size());
  REQUIRE(IntVec(b1set.begin(),b1set.end()) == b2);
  IntVec v1, v2 = pch2->getNodeSet(idx2);
  pch1->getBoundaryNodes(2,v1);
  REQUIRE(v1.size() == v2.size());
  REQUIRE(v1 == v2);
}
