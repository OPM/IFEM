//==============================================================================
//!
//! \file TestASMs1D.C
//!
//! \date Nov 17 2018
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Tests for structured 1D spline FE models.
//!
//==============================================================================

#include "ASMs1D.h"
#include "Vec3.h"
#include "Vec3Oper.h"
#include <sstream>

#include "Catch2Support.h"


class ASMLine : public ASMs1D
{
public:
  explicit ASMLine(int extraKnots = 0, int extraOrder = 0, bool init = true)
  {
    std::stringstream geo("100 1 0 0\n2 0\n2 2\n0 0 1 1\n0 0\n1 0\n");
    REQUIRE(this->read(geo));
    REQUIRE(this->raiseOrder(extraOrder));
    if (extraKnots > 0)
      REQUIRE(this->uniformRefine(extraKnots));
    if (init)
      REQUIRE(this->generateFEMTopology());
  }

  void shiftElmNumbers(int shift)
  {
    for (int& e : myMLGE)
      e += (e == -1 ? 0 : shift);
  }

  virtual ~ASMLine() {}
};


class ASMArch : public ASMs1D
{
public:
  explicit ASMArch(int extraOrder = 0, int extraKnots = 0)
  {
    std::stringstream geo("100 1 0 0\n"
                          "3 1\n3 3\n"
                          "0 0 0 "
                          "0.2561196234760725 0.2561196234760725 0.2561196234760725\n"
                          "17 0 0 1\n"
                          "0 2.171071412473507 0 0.9918115419161011\n"
                          "-17 0 0 1\n");
    REQUIRE(this->read(geo));
    REQUIRE(this->raiseOrder(extraOrder));
    REQUIRE(this->uniformRefine(extraKnots));
    REQUIRE(this->generateFEMTopology());
  }
  virtual ~ASMArch() {}
};


TEST_CASE("TestASMs1D.FindClosestNode")
{
  ASMbase::resetNumbering();

  ASMLine pch1(7);
  std::cout <<"\nTesting a linear patch: ";
  REQUIRE(pch1.write(std::cout));
  REQUIRE(pch1.getNoNodes() == 9);

  Vec3 X;
  double x0  = 0.001;
  double x1  = 0.999;
  double xin = 0.5, u;
  int ipoint = pch1.evalPoint(&xin,&u,X);
  std::cout <<"u = "<< u <<" --> X = "<< X;
  if (ipoint) std::cout <<" (inod = "<< ipoint <<")";
  std::cout << std::endl;

  std::pair<size_t,double> cnod = pch1.findClosestNode(X);
  std::cout <<"Closest node "<< cnod.first <<" distance="<< cnod.second << std::endl;
  REQUIRE(cnod.first == 5);

  pch1.evalPoint(&x0,&u,X);
  REQUIRE(pch1.findClosestNode(X).first == 1);

  pch1.evalPoint(&x1,&u,X);
  REQUIRE(pch1.findClosestNode(X).first == 9);

  int p = 2;
  const char* order[] = { "quadratic","cubic","quartic","quintic","sextic","septic", nullptr };
  for (const char** q = order; *q; p++, q++)
  {
    ASMArch pch(p-2,8-p);
    std::cout <<"\nTesting a "<< *q <<" patch: ";
    REQUIRE(pch.write(std::cout));
    REQUIRE(pch.getNoNodes() == 9);

    ipoint = pch.evalPoint(&xin,&u,X);
    std::cout <<"u = "<< u <<" --> X = "<< X;
    if (ipoint) std::cout <<" (inod = "<< ipoint <<")";
    std::cout << std::endl;

    cnod = pch.findClosestNode(X);
    std::cout <<"Closest node "<< cnod.first <<" distance="<< cnod.second << std::endl;
    REQUIRE(cnod.first == 5);

    pch.evalPoint(&x0,&u,X);
    REQUIRE(pch.findClosestNode(X).first == 1);

    pch.evalPoint(&x1,&u,X);
    REQUIRE(pch.findClosestNode(X).first == 9);
  }
}


TEST_CASE("TestASMs1D.FindElement")
{
  ASMbase::resetNumbering();

  ASMLine pch1(9);
  std::cout <<"\nTesting a linear patch: ";
  REQUIRE(pch1.write(std::cout));
  REQUIRE(pch1.getNoNodes() == 11);
  REQUIRE(pch1.findElement({-0.03,0.0,0.0}).first == 1);
  REQUIRE(pch1.findElement({ 0.02,0.0,0.0}).first == 1);
  REQUIRE(pch1.findElement({ 0.43,0.0,0.0}).first == 5);
  REQUIRE(pch1.findElement({ 0.96,0.0,0.0}).first == 10);
  REQUIRE(pch1.findElement({ 1.16,0.0,0.0}).first == 10);

  ASMLine pch2(9,1);
  std::cout <<"\nTesting a quadratic patch: ";
  REQUIRE(pch2.write(std::cout));
  REQUIRE(pch2.getNoNodes() == 12);
  REQUIRE(pch2.findElement({-0.03,0.0,0.0}).first == 1);
  REQUIRE(pch2.findElement({ 0.02,0.0,0.0}).first == 1);
  REQUIRE(pch2.findElement({ 0.43,0.0,0.0}).first == 5);
  REQUIRE(pch2.findElement({ 0.96,0.0,0.0}).first == 10);
  REQUIRE(pch2.findElement({ 1.16,0.0,0.0}).first == 10);
}


TEST_CASE("TestASMs1D.ElementConnectivities")
{
  ASMbase::resetNumbering();
  const std::array<IntVec,3> refLoc = {{{-1, 1}, {0, 2}, {1, -1}}};

  ASMLine pch1(2);
  const size_t nel = pch1.getNoElms();
  pch1.shiftElmNumbers(nel);
  IntMat neighGlb(nel*2), neighLoc(nel);
  pch1.getElmConnectivities(neighGlb);
  pch1.getElmConnectivities(neighLoc, ASM::GEOMETRY_BASIS);
  REQUIRE(neighLoc.size() == nel);
  REQUIRE(neighGlb.size() == 2*nel);
  for (size_t n = 0; n < neighLoc.size(); ++n) {
    REQUIRE(neighLoc[n].size() == refLoc[n].size());
    REQUIRE(neighGlb[n+nel].size() == refLoc[n].size());
    for (size_t i = 0; i < neighLoc[n].size(); ++i) {
      REQUIRE(neighLoc[n][i] == refLoc[n][i]);
      REQUIRE(neighGlb[n+nel][i] == (refLoc[n][i] > -1 ?
                                     refLoc[n][i] + static_cast<int>(nel) : -1));
    }
  }
}


TEST_CASE("TestASMs1D.BoundaryElements")
{
  ASMbase::resetNumbering();

  ASMLine pch1(2);
  for (size_t i = 1; i <= 2; ++i) {
    IntVec nodes;
    pch1.getBoundaryElms(i, nodes);
    REQUIRE(nodes.size() == 1);
    REQUIRE(nodes.front() == (i == 1 ? 0 : 2));
  }
}


TEST_CASE("TestASMs1D.ElmNodes")
{
  ASMbase::resetNumbering();

  ASMLine pch1(0, 0, false);
  pch1.createProjectionBasis(true);
  pch1.raiseOrder(1);
  pch1.uniformRefine(2);
  pch1.createProjectionBasis(false);
  pch1.uniformRefine(2);
  REQUIRE(pch1.generateFEMTopology());
  const IntMat mnpc = pch1.getElmNodes(1);

  const auto ref = std::array{
      std::array{0,1},
      std::array{1,2},
      std::array{2,3},
  };
  REQUIRE(mnpc.size() == ref.size());
  for (size_t i = 0; i < mnpc.size(); ++i) {
    REQUIRE(mnpc[i].size() == ref[i].size());
    for (size_t j = 0; j < mnpc[i].size(); ++j)
      REQUIRE(mnpc[i][j] == ref[i][j]);
  }

  const auto ref_proj = std::array{
      std::array{0,1,2},
      std::array{1,2,3},
      std::array{2,3,4},
  };
  const IntMat mnpc_proj = pch1.getElmNodes(ASM::PROJECTION_BASIS);
  REQUIRE(mnpc_proj.size() == ref_proj.size());
  for (size_t i = 0; i < mnpc_proj.size(); ++i) {
    REQUIRE(mnpc_proj[i].size() == ref_proj[i].size());
    for (size_t j = 0; j < mnpc_proj[i].size(); ++j)
      REQUIRE(mnpc_proj[i][j] == ref_proj[i][j]);
  }
}
