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

#include "gtest/gtest.h"


class ASMLine : public ASMs1D
{
public:
  explicit ASMLine(int extraKnots = 0, int extraOrder = 0, bool init = true)
  {
    std::stringstream geo("100 1 0 0\n2 0\n2 2\n0 0 1 1\n0 0\n1 0\n");
    EXPECT_TRUE(this->read(geo));
    EXPECT_TRUE(this->raiseOrder(extraOrder));
    if (extraKnots > 0)
      EXPECT_TRUE(this->uniformRefine(extraKnots));
    if (init)
      EXPECT_TRUE(this->generateFEMTopology());
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
    EXPECT_TRUE(this->read(geo));
    EXPECT_TRUE(this->raiseOrder(extraOrder));
    EXPECT_TRUE(this->uniformRefine(extraKnots));
    EXPECT_TRUE(this->generateFEMTopology());
  }
  virtual ~ASMArch() {}
};


TEST(TestASMs1D, findClosestNode)
{
  ASMbase::resetNumbering();

  ASMLine pch1(7);
  std::cout <<"\nTesting a linear patch: ";
  EXPECT_TRUE(pch1.write(std::cout));
  EXPECT_EQ  (pch1.getNoNodes(),9u);

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
  EXPECT_EQ(cnod.first,5u);

  pch1.evalPoint(&x0,&u,X);
  EXPECT_EQ(pch1.findClosestNode(X).first,1u);

  pch1.evalPoint(&x1,&u,X);
  EXPECT_EQ(pch1.findClosestNode(X).first,9u);

  int p = 2;
  const char* order[] = { "quadratic","cubic","quartic","quintic","sextic","septic", nullptr };
  for (const char** q = order; *q; p++, q++)
  {
    ASMArch pch(p-2,8-p);
    std::cout <<"\nTesting a "<< *q <<" patch: ";
    EXPECT_TRUE(pch.write(std::cout));
    EXPECT_EQ  (pch.getNoNodes(),9u);

    ipoint = pch.evalPoint(&xin,&u,X);
    std::cout <<"u = "<< u <<" --> X = "<< X;
    if (ipoint) std::cout <<" (inod = "<< ipoint <<")";
    std::cout << std::endl;

    cnod = pch.findClosestNode(X);
    std::cout <<"Closest node "<< cnod.first <<" distance="<< cnod.second << std::endl;
    EXPECT_EQ(cnod.first,5u);

    pch.evalPoint(&x0,&u,X);
    EXPECT_EQ(pch.findClosestNode(X).first,1u);

    pch.evalPoint(&x1,&u,X);
    EXPECT_EQ(pch.findClosestNode(X).first,9u);
  }
}


TEST(TestASMs1D, findElement)
{
  ASMbase::resetNumbering();

  ASMLine pch1(9);
  std::cout <<"\nTesting a linear patch: ";
  EXPECT_TRUE(pch1.write(std::cout));
  EXPECT_EQ(pch1.getNoNodes(),11u);
  EXPECT_EQ(pch1.findElement({-0.03,0.0,0.0}).first,1);
  EXPECT_EQ(pch1.findElement({ 0.02,0.0,0.0}).first,1);
  EXPECT_EQ(pch1.findElement({ 0.43,0.0,0.0}).first,5);
  EXPECT_EQ(pch1.findElement({ 0.96,0.0,0.0}).first,10);
  EXPECT_EQ(pch1.findElement({ 1.16,0.0,0.0}).first,10);

  ASMLine pch2(9,1);
  std::cout <<"\nTesting a quadratic patch: ";
  EXPECT_TRUE(pch2.write(std::cout));
  EXPECT_EQ(pch2.getNoNodes(),12u);
  EXPECT_EQ(pch2.findElement({-0.03,0.0,0.0}).first,1);
  EXPECT_EQ(pch2.findElement({ 0.02,0.0,0.0}).first,1);
  EXPECT_EQ(pch2.findElement({ 0.43,0.0,0.0}).first,5);
  EXPECT_EQ(pch2.findElement({ 0.96,0.0,0.0}).first,10);
  EXPECT_EQ(pch2.findElement({ 1.16,0.0,0.0}).first,10);
}


TEST(TestASMs1D, ElementConnectivities)
{
  ASMbase::resetNumbering();
  const std::array<IntVec,3> ref = {{{-1, 1}, {0, 2}, {1, -1}}};

  ASMLine pch1(2);
  IntMat neigh(3);
  pch1.getElmConnectivities(neigh);
  ASSERT_EQ(neigh.size(), 3U);
  for (size_t n = 0; n < neigh.size(); ++n) {
    ASSERT_EQ(neigh[n].size(), ref[n].size());
    for (size_t i = 0; i < neigh[n].size(); ++i)
      EXPECT_EQ(neigh[n][i], ref[n][i]);
  }
}


TEST(TestASMs1D, BoundaryElements)
{
  ASMbase::resetNumbering();

  ASMLine pch1(2);
  for (size_t i = 1; i <= 2; ++i) {
    IntVec nodes;
    pch1.getBoundaryElms(i, nodes);
    ASSERT_EQ(nodes.size(), 1u);
    EXPECT_EQ(nodes.front(), i == 1 ? 0 : 2);
  }
}


TEST(TestASMs1D, ElmNodes)
{
  ASMbase::resetNumbering();

  ASMLine pch1(0, 0, false);
  pch1.createProjectionBasis(true);
  pch1.raiseOrder(1);
  pch1.uniformRefine(2);
  pch1.createProjectionBasis(false);
  pch1.uniformRefine(2);
  EXPECT_TRUE(pch1.generateFEMTopology());
  const IntMat mnpc = pch1.getElmNodes(1);

  const auto ref = std::array{
      std::array{0,1},
      std::array{1,2},
      std::array{2,3},
  };
  ASSERT_EQ(mnpc.size(), ref.size());
  for (size_t i = 0; i < mnpc.size(); ++i) {
    EXPECT_EQ(mnpc[i].size(), ref[i].size());
    for (size_t j = 0; j < mnpc[i].size(); ++j)
      EXPECT_EQ(mnpc[i][j], ref[i][j]);
  }

  const auto ref_proj = std::array{
      std::array{0,1,2},
      std::array{1,2,3},
      std::array{2,3,4},
  };
  const IntMat mnpc_proj = pch1.getElmNodes(ASM::PROJECTION_BASIS);
  ASSERT_EQ(mnpc_proj.size(), ref_proj.size());
  for (size_t i = 0; i < mnpc_proj.size(); ++i) {
    EXPECT_EQ(mnpc_proj[i].size(), ref_proj[i].size());
    for (size_t j = 0; j < mnpc_proj[i].size(); ++j)
      EXPECT_EQ(mnpc_proj[i][j], ref_proj[i][j]);
  }
}
