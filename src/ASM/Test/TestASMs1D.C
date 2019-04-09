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
  ASMLine()
  {
    std::stringstream geo("100 1 0 0\n2 0\n2 2\n0 0 1 1\n0 0\n1 0\n");
    EXPECT_TRUE(this->read(geo));
  }
  virtual ~ASMLine() {}
};


class ASMArch : public ASMs1D
{
public:
  ASMArch()
  {
    std::stringstream geo("100 1 0 0\n"
                          "3 1\n3 3\n"
                          "0 0 0 "
                          "0.2561196234760725 0.2561196234760725 0.2561196234760725\n"
                          "17 0 0 1\n"
                          "0 2.171071412473507 0 0.9918115419161011\n"
                          "-17 0 0 1\n");
    EXPECT_TRUE(this->read(geo));
  }
  virtual ~ASMArch() {}
};


TEST(TestASMs1D, findClosestNode)
{
  ASMLine pch1;
  std::cout <<"\nTesting a linear patch: ";
  ASSERT_TRUE(pch1.uniformRefine(7));
  EXPECT_TRUE(pch1.write(std::cout));
  ASSERT_TRUE(pch1.generateFEMTopology());
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
    ASMArch pch;
    std::cout <<"\nTesting a "<< *q <<" patch: ";
    ASSERT_TRUE(pch.raiseOrder(p-2));
    ASSERT_TRUE(pch.uniformRefine(8-p));
    EXPECT_TRUE(pch.write(std::cout));
    ASSERT_TRUE(pch.generateFEMTopology());
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


TEST(TestASMs1D, ElementConnectivities)
{
  ASMLine pch1;
  ASSERT_TRUE(pch1.uniformRefine(2));
  ASSERT_TRUE(pch1.generateFEMTopology());
  IntMat neigh(3);
  pch1.getElmConnectivities(neigh);
  const std::array<std::vector<int>,3> ref = {{{-1, 1}, {0, 2}, {1, -1}}};
  ASSERT_EQ(neigh.size(), 3U);
  for (size_t n = 0; n < neigh.size(); ++n) {
    ASSERT_EQ(neigh[n].size(), ref[n].size());
    for (size_t i = 0; i < neigh[n].size(); ++i)
      EXPECT_EQ(neigh[n][i], ref[n][i]);
  }
}


TEST(TestASMs1D, BoundaryElements)
{
  ASMLine pch1;
  ASSERT_TRUE(pch1.uniformRefine(2));
  ASSERT_TRUE(pch1.generateFEMTopology());
  IntMat neigh(3);
  std::array<IntVec,2> n;
  for (size_t i = 1; i <= 2; ++i) {
    pch1.getBoundaryElms(i, 0, n[i-1]);
    ASSERT_EQ(n[i-1].size(), 1U);
    EXPECT_EQ(n[i-1][0], i == 1 ? 0 : 2);
  }
}
