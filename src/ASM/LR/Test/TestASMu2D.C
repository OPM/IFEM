//==============================================================================
//!
//! \file TestASMu2D.C
//!
//! \date Jul 14 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for driver for assembly of unstructured 2D spline FE models.
//!
//==============================================================================

#include "ASMu2D.h"
#include "IFEM.h"
#include "SIM2D.h"

#include "gtest/gtest.h"


struct EdgeTest {
  int edge;
  int edgeIdx;
  std::array<int,2> c1;
  std::array<int,2> c2;
};


class TestASMu2D : public testing::Test,
                   public testing::WithParamInterface<EdgeTest>
{
};


TEST(TestASMu2D, BoundaryNodesE1)
{
  SIM2D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes.xinp"));
  sim.preprocess();

  int bcode = sim.getUniquePropertyCode("Edge1",0);
  std::vector<int> vec;
  sim.getBoundaryNodes(bcode, vec);
  ASSERT_EQ(vec.size(), 4U);
  for (int i = 0; i < 4; ++i)
    ASSERT_EQ(vec[i], 4*i+1);
}


TEST(TestASMu2D, BoundaryNodesE2)
{
  SIM2D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes.xinp"));
  sim.preprocess();

  int bcode = sim.getUniquePropertyCode("Edge2",0);
  std::vector<int> vec;
  sim.getBoundaryNodes(bcode, vec);
  ASSERT_EQ(vec.size(), 4U);
  for (int i = 0; i < 4; ++i)
    ASSERT_EQ(vec[i], 2+4*i);
}


TEST(TestASMu2D, BoundaryNodesE3)
{
  SIM2D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes.xinp"));
  sim.preprocess();

  int bcode = sim.getUniquePropertyCode("Edge3",0);
  std::vector<int> vec;
  sim.getBoundaryNodes(bcode, vec);
  ASSERT_EQ(vec.size(), 4U);
  for (int i = 0; i < 4; ++i)
    ASSERT_EQ(vec[i], i+1);
}


TEST(TestASMu2D, BoundaryNodesE4)
{
  SIM2D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes.xinp"));
  sim.preprocess();

  int bcode = sim.getUniquePropertyCode("Edge4",0);
  std::vector<int> vec;
  sim.getBoundaryNodes(bcode, vec);
  ASSERT_EQ(vec.size(), 4U);
  for (int i = 0; i < 4; ++i)
    ASSERT_EQ(vec[i], 5+i);
}


TEST_P(TestASMu2D, ConstrainEdge)
{
  SIM2D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes.xinp"));
  sim.preprocess();
  ASMu2D* pch = static_cast<ASMu2D*>(sim.getPatch(1));
  pch->constrainEdge(GetParam().edgeIdx, false, 1, 1, 1);
  std::vector<int> glbNodes;
  pch->getBoundaryNodes(GetParam().edge, glbNodes, 1);
  for (auto& it : glbNodes)
    ASSERT_TRUE(pch->findMPC(it, 1) != nullptr);
}


TEST_P(TestASMu2D, ConstrainEdgeOpen)
{
  SIM2D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes.xinp"));
  sim.preprocess();
  ASMu2D* pch = static_cast<ASMu2D*>(sim.getPatch(1));
  pch->constrainEdge(GetParam().edgeIdx, true, 1, 1, 1);
  std::vector<int> glbNodes;
  pch->getBoundaryNodes(GetParam().edge, glbNodes, 1);
  int crn = pch->getCorner(GetParam().c1[0], GetParam().c1[1], 1);
  ASSERT_TRUE(pch->findMPC(crn, 1) == nullptr);
  glbNodes.erase(std::find(glbNodes.begin(), glbNodes.end(), crn));
  crn = pch->getCorner(GetParam().c2[0], GetParam().c2[1], 1);
  ASSERT_TRUE(pch->findMPC(crn, 1) == nullptr);
  glbNodes.erase(std::find(glbNodes.begin(), glbNodes.end(), crn));
  for (auto& it : glbNodes)
    ASSERT_TRUE(pch->findMPC(it, 1) != nullptr);
}


static const std::vector<EdgeTest> edgeTestData =
        {{1, -1, {-1, -1}, {-1 , 1}},
         {2,  1, { 1, -1}, { 1 , 1}},
         {3, -2, {-1, -1}, { 1, -1}},
         {4,  2, {-1,  1}, { 1,  1}}};



INSTANTIATE_TEST_CASE_P(TestASMu2D, TestASMu2D,
                        testing::ValuesIn(edgeTestData));
