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
#include "SIM2D.h"
#include "LRSpline/LRSplineSurface.h"

#include "gtest/gtest.h"


struct EdgeTest
{
  int edge;
  int edgeIdx;
  std::array<int,2> c1;
  std::array<int,2> c2;
};


class TestASMu2D : public testing::Test,
                   public testing::WithParamInterface<EdgeTest>
{
};


static ASMu2D* getPatch (SIMinput& sim)
{
  sim.opt.discretization = ASM::LRSpline;
  EXPECT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes.xinp"));
  EXPECT_TRUE(sim.createFEMmodel());
  return static_cast<ASMu2D*>(sim.getPatch(1));
}


TEST_P(TestASMu2D, ConstrainEdge)
{
  SIM2D sim(1);
  ASMu2D* pch = getPatch(sim);
  ASSERT_TRUE(pch != nullptr);
  pch->constrainEdge(GetParam().edgeIdx, false, 1, 1, 1);
  std::vector<int> glbNodes;
  pch->getBoundaryNodes(GetParam().edge, glbNodes, 1, 1, 0);
  for (int& it : glbNodes)
    ASSERT_TRUE(pch->findMPC(it, 1) != nullptr);
}


TEST_P(TestASMu2D, ConstrainEdgeOpen)
{
  SIM2D sim(1);
  ASMu2D* pch = getPatch(sim);
  ASSERT_TRUE(pch != nullptr);
  pch->constrainEdge(GetParam().edgeIdx, true, 1, 1, 1);
  std::vector<int> glbNodes;
  pch->getBoundaryNodes(GetParam().edge, glbNodes, 1, 1, 0);
  int crn = pch->getCorner(GetParam().c1[0], GetParam().c1[1], 1);
  ASSERT_TRUE(pch->findMPC(crn, 1) == nullptr);
  glbNodes.erase(std::find(glbNodes.begin(), glbNodes.end(), crn));
  crn = pch->getCorner(GetParam().c2[0], GetParam().c2[1], 1);
  ASSERT_TRUE(pch->findMPC(crn, 1) == nullptr);
  glbNodes.erase(std::find(glbNodes.begin(), glbNodes.end(), crn));
  for (int& it : glbNodes)
    ASSERT_TRUE(pch->findMPC(it, 1) != nullptr);
}


static const std::vector<EdgeTest> edgeTestData =
        {{{1, -1, {-1, -1}, {-1 , 1}},
          {2,  1, { 1, -1}, { 1 , 1}},
          {3, -2, {-1, -1}, { 1, -1}},
          {4,  2, {-1,  1}, { 1,  1}}}};


INSTANTIATE_TEST_CASE_P(TestASMu2D, TestASMu2D,
                        testing::ValuesIn(edgeTestData));


TEST(TestASMu2D, InterfaceChecker)
{
  SIM2D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  sim.createDefaultModel();
  ASMu2D* pch = static_cast<ASMu2D*>(sim.getPatch(1));
  pch->raiseOrder(1,1);
  pch->uniformRefine(0, 3);
  pch->uniformRefine(1, 3);
  LR::LRSplineSurface* lr = pch->getSurface();
  lr->insert_const_u_edge(0.5,   0, 1, 2);
  lr->insert_const_v_edge(0.125, 0.5, 1, 2);
  lr->insert_const_v_edge(0.25,  0.5, 1, 2);
  lr->insert_const_v_edge(0.875, 0.5, 1, 1);
  EXPECT_TRUE(sim.createFEMmodel());

  static std::vector<int> elem_flags = {{10, 11, 11,  9,
                                         14, 15, 15, 13,
                                         14, 15, 15, 13,
                                          6,  7, 15, 13,
                                         15, 13,
                                          7,  5}};

  static std::vector<std::array<std::vector<double>,4>> elem_pts =
  {{{{    {},        {0.25},     {}, {0.25}}},
   {{ {0.25}, {0.125, 0.25},     {},  {0.5}}},
   {{{0.125},       {0.125},     {}, {0.75}}},
   {{{0.125},            {},     {},  {1.0}}},
   {{     {},         {0.5}, {0.25}, {0.25}}},
   {{  {0.5},         {0.5},  {0.5},  {0.5}}},
   {{  {0.5},         {0.5}, {0.75}, {0.75}}},
   {{  {0.5},         {0.5},  {1.0},  {1.0}}},
   {{     {},        {0.75}, {0.25}, {0.25}}},
   {{ {0.75},        {0.75},  {0.5},  {0.5}}},
   {{ {0.75},        {0.75}, {0.75}, {0.75}}},
   {{ {0.75},            {},  {1.0},  {1.0}}},
   {{     {},         {1.0}, {0.25},     {}}},
   {{  {1.0},  {0.875, 1.0},  {0.5},  {0.5}}},
   {{{0.875},       {0.875}, {0.75}, {0.75}}},
   {{{0.875},            {},  {1.0},  {1.0}}},
   {{ {0.25},        {0.25}, {0.75}, {0.75}}},
   {{ {0.25},            {},  {1.0},  {1.0}}},
   {{  {1.0},         {1.0}, {0.75},     {}}},
   {{  {1.0},            {},  {1.0},     {}}}}};

  EXPECT_EQ(elem_flags.size(), pch->getNoElms());
  EXPECT_EQ(elem_pts.size(), pch->getNoElms());
  ASMu2D::InterfaceChecker iChk(*pch);
  for (size_t i = 0; i < pch->getNoElms(); ++i) {
    EXPECT_EQ(iChk.hasContribution(i+1), elem_flags[i]);
    for (size_t edge = 1; edge <= 4; ++edge) {
      if (elem_flags[i] & (1 << (edge-1))) {
        auto val = iChk.getIntersections(i+1, edge);
        EXPECT_EQ(val.size(), elem_pts[i][edge-1].size());
        for (size_t j = 0; j < val.size(); ++j)
          EXPECT_FLOAT_EQ(val[j], elem_pts[i][edge-1][j]);
      }
    }
  }
}
