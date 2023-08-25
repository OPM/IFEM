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

#include "ASMuSquare.h"
#include "ASMs2D.h"
#include "GaussQuadrature.h"
#include "SIM2D.h"
#include "LRSpline/LRSplineSurface.h"

#include "gtest/gtest.h"
#include <fstream>
#include <numeric>


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


class TestuSIM2D : public SIM2D
{
public:
  TestuSIM2D() : SIM2D(1)
  {
    opt.discretization = ASM::LRSpline;
    EXPECT_TRUE(this->read("src/ASM/LR/Test/refdata/boundary_nodes.xinp"));
    EXPECT_TRUE(this->createFEMmodel());
  }
  virtual ~TestuSIM2D() {}
};


TEST_P(TestASMu2D, ConstrainEdge)
{
  TestuSIM2D sim;
  ASMu2D* pch = static_cast<ASMu2D*>(sim.getPatch(1));
  ASSERT_TRUE(pch != nullptr);

  pch->constrainEdge(GetParam().edgeIdx, false, 1, 1, 1);

  std::vector<int> glbNodes;
  pch->getBoundaryNodes(GetParam().edge, glbNodes, 1);

  for (int node : glbNodes)
    EXPECT_TRUE(pch->findMPC(node,1) != nullptr);
}


TEST_P(TestASMu2D, ConstrainEdgeOpen)
{
  TestuSIM2D sim;
  ASMu2D* pch = static_cast<ASMu2D*>(sim.getPatch(1));
  ASSERT_TRUE(pch != nullptr);

  pch->constrainEdge(GetParam().edgeIdx, true, 1, 1, 1);

  std::vector<int> glbNodes;
  pch->getBoundaryNodes(GetParam().edge, glbNodes, 1);

  int crn = pch->getCorner(GetParam().c1[0], GetParam().c1[1], 1);
  EXPECT_TRUE(pch->findMPC(crn,1) == nullptr);
  glbNodes.erase(std::find(glbNodes.begin(), glbNodes.end(), crn));
  crn = pch->getCorner(GetParam().c2[0], GetParam().c2[1], 1);
  EXPECT_TRUE(pch->findMPC(crn,1) == nullptr);
  glbNodes.erase(std::find(glbNodes.begin(), glbNodes.end(), crn));

  for (int node : glbNodes)
    EXPECT_TRUE(pch->findMPC(node,1) != nullptr);
}


static const std::vector<EdgeTest> edgeTestData =
        {{{1, -1, {{-1, -1}}, {{-1 , 1}}},
          {2,  1, {{ 1, -1}}, {{ 1 , 1}}},
          {3, -2, {{-1, -1}}, {{ 1, -1}}},
          {4,  2, {{-1,  1}}, {{ 1,  1}}}}};


INSTANTIATE_TEST_CASE_P(TestASMu2D, TestASMu2D,
                        testing::ValuesIn(edgeTestData));


TEST(TestASMu2D, InterfaceChecker)
{
  SIM2D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASMu2D* pch = static_cast<ASMu2D*>(sim.createDefaultModel());
  pch->raiseOrder(1,1);
  pch->uniformRefine(0, 3);
  pch->uniformRefine(1, 3);
  LR::LRSplineSurface* lr = pch->getBasis();
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

  static std::vector<std::array<int,4>> line_cont = {{{{-1,  1, -1,  1}},
                                                      {{ 1,  0, -1,  1}},
                                                      {{ 0,  1, -1,  0}},
                                                      {{ 1, -1, -1,  0}},
                                                      {{-1,  1,  1,  1}},
                                                      {{ 1,  0,  1,  1}},
                                                      {{ 0,  1,  0,  1}},
                                                      {{ 1, -1,  0,  1}},
                                                      {{-1,  1,  1,  1}},
                                                      {{ 1,  0,  1,  1}},
                                                      {{ 0,  1,  1,  1}},
                                                      {{ 1, -1,  1,  1}},
                                                      {{-1,  1,  1, -1}},
                                                      {{ 1,  0,  1, -1}},
                                                      {{ 0,  1,  1,  1}},
                                                      {{ 1, -1,  1,  1}},
                                                      {{ 0,  1,  0,  0}},
                                                      {{ 1, -1,  0,  0}},
                                                      {{ 0,  1,  1, -1}},
                                                      {{ 1, -1,  1, -1}}}};


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
        int continuity;
        auto val = iChk.getIntersections(i+1, edge, &continuity);
        EXPECT_EQ(line_cont[i][edge-1], continuity);
        EXPECT_EQ(val.size(), elem_pts[i][edge-1].size());
        for (size_t j = 0; j < val.size(); ++j)
          EXPECT_FLOAT_EQ(val[j], elem_pts[i][edge-1][j]);
      }
    }
  }
}


TEST(TestASMu2D, TransferGaussPtVars)
{
  SIM2D sim(1);
  sim.opt.discretization = ASM::LRSpline;

  ASMu2D* pch = static_cast<ASMu2D*>(sim.createDefaultModel());
  LR::LRSplineSurface* lr = pch->getBasis(1);
  lr->generateIDs();

  RealArray oldAr(9), newAr;
  const double* xi = GaussQuadrature::getCoord(3);
  size_t id[2];
  for (size_t idx = 0; idx < 2; ++idx) {
    SIM2D simNew(1);
    simNew.opt.discretization = ASM::LRSpline;
    ASMu2D* pchNew = static_cast<ASMu2D*>(simNew.createDefaultModel());
    pchNew->uniformRefine(idx, 1);
    LR::LRSplineSurface* lrNew = pchNew->getBasis(1);
    ASSERT_TRUE(lrNew != nullptr);
    lrNew->generateIDs();
    for (id[1] = 0; id[1] < 3; ++id[1])
      for (id[0] = 0; id[0] < 3; ++id[0])
        oldAr[id[0]+id[1]*3] = (1.0 + xi[id[idx]]) / 2.0;
    pchNew->transferGaussPtVars(lr, oldAr, newAr, 3);
    size_t k = 0;
    for (size_t iEl = 0; iEl < 2; ++iEl)
      for (id[1] = 0; id[1] < 3; ++id[1])
        for (id[0] = 0; id[0] < 3; ++id[0], ++k)
          EXPECT_FLOAT_EQ(newAr[k], 0.5*iEl + 0.5*(xi[id[idx]] + 1.0) / 2.0);
  }
}


TEST(TestASMu2D, TransferGaussPtVarsN)
{
  SIM2D sim(1), sim2(1);
  sim.opt.discretization = sim2.opt.discretization = ASM::LRSpline;

  ASMu2D* pch = static_cast<ASMu2D*>(sim.createDefaultModel());
  LR::LRSplineSurface* lr = pch->getBasis(1);
  lr->generateIDs();

  ASMu2D* pchNew = static_cast<ASMu2D*>(sim2.createDefaultModel());
  pchNew->uniformRefine(0, 1);
  LR::LRSplineSurface* lrNew = pchNew->getBasis(1);
  ASSERT_TRUE(lrNew != nullptr);
  lrNew->generateIDs();

  RealArray oldAr(9), newAr;
  std::iota(oldAr.begin(), oldAr.end(), 1);

  pchNew->transferGaussPtVarsN(lr, oldAr, newAr, 3);
  static RealArray refAr = {{1.0, 1.0, 2.0,
                             4.0, 4.0, 5.0,
                             7.0, 7.0, 8.0,
                             2.0, 3.0, 3.0,
                             5.0, 6.0, 6.0,
                             8.0, 9.0, 9.0}};
  EXPECT_EQ(refAr.size(), newAr.size());
  for (size_t i = 0; i < refAr.size(); ++i)
    EXPECT_FLOAT_EQ(refAr[i], newAr[i]);
}


TEST(TestASMu2D, Connect)
{
  SIM2D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient0.xinp"));
  ASSERT_TRUE(sim.createFEMmodel());

  SIM2D sim2(1);
  sim2.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim2.read("src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient1.xinp"));
  ASSERT_TRUE(sim2.createFEMmodel());
}


TEST(TestASMu2D, ElementConnectivities)
{
  ASMuSquare pch1;
  ASMbase::resetNumbering();
  ASSERT_TRUE(pch1.uniformRefine(0,1));
  ASSERT_TRUE(pch1.uniformRefine(1,1));
  ASSERT_TRUE(pch1.generateFEMTopology());
  pch1.getBasis(1)->generateIDs();
  IntMat neigh(4);
  pch1.getElmConnectivities(neigh);
  const std::array<std::vector<int>,4> ref = {{{1, 2},
                                               {0, 3},
                                               {3, 0},
                                               {2, 1}}};
  ASSERT_EQ(neigh.size(), 4U);
  for (size_t n = 0; n < neigh.size(); ++n) {
    ASSERT_EQ(neigh[n].size(), ref[n].size());
    for (size_t i = 0; i < neigh[n].size(); ++i)
      EXPECT_EQ(neigh[n][i], ref[n][i]);
  }
}

TEST(TestASMu2D, EvalPointNurbs)
{
  ASMu2D u2Dpch(2, 1);
  std::ifstream g2file("src/ASM/LR/Test/refdata/hole2D.g2");
  std::ifstream g2file2("src/ASM/LR/Test/refdata/hole2D.g2");
  ASMs2D s2Dpch(2,1);
  ASSERT_TRUE(u2Dpch.read(g2file));
  ASSERT_TRUE(s2Dpch.read(g2file2));
  ASSERT_TRUE(u2Dpch.generateFEMTopology());
  ASSERT_TRUE(s2Dpch.generateFEMTopology());

  double xi[2] = {0.1, 0.1};
  double param1[2], param2[2];
  Vec3 x1, x2;
  s2Dpch.evalPoint(xi,param1,x1);
  u2Dpch.evalPoint(xi,param2,x2);
  EXPECT_FLOAT_EQ(param1[0], param2[0]);
  EXPECT_FLOAT_EQ(param1[1], param2[1]);
  EXPECT_FLOAT_EQ(x1[0], x2[0]);
  EXPECT_FLOAT_EQ(x1[1], x2[1]);
}


TEST(TestASMu2D, Write)
{
  ASMbase::resetNumbering();
  ASMuSquare pch1;
  EXPECT_TRUE(pch1.generateFEMTopology());

  std::stringstream str;
  EXPECT_TRUE(pch1.write(str, 1));
  EXPECT_EQ(str.str(), ASMuSquare::square);

  EXPECT_FALSE(pch1.write(str, 2));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::GEOMETRY_BASIS));
  EXPECT_EQ(str.str(), ASMuSquare::square);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::PROJECTION_BASIS));
  EXPECT_EQ(str.str(), ASMuSquare::square);

  EXPECT_FALSE(pch1.write(str, ASM::PROJECTION_BASIS_2));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::REFINEMENT_BASIS));
  EXPECT_EQ(str.str(), ASMuSquare::square);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::INTEGRATION_BASIS));
  EXPECT_EQ(str.str(), ASMuSquare::square);
}
