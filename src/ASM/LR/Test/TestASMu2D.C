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
#include "ASMuTrap.h"
#include "ASMs2D.h"
#include "GaussQuadrature.h"
#include "SIM2D.h"
#include "Vec3Oper.h"

#include "Catch2Support.h"
#include <catch2/matchers/catch_matchers_range_equals.hpp>

#include <LRSpline/LRSplineSurface.h>

#include <fstream>
#include <numeric>


namespace {

struct EdgeTest
{
  int edge;
  int edgeIdx;
  std::array<int,2> c1;
  std::array<int,2> c2;
};


class TestuSIM2D : public SIM2D
{
public:
  TestuSIM2D() : SIM2D(1)
  {
    opt.discretization = ASM::LRSpline;
    REQUIRE(this->read("src/ASM/LR/Test/refdata/boundary_nodes.xinp"));
    REQUIRE(this->createFEMmodel());
  }
  virtual ~TestuSIM2D() {}
};

}


TEST_CASE("TestASMu2D.ConstrainEdge")
{
  const EdgeTest param = GENERATE(
    EdgeTest{1, -1, {{-1, -1}}, {{-1 , 1}}},
    EdgeTest{2,  1, {{ 1, -1}}, {{ 1 , 1}}},
    EdgeTest{3, -2, {{-1, -1}}, {{ 1, -1}}},
    EdgeTest{4,  2, {{-1,  1}}, {{ 1,  1}}}
  );

  SECTION("Edge" + std::to_string(param.edge)) {
    TestuSIM2D sim;
    ASMu2D* pch = static_cast<ASMu2D*>(sim.getPatch(1));
    REQUIRE(pch != nullptr);

    pch->constrainEdge(param.edgeIdx, false, 1, 1, 1);

    std::vector<int> glbNodes;
    pch->getBoundaryNodes(param.edge, glbNodes, 1);

    for (int node : glbNodes)
      REQUIRE(pch->findMPC(node,1) != nullptr);
  }
}


TEST_CASE("TestASMu2D.ConstrainEdgeOpen")
{
  const EdgeTest param = GENERATE(
    EdgeTest{1, -1, {{-1, -1}}, {{-1 , 1}}},
    EdgeTest{2,  1, {{ 1, -1}}, {{ 1 , 1}}},
    EdgeTest{3, -2, {{-1, -1}}, {{ 1, -1}}},
    EdgeTest{4,  2, {{-1,  1}}, {{ 1,  1}}}
  );

  SECTION("Edge" + std::to_string(param.edge)) {
    TestuSIM2D sim;
    ASMu2D* pch = static_cast<ASMu2D*>(sim.getPatch(1));
    REQUIRE(pch != nullptr);

    pch->constrainEdge(param.edgeIdx, true, 1, 1, 1);

    std::vector<int> glbNodes;
    pch->getBoundaryNodes(param.edge, glbNodes, 1);

    int crn = pch->getCorner(param.c1[0], param.c1[1], 1);
    REQUIRE(pch->findMPC(crn,1) == nullptr);
    glbNodes.erase(std::find(glbNodes.begin(), glbNodes.end(), crn));
    crn = pch->getCorner(param.c2[0], param.c2[1], 1);
    REQUIRE(pch->findMPC(crn,1) == nullptr);
    glbNodes.erase(std::find(glbNodes.begin(), glbNodes.end(), crn));

    for (int node : glbNodes)
      REQUIRE(pch->findMPC(node,1) != nullptr);
  }
}


TEST_CASE("TestASMu2D.InterfaceChecker")
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
  REQUIRE(sim.createFEMmodel());

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

  REQUIRE(elem_flags.size() == pch->getNoElms());
  REQUIRE(elem_pts.size() == pch->getNoElms());
  ASMu2D::InterfaceChecker iChk(*pch);
  for (size_t i = 0; i < pch->getNoElms(); ++i) {
    REQUIRE(iChk.hasContribution(i+1) == elem_flags[i]);
    for (size_t edge = 1; edge <= 4; ++edge) {
      if (elem_flags[i] & (1 << (edge-1))) {
        int continuity;
        auto val = iChk.getIntersections(i+1, edge, &continuity);
        REQUIRE(line_cont[i][edge-1] == continuity);
        REQUIRE(val.size() == elem_pts[i][edge-1].size());
        for (size_t j = 0; j < val.size(); ++j)
          REQUIRE_THAT(val[j], WithinRel(elem_pts[i][edge-1][j]));
      }
    }
  }
}


TEST_CASE("TestASMu2D.TransferGaussPtVars")
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
    REQUIRE(lrNew != nullptr);
    lrNew->generateIDs();
    for (id[1] = 0; id[1] < 3; ++id[1])
      for (id[0] = 0; id[0] < 3; ++id[0])
        oldAr[id[0]+id[1]*3] = (1.0 + xi[id[idx]]) / 2.0;
    pchNew->transferGaussPtVars(lr, oldAr, newAr, 3);
    size_t k = 0;
    for (size_t iEl = 0; iEl < 2; ++iEl)
      for (id[1] = 0; id[1] < 3; ++id[1])
        for (id[0] = 0; id[0] < 3; ++id[0], ++k)
          REQUIRE_THAT(newAr[k], WithinRel(0.5*iEl + 0.5*(xi[id[idx]] + 1.0) / 2.0));
  }
}


TEST_CASE("TestASMu2D.TransferGaussPtVarsN")
{
  SIM2D sim(1), sim2(1);
  sim.opt.discretization = sim2.opt.discretization = ASM::LRSpline;

  ASMu2D* pch = static_cast<ASMu2D*>(sim.createDefaultModel());
  LR::LRSplineSurface* lr = pch->getBasis(1);
  lr->generateIDs();

  ASMu2D* pchNew = static_cast<ASMu2D*>(sim2.createDefaultModel());
  pchNew->uniformRefine(0, 1);
  LR::LRSplineSurface* lrNew = pchNew->getBasis(1);
  REQUIRE(lrNew != nullptr);
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
  REQUIRE(refAr.size() == newAr.size());
  for (size_t i = 0; i < refAr.size(); ++i)
    REQUIRE_THAT(refAr[i], WithinRel(newAr[i]));
}


TEST_CASE("TestASMu2D.Connect")
{
  const int param = GENERATE(0,1);
  SECTION("Orient " + std::to_string(param)) {
    SIM2D sim(1);
    sim.opt.discretization = ASM::LRSpline;
    REQUIRE(sim.read(("src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient" +
                      std::to_string(param) + ".xinp").c_str()));
    REQUIRE(sim.createFEMmodel());
  }
}


TEST_CASE("TestASMu2D.ElementConnectivities")
{
  ASMuSquare pch1;
  ASMbase::resetNumbering();
  REQUIRE(pch1.uniformRefine(0,1));
  REQUIRE(pch1.uniformRefine(1,1));
  REQUIRE(pch1.generateFEMTopology());
  pch1.getBasis(1)->generateIDs();
  const size_t nel = pch1.getNoElms();
  pch1.shiftElemNumbers(nel);
  IntMat neighGlb(2*nel), neighLoc(nel);
  pch1.getElmConnectivities(neighGlb);
  pch1.getElmConnectivities(neighLoc, ASM::GEOMETRY_BASIS);
  const std::array<std::vector<int>,4> ref = {{{1, 2},
                                               {0, 3},
                                               {3, 0},
                                               {2, 1}}};
  REQUIRE(neighLoc.size() == nel);
  REQUIRE(neighGlb.size() == 2*nel);
  for (size_t n = 0; n < neighLoc.size(); ++n) {
    REQUIRE(neighLoc[n].size() == ref[n].size());
    REQUIRE(neighGlb[n+nel].size() == ref[n].size());
    for (size_t i = 0; i < neighLoc[n].size(); ++i) {
      REQUIRE(neighLoc[n][i] == ref[n][i]);
      REQUIRE(neighGlb[n+nel][i] == (ref[n][i] > -1 ?
                                      ref[n][i] + static_cast<int>(nel) : -1));
    }
  }
}


TEST_CASE("TestASMu2D.EvalPointNurbs")
{
  ASMu2D u2Dpch(2, 1);
  std::ifstream g2file("src/ASM/LR/Test/refdata/hole2D.g2");
  std::ifstream g2file2("src/ASM/LR/Test/refdata/hole2D.g2");
  ASMs2D s2Dpch(2,1);
  REQUIRE(u2Dpch.read(g2file));
  REQUIRE(s2Dpch.read(g2file2));
  REQUIRE(u2Dpch.generateFEMTopology());
  REQUIRE(s2Dpch.generateFEMTopology());

  double xi[2] = {0.1, 0.1};
  double param1[2], param2[2];
  Vec3 x1, x2;
  s2Dpch.evalPoint(xi,param1,x1);
  u2Dpch.evalPoint(xi,param2,x2);
  REQUIRE_THAT(param1[0], WithinRel(param2[0]));
  REQUIRE_THAT(param1[1], WithinRel(param2[1]));
  REQUIRE_THAT(x1[0], WithinRel(x2[0]));
  REQUIRE_THAT(x1[1], WithinRel(x2[1]));
}


TEST_CASE("TestASMu2D.Write")
{
  ASMbase::resetNumbering();
  ASMuSquare pch1;
  REQUIRE(pch1.generateFEMTopology());

  std::stringstream str;
  REQUIRE(pch1.write(str, 1));
  REQUIRE(str.str() == ASMuSquare::square);

  REQUIRE(!pch1.write(str, 2));

  str.str("");
  REQUIRE(pch1.write(str, ASM::GEOMETRY_BASIS));
  REQUIRE(str.str() == ASMuSquare::square);

  str.str("");
  REQUIRE(pch1.write(str, ASM::PROJECTION_BASIS));
  REQUIRE(str.str() == ASMuSquare::square);

  REQUIRE(!pch1.write(str, ASM::PROJECTION_BASIS_2));

  str.str("");
  REQUIRE(pch1.write(str, ASM::REFINEMENT_BASIS));
  REQUIRE(str.str() == ASMuSquare::square);

  str.str("");
  REQUIRE(pch1.write(str, ASM::INTEGRATION_BASIS));
  REQUIRE(str.str() == ASMuSquare::square);
}


TEST_CASE("TestASMu2D.ElmNodes")
{
  ASMbase::resetNumbering();

  ASMuSquare pch1;
  pch1.createProjectionBasis(true);
  pch1.raiseOrder(1,1);
  pch1.uniformRefine(0,1);
  pch1.uniformRefine(1,1);
  pch1.createProjectionBasis(false);
  REQUIRE(pch1.uniformRefine(0,1));
  REQUIRE(pch1.uniformRefine(1,1));
  REQUIRE(pch1.generateFEMTopology());

  const IntMat mnpc = pch1.getElmNodes(1);

  const auto ref = std::array{
      std::array{0,2,6,8},
      std::array{1,2,7,8},
      std::array{3,5,6,8},
      std::array{4,5,7,8},
  };
  REQUIRE(mnpc.size() == ref.size());
  for (size_t i = 0; i < mnpc.size(); ++i) {
    REQUIRE(mnpc[i].size() == ref[i].size());
    for (size_t j = 0; j < mnpc[i].size(); ++j)
      REQUIRE(mnpc[i][j] == ref[i][j]);
  }

  const auto ref_proj = std::array{
      std::array{0,2,3,8,10,11,12,14,15},
      std::array{1,2,3,9,10,11,13,14,15},
      std::array{4,6,7,8,10,11,12,14,15},
      std::array{5,6,7,9,10,11,13,14,15},
  };
  const IntMat mnpc_proj = pch1.getElmNodes(ASM::PROJECTION_BASIS);
  REQUIRE(mnpc_proj.size() == ref_proj.size());
  for (size_t i = 0; i < mnpc_proj.size(); ++i) {
    REQUIRE(mnpc_proj[i].size() == ref_proj[i].size());
    for (size_t j = 0; j < mnpc_proj[i].size(); ++j)
      REQUIRE(mnpc_proj[i][j] == ref_proj[i][j]);
  }
}


TEST_CASE("TestASMu2D.GetElementCorners")
{
  ASMuTrap<ASMu2D> pch;
  REQUIRE(pch.raiseOrder(1,1));
  REQUIRE(pch.uniformRefine(0, 1));
  REQUIRE(pch.uniformRefine(1, 1));
  REQUIRE(pch.generateFEMTopology());

  struct Ref {
    int iel, c1, c2;
    std::array<double, 2> u, v;
    std::array<Vec3,4> XC;
  };

  auto makePtArray = [](const std::array<double,2>& u,
                        const std::array<double,2>& v)
  {
    return std::array{
      u[0], v[0],
      u[1], v[0],
      u[0], v[1],
      u[1], v[1]
    };
  };

  const Ref param = GENERATE(
    Ref{
      1, 4, 1,
      {0.0, 0.5}, {0.0, 0.5},
      {Vec3{0.0, 0.0}, Vec3{3.0, 0.0}, Vec3{0.5, 1.5}, Vec3{2.75, 1.5}},
    },
    Ref{
      2, 3, 2,
      {0.5, 1.0}, {0.0, 0.5},
      {Vec3{3.0, 0.0}, Vec3{6.0, 0.0}, Vec3{2.75, 1.5}, Vec3{5.0, 1.5}},
    },
    Ref{
      3, 4, 1,
      {0.0, 0.5}, {0.5, 1.0},
      {Vec3{0.5, 1.5}, Vec3{2.75, 1.5}, Vec3{1.0, 3.0}, Vec3{2.5, 3.0}},
    },
    Ref{
      4, 3, 2,
      {0.5, 1.0}, {0.5, 1.0},
      {Vec3{2.75, 1.5}, Vec3{5.0, 1.5}, Vec3{2.5, 3.0}, Vec3{4.0, 3.0}},
    }
  );

  SECTION("Elem " + std::to_string(param.iel))
  {
    const auto& [size, XC, prm] = pch.getElementCorners(param.iel);
    REQUIRE(XC.size() == 4);
    REQUIRE(prm.size() == 8);
    using Catch::Matchers::RangeEquals;
    REQUIRE_THAT(size, WithinRel((XC[param.c1-1] - XC[param.c2-1]).length()));
    REQUIRE_THAT(prm, RangeEquals(makePtArray(param.u, param.v)));
    REQUIRE_THAT(XC, RangeEquals(param.XC));
  }
}
