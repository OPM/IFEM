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
#include "GlobalIntegral.h"
#include "LaplacianIntegrand.h"
#include "TimeDomain.h"
#include "Test/ASM2DTests.h"

#include "Catch2Support.h"

#include <LRSpline/LRSplineSurface.h>

#include <fstream>
#include <numeric>
#include <sstream>


TEST_CASE("TestASMu2D.BoundaryElements")
{
  ASM2DTests<ASMuSquare>::BoundaryElements();
}


TEST_CASE("TestASMu2D.Connect")
{
  ASM2DTests<ASMuSquare>::Connect(ASM::LRSpline);
}


TEST_CASE("TestASMu2D.ConstrainEdge")
{
  ASM2DTests<ASMuSquare>::ConstrainEdge();
}


TEST_CASE("TestASMu2D.ConstrainEdgeOpen")
{
  ASM2DTests<ASMuSquare>::ConstrainEdgeOpen();
}


TEST_CASE("TestASMu2D.ElementConnectivities")
{
  const auto ref = std::array{
    std::array{1, 2},
    std::array{0, 3},
    std::array{3, 0},
    std::array{2, 1},
  };

  ASM2DTests<ASMuSquare>::GetElementConnectivities(ref);
}


TEST_CASE("TestASMu2D.ElmNodes")
{
  const auto ref = std::array{
    std::array{0,2,6,8},
    std::array{1,2,7,8},
    std::array{3,5,6,8},
    std::array{4,5,7,8},
  };
  const auto ref_proj = std::array{
    std::array{0,2,3,8,10,11,12,14,15},
    std::array{1,2,3,9,10,11,13,14,15},
    std::array{4,6,7,8,10,11,12,14,15},
    std::array{5,6,7,9,10,11,13,14,15},
  };

  ASM2DTests<ASMuSquare>::ElmNodes(ref, ref_proj);
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


TEST_CASE("TestASMu2D.SecondDerivatives")
{
  // A curved biquadratic surface, as a tensor-product and as an LR B-spline
  // patch, which then have the same second derivatives
  const char* g2 = "200 1 0 0\n2 0\n"
                   "4 3\n0 0 0 0.4 1 1 1\n"
                   "4 3\n0 0 0 0.7 1 1 1\n"
                   "0 0\n0.3 0.1\n0.7 -0.1\n1 0\n"
                   "0 0.3\n0.4 0.4\n0.6 0.2\n1.1 0.3\n"
                   "-0.1 0.7\n0.3 0.6\n0.8 0.8\n1 0.7\n"
                   "0 1\n0.4 1.1\n0.7 0.9\n1 1\n";
  ASMs2D spch(2,1);
  ASMu2D upch(2,1);
  for (ASMbase* pch : {static_cast<ASMbase*>(&spch), static_cast<ASMbase*>(&upch)})
  {
    std::istringstream is(g2);
    REQUIRE(pch->read(is));
    ASMbase::resetNumbering();
    REQUIRE(pch->generateFEMTopology());
  }

  // The integral of the square of the Laplacian of a field
  LaplacianIntegrand sitg(spch), uitg(upch);
  GlobalIntegral dummy;
  static_cast<ASMbase&>(spch).generateThreadGroups(sitg,true,false);
  static_cast<ASMbase&>(upch).generateThreadGroups(uitg,true,false);
  REQUIRE(spch.integrate(sitg,dummy,TimeDomain()));
  REQUIRE(upch.integrate(uitg,dummy,TimeDomain()));
  CHECK_THAT(uitg.result, WithinRel(sitg.result, 1e-12));

  // The second derivatives of the field at some points
  RealArray gpar[2] = {{0.1, 0.35, 0.6, 0.9}, {0.2, 0.5, 0.75, 0.95}};
  Vector ssol(spch.getNoNodes()), usol(upch.getNoNodes());
  for (size_t i = 1; i <= ssol.size(); i++)
  {
    const Vec3 X = spch.getCoord(i);
    ssol(i) = std::sin(X.x) + X.x*X.y*X.y;
  }
  for (size_t i = 1; i <= usol.size(); i++)
  {
    const Vec3 X = upch.getCoord(i);
    usol(i) = std::sin(X.x) + X.x*X.y*X.y;
  }
  Matrix sder, uder;
  REQUIRE(spch.evalSolution(sder,ssol,gpar,false,2));
  REQUIRE(upch.evalSolution(uder,usol,gpar,false,2));
  REQUIRE(uder.rows() == sder.rows());
  REQUIRE(uder.cols() == sder.cols());
  for (size_t i = 1; i <= sder.rows(); i++)
    for (size_t j = 1; j <= sder.cols(); j++)
      CHECK_THAT(uder(i,j), WithinAbs(sder(i,j), 1e-10));
}


TEST_CASE("TestASMu2D.GetElementCorners")
{
  ASM2DTests<ASMuTrap<ASMu2D>>::GetElementCorners();
}


TEST_CASE("TestASMu2D.InterfaceChecker")
{
  const auto elem_flags = std::array{
    10, 11, 11,  9,
    14, 15, 15, 13,
    14, 15, 15, 13,
     6,  7, 15, 13,
    15, 13,
     7,  5
  };
  const auto line_cont = std::array{
    std::array{-1,  1, -1,  1},
    std::array{ 1,  0, -1,  1},
    std::array{ 0,  1, -1,  0},
    std::array{ 1, -1, -1,  0},
    std::array{-1,  1,  1,  1},
    std::array{ 1,  0,  1,  1},
    std::array{ 0,  1,  0,  1},
    std::array{ 1, -1,  0,  1},
    std::array{-1,  1,  1,  1},
    std::array{ 1,  0,  1,  1},
    std::array{ 0,  1,  1,  1},
    std::array{ 1, -1,  1,  1},
    std::array{-1,  1,  1, -1},
    std::array{ 1,  0,  1, -1},
    std::array{ 0,  1,  1,  1},
    std::array{ 1, -1,  1,  1},
    std::array{ 0,  1,  0,  0},
    std::array{ 1, -1,  0,  0},
    std::array{ 0,  1,  1, -1},
    std::array{ 1, -1,  1, -1},
  };
  using RA = RealArray;
  const auto elem_pts = std::array{
    std::array{RA{},      RA{0.25},        RA{},     RA{0.25}},
    std::array{RA{0.25},  RA{0.125, 0.25}, RA{},     RA{0.5}},
    std::array{RA{0.125}, RA{0.125},       RA{},     RA{0.75}},
    std::array{RA{0.125}, RA{},            RA{},     RA{1.0}},
    std::array{RA{},      RA{0.5},         RA{0.25}, RA{0.25}},
    std::array{RA{0.5},   RA{0.5},         RA{0.5},  RA{0.5}},
    std::array{RA{0.5},   RA{0.5},         RA{0.75}, RA{0.75}},
    std::array{RA{0.5},   RA{0.5},         RA{1.0},  RA{1.0}},
    std::array{RA{},      RA{0.75},        RA{0.25}, RA{0.25}},
    std::array{RA{0.75},  RA{0.75},        RA{0.5},  RA{0.5}},
    std::array{RA{0.75},  RA{0.75},        RA{0.75}, RA{0.75}},
    std::array{RA{0.75},  RA{},            RA{1.0},  RA{1.0}},
    std::array{RA{},      RA{1.0},         RA{0.25}, RA{}},
    std::array{RA{1.0},   RA{0.875, 1.0},  RA{0.5},  RA{0.5}},
    std::array{RA{0.875}, RA{0.875},       RA{0.75}, RA{0.75}},
    std::array{RA{0.875}, RA{},            RA{1.0},  RA{1.0}},
    std::array{RA{0.25},  RA{0.25},        RA{0.75}, RA{0.75}},
    std::array{RA{0.25},  RA{},            RA{1.0},  RA{1.0}},
    std::array{RA{1.0},   RA{1.0},         RA{0.75}, RA{}},
    std::array{RA{1.0},   RA{},            RA{1.0},  RA{}},
  };

  ASMuSquare pch;
  REQUIRE(pch.raiseOrder(1,1));
  REQUIRE(pch.uniformRefine(0, 3));
  REQUIRE(pch.uniformRefine(1, 3));
  LR::LRSplineSurface* lr = pch.getBasis();
  lr->insert_const_u_edge(0.5,   0, 1, 2);
  lr->insert_const_v_edge(0.125, 0.5, 1, 2);
  lr->insert_const_v_edge(0.25,  0.5, 1, 2);
  lr->insert_const_v_edge(0.875, 0.5, 1, 1);
  REQUIRE(pch.generateFEMTopology());

  REQUIRE(elem_flags.size() == pch.getNoElms());
  REQUIRE(elem_pts.size() == pch.getNoElms());
  ASMu2D::InterfaceChecker iChk(pch);
  for (size_t i = 0; i < pch.getNoElms(); ++i) {
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


TEST_CASE("TestASMu2D.TransferCtrlPtVars")
{
  // x*y + 0.5
  const RealArray orig_scalar{
    0.0 * 0.0 + 0.5,
    1.0 * 0.0 + 0.5,
    0.0 * 1.0 + 0.5,
    1.0 * 1.0 + 0.5,
  };

  // x*y + 0.5, x*y - 0.5
  const RealArray orig_vec{
    0.0 * 0.0 + 0.5,
    0.0 * 0.0 - 0.5,
    1.0 * 0.0 + 0.5,
    1.0 * 0.0 - 0.5,
    0.0 * 1.0 + 0.5,
    0.0 * 1.0 - 0.5,
    1.0 * 1.0 + 0.5,
    1.0 * 1.0 - 0.5,
  };

  ASMuSquare pch;
  REQUIRE(pch.generateFEMTopology());

  const Real*  xi = GaussQuadrature::getCoord(2);
  for (size_t idx = 0; idx < 2; ++idx) {
    ASMuSquare pchNew;
    REQUIRE(pchNew.uniformRefine(idx, 2));
    REQUIRE(pchNew.generateFEMTopology());
    RealArray new_scalar, new_vec;
    pchNew.transferCntrlPtVars(pch.getBasis(), orig_scalar, new_scalar, 2, 1);
    pchNew.transferCntrlPtVars(pch.getBasis(), orig_vec, new_vec, 2, 2);
    REQUIRE(new_scalar.size() == pchNew.getNoElms() * 4);
    REQUIRE(new_vec.size() == pchNew.getNoElms() * 4 * 2);

    size_t pt = 0;
    for (size_t iel = 1; iel <= pchNew.getNoElms(); ++iel) {
      std::array<RealArray, 2> uGP;
      LR::getGaussPointParameters(pchNew.getBasis(), uGP[0], 0, 2, iel, xi);
      LR::getGaussPointParameters(pchNew.getBasis(), uGP[1], 1, 2, iel, xi);
      for (size_t j = 0; j < 2; ++j)
        for (size_t i = 0; i < 2; ++i, ++pt) {
          double prm[2] = {uGP[0][i], uGP[1][j]};
          Vec3 pos;
          pchNew.evalPoint(prm, nullptr, pos);
          REQUIRE_THAT(new_scalar[pt], WithinRel(pos.x * pos.y + 0.5));
          REQUIRE_THAT(new_vec[2*pt], WithinRel(pos.x * pos.y + 0.5));
          REQUIRE_THAT(new_vec[2*pt+1], WithinRel(pos.x * pos.y - 0.5));
        }
    }
  }
}


TEST_CASE("TestASMu2D.TransferGaussPtVars")
{
  ASMuSquare pch;
  REQUIRE(pch.generateFEMTopology());

  RealArray oldAr(9), newAr;
  const double* xi = GaussQuadrature::getCoord(3);
  size_t id[2];
  for (size_t idx = 0; idx < 2; ++idx) {
    ASMuSquare pchNew;
    REQUIRE(pchNew.uniformRefine(idx, 1));
    REQUIRE(pchNew.generateFEMTopology());
    for (id[1] = 0; id[1] < 3; ++id[1])
      for (id[0] = 0; id[0] < 3; ++id[0])
        oldAr[id[0]+id[1]*3] = (1.0 + xi[id[idx]]) / 2.0;
    pchNew.transferGaussPtVars(pch.getBasis(1), oldAr, newAr, 3);
    size_t k = 0;
    for (size_t iEl = 0; iEl < 2; ++iEl)
      for (id[1] = 0; id[1] < 3; ++id[1])
        for (id[0] = 0; id[0] < 3; ++id[0], ++k)
          REQUIRE_THAT(newAr[k], WithinRel(0.5*iEl + 0.5*(xi[id[idx]] + 1.0) / 2.0));
  }
}


TEST_CASE("TestASMu2D.TransferGaussPtVarsN")
{
  ASMuSquare pch, pchNew;
  REQUIRE(pch.generateFEMTopology());
  REQUIRE(pchNew.uniformRefine(0,1));
  REQUIRE(pchNew.generateFEMTopology());

  RealArray oldAr(9), newAr;
  std::iota(oldAr.begin(), oldAr.end(), 1);

  pchNew.transferGaussPtVarsN(pch.getBasis(1), oldAr, newAr, 3);
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


TEST_CASE("TestASMu2D.Write")
{
  ASM2DTests<ASMuSquare>::Write(true);
}
