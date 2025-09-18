//==============================================================================
//!
//! \file TestASMu3D.C
//!
//! \date Jul 14 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for driver for assembly of unstructured 3D spline FE models.
//!
//==============================================================================

#include "ASMuCube.h"
#include "SIM3D.h"
#include "GaussQuadrature.h"
#include "LRSpline/LRSplineVolume.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <numeric>
#include <sstream>

using Catch::Matchers::WithinRel;


TEST_CASE("TestASMu3D.BoundaryNodes")
{
  const int param = GENERATE(1,2,3,4,5,6);

  SECTION("Face " + std::to_string(param)) {
    SIM3D sim(1);
    sim.opt.discretization = ASM::LRSpline;
    REQUIRE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes_3d.xinp"));
    REQUIRE(sim.createFEMmodel());

    std::stringstream str;
    str << "Face" << param;
    int bcode = sim.getUniquePropertyCode(str.str(),0);
    std::vector<int> vec;
    sim.getBoundaryNodes(bcode,vec);
    REQUIRE(vec.size() == 16);
    int i, j, k = 0;
    for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j, ++k)
        switch (param) {
        case 1: REQUIRE(vec[k] == 1+4*(4*i+j)); break;
        case 2: REQUIRE(vec[k] == 2+4*(4*i+j)); break;
        case 3: REQUIRE(vec[k] == 16*i+j+1); break;
        case 4: REQUIRE(vec[k] == 5+16*i+j); break;
        case 5: REQUIRE(vec[k] == 4*i+j+1); break;
        case 6: REQUIRE(vec[k] == 17+4*i+j); break;
        }
    }
}


TEST_CASE("TestASMu3D.Connect")
{
  const int param = GENERATE(0,1,2,3,4,5,6,7);

  SECTION("Orient " + std::to_string(param)) {
    SIM3D sim(3);
    sim.opt.discretization = ASM::LRSpline;
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
    str << param << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.createFEMmodel());
  }
}


TEST_CASE("TestASMu3D.ConnectUneven")
{
  const int param = GENERATE(0,1,2,3,4,5,6,7,8,9,10,11,12);

  SECTION("Uneven " + std::to_string(param)) {
    SIM3D sim(1);
    sim.opt.discretization = ASM::LRSpline;
    std::stringstream str;
    str << "src/ASM/Test/refdata/3d_uneven";
    int idx = param-1;
    if (idx >= 0)
      str << "_" << idx/3 << idx%3;
    str << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.createFEMmodel());
  }
}


TEST_CASE("TestASMu3D.TransferGaussPtVars")
{
  ASMuCube pch;
  LR::LRSplineVolume* lr = pch.getBasis(1);
  REQUIRE(lr != nullptr);
  lr->generateIDs();

  RealArray oldAr(3*3*3), newAr;
  const double* xi = GaussQuadrature::getCoord(3);
  size_t id[3];
  for (size_t idx = 0; idx < 3; ++idx) {
    ASMuCube pchNew;
    pchNew.uniformRefine(idx,1);
    pchNew.getBasis(1)->generateIDs();
    for (id[2] = 0; id[2] < 3; ++id[2])
      for (id[1] = 0; id[1] < 3; ++id[1])
        for (id[0] = 0; id[0] < 3; ++id[0])
          oldAr[id[0]+(id[1]+id[2]*3)*3] = (1.0 + xi[id[idx]]) / 2.0;
    pchNew.transferGaussPtVars(lr, oldAr, newAr, 3);
    size_t k = 0;
    for (size_t iEl = 0; iEl < 2; ++iEl)
      for (id[2] = 0; id[2] < 3; ++id[2])
        for (id[1] = 0; id[1] < 3; ++id[1])
          for (id[0] = 0; id[0] < 3; ++id[0], ++k)
            REQUIRE_THAT(newAr[k], WithinRel(0.5*iEl + 0.5*(xi[id[idx]] + 1.0) / 2.0));
  }
}


TEST_CASE("TestASMu3D.TransferGaussPtVarsN")
{
  ASMuCube pch, pchNew;
  LR::LRSplineVolume* lr = pch.getBasis(1);
  REQUIRE(lr != nullptr);
  lr->generateIDs();

  pchNew.uniformRefine(0,1);
  pchNew.getBasis(1)->generateIDs();

  RealArray oldAr(3*3*3), newAr;
  std::iota(oldAr.begin(), oldAr.end(), 1);

  pchNew.transferGaussPtVarsN(lr, oldAr, newAr, 3);
  static RealArray refAr = {{ 1.0,  1.0,  2.0,
                              4.0,  4.0,  5.0,
                              7.0,  7.0,  8.0,
                             10.0, 10.0, 11.0,
                             13.0, 13.0, 14.0,
                             16.0, 16.0, 17.0,
                             19.0, 19.0, 20.0,
                             22.0, 22.0, 23.0,
                             25.0, 25.0, 26.0,
                              2.0,  3.0,  3.0,
                              5.0,  6.0,  6.0,
                              8.0,  9.0,  9.0,
                             11.0, 12.0, 12.0,
                             14.0, 15.0, 15.0,
                             17.0, 18.0, 18.0,
                             20.0, 21.0, 21.0,
                             23.0, 24.0, 24.0,
                             26.0, 27.0, 27.0}};
  REQUIRE(refAr.size() == newAr.size());
  for (size_t i = 0; i < refAr.size(); ++i)
    REQUIRE_THAT(refAr[i], WithinRel(newAr[i]));
}


TEST_CASE("TestASMu3D.ElementConnectivities")
{
  ASMuCube pch1;
  ASMbase::resetNumbering();
  REQUIRE(pch1.uniformRefine(0,1));
  REQUIRE(pch1.uniformRefine(1,1));
  REQUIRE(pch1.uniformRefine(2,1));
  REQUIRE(pch1.generateFEMTopology());
  const size_t nel = pch1.getNoElms();
  IntMat neighGlb(2*nel), neighLoc(nel);
  pch1.shiftElemNumbers(nel);
  pch1.getElmConnectivities(neighGlb);
  pch1.getElmConnectivities(neighLoc, ASM::GEOMETRY_BASIS);
  const std::array<std::vector<int>,8> ref = {{{1, 2, 4},
                                               {0, 3, 5},
                                               {3, 0, 6},
                                               {2, 1, 7},
                                               {5, 6, 0},
                                               {4, 7, 1},
                                               {7, 4, 2},
                                               {6, 5, 3}}};
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


class ASMu3DTest : public ASMuCube
{
public:
  void getFaceCorners(std::array<int,4>& corners, int dir)
  {
    DirichletFace df(this->getBasis(),dir);
    for (int i = 0; i < 4; i++)
      corners[i] = df.corners[i];
  }
};


TEST_CASE("TestASMu3D.DirichletFace")
{
  ASMu3DTest pch1;
  ASMbase::resetNumbering();
  REQUIRE(pch1.generateFEMTopology());

  const std::array<std::array<int,4>,7> refArr = {{
      {{ 1, 2, 3, 4 }}, // Bottom
      {{ 1, 2, 5, 6 }}, // South
      {{ 1, 3, 5, 7 }}, // West
      {{ 0, 0, 0, 0 }}, // None
      {{ 2, 4, 6, 8 }}, // East
      {{ 3, 4, 7, 8 }}, // North
      {{ 5, 6, 7, 8 }}  // Top
    }};

  std::array<int,4> corners;
  for (int dir = -3; dir <= 3; dir++)
  {
    pch1.getFaceCorners(corners,dir);
    REQUIRE(corners == refArr[dir+3]);
  }
}


TEST_CASE("TestASMu3D.Write")
{
  ASMbase::resetNumbering();
  ASMuCube pch1;
  REQUIRE(pch1.generateFEMTopology());

  std::stringstream str;
  REQUIRE(pch1.write(str, 1));
  REQUIRE(str.str() == ASMuCube::cube);

  REQUIRE(!pch1.write(str, 2));

  str.str("");
  REQUIRE(pch1.write(str, ASM::GEOMETRY_BASIS));
  REQUIRE(str.str() == ASMuCube::cube);

  str.str("");
  REQUIRE(pch1.write(str, ASM::PROJECTION_BASIS));
  REQUIRE(str.str() == ASMuCube::cube);

  REQUIRE(!pch1.write(str, ASM::PROJECTION_BASIS_2));

  str.str("");
  REQUIRE(pch1.write(str, ASM::REFINEMENT_BASIS));
  REQUIRE(str.str() == ASMuCube::cube);

  str.str("");
  REQUIRE(pch1.write(str, ASM::INTEGRATION_BASIS));
  REQUIRE(str.str() == ASMuCube::cube);
}


TEST_CASE("TestASMu3D.ElmNodes")
{
  ASMbase::resetNumbering();

  ASMuCube pch1;
  pch1.createProjectionBasis(true);
  pch1.raiseOrder(1,1,1,false);
  pch1.uniformRefine(0,1);
  pch1.uniformRefine(1,1);
  pch1.uniformRefine(2,1);
  pch1.createProjectionBasis(false);
  REQUIRE(pch1.uniformRefine(0,1));
  REQUIRE(pch1.uniformRefine(1,1));
  REQUIRE(pch1.uniformRefine(2,1));
  REQUIRE(pch1.generateFEMTopology());

  const IntMat mnpc = pch1.getElmNodes(1);

  const auto ref = std::array{
      std::array{0,2,6,8,18,20,24,26},
      std::array{1,2,7,8,19,20,25,26},
      std::array{3,5,6,8,21,23,24,26},
      std::array{4,5,7,8,22,23,25,26},
      std::array{9,11,15,17,18,20,24,26},
      std::array{10,11,16,17,19,20,25,26},
      std::array{12,14,15,17,21,23,24,26},
      std::array{13,14,16,17,22,23,25,26},
  };
  REQUIRE(mnpc.size() == ref.size());
  for (size_t i = 0; i < mnpc.size(); ++i) {
    REQUIRE(mnpc[i].size() == ref[i].size());
    for (size_t j = 0; j < mnpc[i].size(); ++j)
      REQUIRE(mnpc[i][j] == ref[i][j]);
  }

  const auto ref_proj = std::array{
      std::array{0,1,2,4,5,6,8,9,10,16,17,18,20,21,22,24,25,26,32,33,34,36,37,38,40,41,42},
      std::array{1,2,3,5,6,7,9,10,11,17,18,19,21,22,23,25,26,27,33,34,35,37,38,39,41,42,43},
      std::array{4,5,6,8,9,10,12,13,14,20,21,22,24,25,26,28,29,30,36,37,38,40,41,42,44,45,46},
      std::array{5,6,7,9,10,11,13,14,15,21,22,23,25,26,27,29,30,31,37,38,39,41,42,43,45,46,47},
      std::array{16,17,18,20,21,22,24,25,26,32,33,34,36,37,38,40,41,42,48,49,50,52,53,54,56,57,58},
      std::array{17,18,19,21,22,23,25,26,27,33,34,35,37,38,39,41,42,43,49,50,51,53,54,55,57,58,59},
      std::array{20,21,22,24,25,26,28,29,30,36,37,38,40,41,42,44,45,46,52,53,54,56,57,58,60,61,62},
      std::array{21,22,23,25,26,27,29,30,31,37,38,39,41,42,43,45,46,47,53,54,55,57,58,59,61,62,63},
  };
  const IntMat mnpc_proj = pch1.getElmNodes(ASM::PROJECTION_BASIS);
  REQUIRE(mnpc_proj.size() == ref_proj.size());
  for (size_t i = 0; i < mnpc_proj.size(); ++i) {
    REQUIRE(mnpc_proj[i].size() == ref_proj[i].size());
    for (size_t j = 0; j < mnpc_proj[i].size(); ++j)
      REQUIRE(mnpc_proj[i][j] == ref_proj[i][j]);
  }
}
