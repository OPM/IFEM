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
#include "ASMuTrap3.h"
#include "GaussQuadrature.h"
#include "Test/ASM3DTests.h"
#include "Vec3Oper.h"

#include "Catch2Support.h"

#include <LRSpline/LRSplineVolume.h>

#include <functional>
#include <numeric>


TEST_CASE("TestASMu3D.BoundaryElements")
{
  ASM3DTests<ASMuCube>::BoundaryElements();
}


TEST_CASE("TestASMu3D.BoundaryNodes")
{
  using Func = std::function<int(int,int)>;
  const auto ref = std::array{
    Func{[](int i, int j) { return 1+4*(4*j+i); }},
    Func{[](int i, int j) { return 2+4*(4*j+i); }},
    Func{[](int i, int j) { return 16*j+i+1; }},
    Func{[](int i, int j) { return 5+16*j+i; }},
    Func{[](int i, int j) { return 4*j+i+1; }},
    Func{[](int i, int j) { return 17+4*j+i; }},
  };

  ASM3DTests<ASMuCube>::BoundaryNodes(ref);
}


TEST_CASE("TestASMu3D.Connect")
{
  ASM3DTests<ASMuCube>::Connect(ASM::LRSpline);
}


TEST_CASE("TestASMu3D.ConnectUneven")
{
  ASM3DTests<ASMuCube>::ConnectUneven(ASM::LRSpline);
}


TEST_CASE("TestASMu3D.ConstrainFace")
{
  ASM3DTests<ASMuCube>::ConstrainFace();
}


TEST_CASE("TestASMu3D.ElementConnectivities")
{
  const auto ref = std::array{
    std::array{1, 2, 4},
    std::array{0, 3, 5},
    std::array{3, 0, 6},
    std::array{2, 1, 7},
    std::array{5, 6, 0},
    std::array{4, 7, 1},
    std::array{7, 4, 2},
    std::array{6, 5, 3},
  };

  ASM3DTests<ASMuCube>::ElementConnectivities(ref);
}


TEST_CASE("TestASMu3D.ElmNodes")
{
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

  ASM3DTests<ASMuCube>::ElmNodes(ref, ref_proj);
}


TEST_CASE("TestASMu3D.FaceCorners")
{
  ASM3DTests<ASMuCube>::FaceCorners();
}


TEST_CASE("TestASMu3D.FaceIntegrate")
{
  ASM3DTests<ASMuCube>::FaceIntegrate();
}


TEST_CASE("TestASMu3D.GetElementCorners")
{
  ASM3DTests<ASMuTrap3<ASMu3D>>::GetElementCorners();
}


TEST_CASE("TestASMu3D.TransferCtrlPtVars")
{
  // x*y*z + 0.5
  const auto orig_scalar = std::vector{
    0.0 * 0.0 * 0.0 + 0.5,
    1.0 * 0.0 * 0.0 + 0.5,
    0.0 * 1.0 * 0.0 + 0.5,
    1.0 * 1.0 * 0.0 + 0.5,
    0.0 * 0.0 * 1.0 + 0.5,
    1.0 * 0.0 * 1.0 + 0.5,
    0.0 * 1.0 * 1.0 + 0.5,
    1.0 * 1.0 * 1.0 + 0.5,
  };

  // x*y*z + 0.5, x*y*z - 0.5
  const auto orig_vec = std::vector{
    0.0 * 0.0 * 0.0 + 0.5,
    0.0 * 0.0 * 0.0 - 0.5,
    1.0 * 0.0 * 0.0 + 0.5,
    1.0 * 0.0 * 0.0 - 0.5,
    0.0 * 1.0 * 0.0 + 0.5,
    0.0 * 1.0 * 0.0 - 0.5,
    1.0 * 1.0 * 0.0 + 0.5,
    1.0 * 1.0 * 0.0 - 0.5,
    0.0 * 0.0 * 1.0 + 0.5,
    0.0 * 0.0 * 1.0 - 0.5,
    1.0 * 0.0 * 1.0 + 0.5,
    1.0 * 0.0 * 1.0 - 0.5,
    0.0 * 1.0 * 1.0 + 0.5,
    0.0 * 1.0 * 1.0 - 0.5,
    1.0 * 1.0 * 1.0 + 0.5,
    1.0 * 1.0 * 1.0 - 0.5,
  };

  ASMuCube pch;
  REQUIRE(pch.generateFEMTopology());

  const Real*  xi = GaussQuadrature::getCoord(2);
  for (size_t idx = 0; idx < 3; ++idx) {
    ASMuCube pchNew;
    REQUIRE(pchNew.uniformRefine(idx, 2));
    REQUIRE(pchNew.generateFEMTopology());
    RealArray new_scalar, new_vec;
    REQUIRE(pchNew.transferCntrlPtVars(pch.getBasis(), orig_scalar, new_scalar, 2, 1));
    REQUIRE(pchNew.transferCntrlPtVars(pch.getBasis(), orig_vec, new_vec, 2, 2));
    REQUIRE(new_scalar.size() == pchNew.getNoElms() * 8);
    REQUIRE(new_vec.size() == pchNew.getNoElms() * 8 * 2);

    size_t pt = 0;
    for (size_t iel = 1; iel <= pchNew.getNoElms(); ++iel) {
      std::array<RealArray, 3> uGP;
      LR::getGaussPointParameters(pchNew.getBasis(), uGP[0], 0, 2, iel, xi);
      LR::getGaussPointParameters(pchNew.getBasis(), uGP[1], 1, 2, iel, xi);
      LR::getGaussPointParameters(pchNew.getBasis(), uGP[2], 2, 2, iel, xi);
      for (size_t k = 0; k < 2; ++k)
        for (size_t j = 0; j < 2; ++j)
          for (size_t i = 0; i < 2; ++i, ++pt) {
            double prm[3] = {uGP[0][i], uGP[1][j], uGP[2][k]};
            Vec3 pos;
            pchNew.evalPoint(prm, nullptr, pos);
            REQUIRE_THAT(new_scalar[pt],  WithinRel(pos.x * pos.y * pos.z + 0.5));
            REQUIRE_THAT(new_vec[2*pt],   WithinRel(pos.x * pos.y * pos.z + 0.5));
            REQUIRE_THAT(new_vec[2*pt+1], WithinRel(pos.x * pos.y * pos.z - 0.5));
        }
    }
  }
}


TEST_CASE("TestASMu3D.TransferGaussPtVars")
{
  ASMuCube pch;
  REQUIRE(pch.generateFEMTopology());

  RealArray oldAr(3*3*3), newAr;
  const double* xi = GaussQuadrature::getCoord(3);
  size_t id[3];
  for (size_t idx = 0; idx < 3; ++idx) {
    ASMuCube pchNew;
    REQUIRE(pchNew.uniformRefine(idx,1));
    REQUIRE(pchNew.generateFEMTopology());
    for (id[2] = 0; id[2] < 3; ++id[2])
      for (id[1] = 0; id[1] < 3; ++id[1])
        for (id[0] = 0; id[0] < 3; ++id[0])
          oldAr[id[0]+(id[1]+id[2]*3)*3] = (1.0 + xi[id[idx]]) / 2.0;
    pchNew.transferGaussPtVars(pch.getBasis(1), oldAr, newAr, 3);
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
  const auto refAr = std::array{
    1.0,  1.0,  2.0,
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
    26.0, 27.0, 27.0,
  };

  ASMuCube pch, pchNew;
  REQUIRE(pch.generateFEMTopology());
  REQUIRE(pchNew.uniformRefine(0,1));
  REQUIRE(pchNew.generateFEMTopology());

  RealArray oldAr(3*3*3), newAr;
  std::iota(oldAr.begin(), oldAr.end(), 1);

  REQUIRE(pchNew.transferGaussPtVarsN(pch.getBasis(1), oldAr, newAr, 3));
  REQUIRE(refAr.size() == newAr.size());
  for (size_t i = 0; i < refAr.size(); ++i)
    REQUIRE_THAT(refAr[i], WithinRel(newAr[i]));
}


TEST_CASE("TestASMu3D.Write")
{
  ASM3DTests<ASMuCube>::Write(true);
}
