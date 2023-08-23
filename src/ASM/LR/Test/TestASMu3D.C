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

#include "ASMu3D.h"
#include "SIM3D.h"
#include "GaussQuadrature.h"
#include "LRSpline/LRSplineVolume.h"

#include "gtest/gtest.h"
#include <numeric>

class TestASMu3D : public testing::Test,
                   public testing::WithParamInterface<int>
{
};


TEST_P(TestASMu3D, BoundaryNodes)
{
  if (GetParam() < 1 || GetParam() > 6)
    return;

  SIM3D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes_3d.xinp"));
  ASSERT_TRUE(sim.createFEMmodel());

  std::stringstream str;
  str << "Face" << GetParam();
  int bcode = sim.getUniquePropertyCode(str.str(),0);
  std::vector<int> vec;
  sim.getBoundaryNodes(bcode,vec);
  ASSERT_EQ(vec.size(), 16U);
  int i, j, k = 0;
  for (i = 0; i < 4; ++i)
    for (j = 0; j < 4; ++j, ++k)
      switch (GetParam()) {
      case 1: EXPECT_EQ(vec[k], 1+4*(4*i+j)); break;
      case 2: EXPECT_EQ(vec[k], 2+4*(4*i+j)); break;
      case 3: EXPECT_EQ(vec[k], 16*i+j+1); break;
      case 4: EXPECT_EQ(vec[k], 5+16*i+j); break;
      case 5: EXPECT_EQ(vec[k], 4*i+j+1); break;
      case 6: EXPECT_EQ(vec[k], 17+4*i+j); break;
      }
}


TEST_P(TestASMu3D, Connect)
{
  if (GetParam() > 7)
    return;

  SIM3D sim(3);
  sim.opt.discretization = ASM::LRSpline;
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.createFEMmodel());
}


TEST_P(TestASMu3D, ConnectUneven)
{
  SIM3D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  std::stringstream str;
  str << "src/ASM/Test/refdata/3d_uneven";
  int idx = GetParam()-1;
  if (idx >= 0)
    str << "_" << idx/3 << idx%3;
  str << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.createFEMmodel());
}


const std::vector<int> tests = {0,1,2,3,4,5,6,7,8,9,10,11,12};
INSTANTIATE_TEST_CASE_P(TestASMu3D,
                        TestASMu3D,
                        testing::ValuesIn(tests));


class ASMuCube : public ASMu3D
{
public:
  ASMuCube()
  {
    std::stringstream geo("700 1 0 0\n3 0\n"
                          "2 2\n0 0 1 1\n"
                          "2 2\n0 0 1 1\n"
                          "2 2\n0 0 1 1\n"
                          "0 0 0\n1 0 0\n0 1 0\n1 1 0\n"
                          "0 0 1\n1 0 1\n0 1 1\n1 1 1\n");
    EXPECT_TRUE(this->read(geo));
  }
  virtual ~ASMuCube() {}
};


TEST(TestASMu3D, TransferGaussPtVars)
{
  ASMuCube pch;
  LR::LRSplineVolume* lr = pch.getBasis(1);
  ASSERT_TRUE(lr != nullptr);
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
            EXPECT_FLOAT_EQ(newAr[k], 0.5*iEl + 0.5*(xi[id[idx]] + 1.0) / 2.0);
  }
}


TEST(TestASMu3D, TransferGaussPtVarsN)
{
  ASMuCube pch, pchNew;
  LR::LRSplineVolume* lr = pch.getBasis(1);
  ASSERT_TRUE(lr != nullptr);
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
  ASSERT_EQ(refAr.size(), newAr.size());
  for (size_t i = 0; i < refAr.size(); ++i)
    EXPECT_FLOAT_EQ(refAr[i], newAr[i]);
}


TEST(TestASMu3D, ElementConnectivities)
{
  ASMuCube pch1;
  ASMbase::resetNumbering();
  ASSERT_TRUE(pch1.uniformRefine(0,1));
  ASSERT_TRUE(pch1.uniformRefine(1,1));
  ASSERT_TRUE(pch1.uniformRefine(2,1));
  ASSERT_TRUE(pch1.generateFEMTopology());
  IntMat neigh(8);
  pch1.getElmConnectivities(neigh);
  const std::array<std::vector<int>,8> ref = {{{1, 2, 4},
                                               {0, 3, 5},
                                               {3, 0, 6},
                                               {2, 1, 7},
                                               {5, 6, 0},
                                               {4, 7, 1},
                                               {7, 4, 2},
                                               {6, 5, 3}}};
  ASSERT_EQ(neigh.size(), 8U);
  for (size_t n = 0; n < neigh.size(); ++n) {
    ASSERT_EQ(neigh[n].size(), ref[n].size());
    for (size_t i = 0; i < neigh[n].size(); ++i)
      EXPECT_EQ(neigh[n][i], ref[n][i]);
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


TEST(TestASMu3D, DirichletFace)
{
  ASMu3DTest pch1;
  ASMbase::resetNumbering();
  ASSERT_TRUE(pch1.generateFEMTopology());

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
    EXPECT_EQ(corners,refArr[dir+3]);
  }
}
