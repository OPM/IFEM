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
  if (GetParam() == 0 || GetParam() > 6)
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
  auto it = vec.begin();
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      if (GetParam() == 1)
        EXPECT_EQ(*it++, 1+4*(4*i+j));
      else if (GetParam() == 2)
        EXPECT_EQ(*it++, 2+4*(4*i+j));
      else if (GetParam() == 3)
        EXPECT_EQ(*it++, 16*i+j+1);
      else if (GetParam() == 4)
        EXPECT_EQ(*it++, 5+16*i+j);
      else if (GetParam() == 5)
        EXPECT_EQ(*it++, 4*i+j+1);
      else if (GetParam() == 6)
        EXPECT_EQ(*it++, 17+4*i+j);
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
  if (GetParam() > 0)
    str << "_" << (GetParam()-1)/3 << (GetParam()-1) % 3;
  str << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.createFEMmodel());
}


const std::vector<int> tests = {0,1,2,3,4,5,6,7,8,9,10,11,12};
INSTANTIATE_TEST_CASE_P(TestASMu3D,
                        TestASMu3D,
                        testing::ValuesIn(tests));


TEST(TestASMu3D, TransferGaussPtVars)
{
  SIM3D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  sim.createDefaultModel();
  ASMu3D* pch = static_cast<ASMu3D*>(sim.getPatch(1));
  size_t id[3];
  const double* xi = GaussQuadrature::getCoord(3);
  RealArray oldAr(3*3*3), newAr;
  for (size_t idx = 0; idx < 3; ++idx) {
    SIM3D simNew(1);
    simNew.opt.discretization = ASM::LRSpline;
    simNew.createDefaultModel();
    ASMu3D* pchNew = static_cast<ASMu3D*>(simNew.getPatch(1));
    pchNew->uniformRefine(idx, 1);

    for (id[2] = 0; id[2] < 3; ++id[2])
      for (id[1] = 0; id[1] < 3; ++id[1])
        for (id[0] = 0; id[0] < 3; ++id[0])
          oldAr[id[0]+(id[1]+id[2]*3)*3] = (1.0 + xi[id[idx]]) / 2.0;

    pchNew->transferGaussPtVars(pch->getVolume(), oldAr, newAr, 3);
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
  SIM3D sim(1), sim2(1);
  sim.opt.discretization = sim2.opt.discretization = ASM::LRSpline;
  sim.createDefaultModel();
  sim2.createDefaultModel();
  ASMu3D* pch = static_cast<ASMu3D*>(sim.getPatch(1));
  ASMu3D* pchNew = static_cast<ASMu3D*>(sim2.getPatch(1));
  pchNew->uniformRefine(0, 1);
  RealArray oldAr(3*3*3), newAr;
  std::iota(oldAr.begin(), oldAr.end(), 1);

  pchNew->transferGaussPtVarsN(pch->getVolume(), oldAr, newAr, 3);
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
  EXPECT_EQ(refAr.size(), newAr.size());
  for (size_t i = 0; i < refAr.size(); ++i)
    EXPECT_FLOAT_EQ(refAr[i], newAr[i]);
}


class ASMuCube : public ASMu3D
{
public:
  ASMuCube()
  {
    std::stringstream geo("700 1 0 0\n3 0\n2 2\n0 0 1 1\n2 2\n0 0 1 1\n2 2\n0 0 1 1\n0 0 0\n1 0 0\n0 1 0\n1 1 0\n0 0 1\n1 0 0\n0 1 1\n1 1 1\n");
    EXPECT_TRUE(this->read(geo));
  }
  virtual ~ASMuCube() {}
};


TEST(TestASMu3D, ElementConnectivities)
{
  ASMuCube pch1;
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
