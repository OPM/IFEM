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
#include "GaussQuadrature.h"
#include "SIM3D.h"
#include "LRSpline/LRSplineVolume.h"

#include "gtest/gtest.h"

class TestASMu3D :
  public testing::Test,
  public testing::WithParamInterface<int>
{
};


TEST_P(TestASMu3D, BoundaryNodes)
{
  SIM3D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes_3d.xinp"));
  sim.preprocess();

  std::stringstream str;
  str << "Face" << GetParam();
  int bcode = sim.getUniquePropertyCode(str.str(),0);
  std::vector<int> vec;
  sim.getBoundaryNodes(bcode,vec);
  ASSERT_EQ(vec.size(), 16U);
  auto it = vec.begin();
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      if (GetParam() == 1)
        ASSERT_EQ(*it++, 1+4*(4*i+j));
      else if (GetParam() == 2)
        ASSERT_EQ(*it++, 2+4*(4*i+j));
      else if (GetParam() == 3)
        ASSERT_EQ(*it++, 16*i+j+1);
      else if (GetParam() == 4)
        ASSERT_EQ(*it++, 5+16*i+j);
      else if (GetParam() == 5)
        ASSERT_EQ(*it++, 4*i+j+1);
      else if (GetParam() == 6)
        ASSERT_EQ(*it++, 17+4*i+j);
    }
  }
}


const std::vector<int> tests = {1,2,3,4,5,6};
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
