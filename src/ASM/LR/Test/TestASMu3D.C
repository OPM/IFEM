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
