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
#include "IFEM.h"
#include "SIM2D.h"

#include "gtest/gtest.h"

TEST(TestASMu2D, BoundaryNodesE1)
{
  SIM2D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes.xinp"));
  sim.preprocess();

  int bcode = sim.getUniquePropertyCode("Edge1",0);
  std::vector<int> vec;
  sim.getBoundaryNodes(bcode, vec);
  ASSERT_EQ(vec.size(), 4U);
  for (int i = 0; i < 4; ++i)
    ASSERT_EQ(vec[i], 4*i+1);
}


TEST(TestASMu2D, BoundaryNodesE2)
{
  SIM2D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes.xinp"));
  sim.preprocess();

  int bcode = sim.getUniquePropertyCode("Edge2",0);
  std::vector<int> vec;
  sim.getBoundaryNodes(bcode, vec);
  ASSERT_EQ(vec.size(), 4U);
  for (int i = 0; i < 4; ++i)
    ASSERT_EQ(vec[i], 2+4*i);
}


TEST(TestASMu2D, BoundaryNodesE3)
{
  SIM2D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes.xinp"));
  sim.preprocess();

  int bcode = sim.getUniquePropertyCode("Edge3",0);
  std::vector<int> vec;
  sim.getBoundaryNodes(bcode, vec);
  ASSERT_EQ(vec.size(), 4U);
  for (int i = 0; i < 4; ++i)
    ASSERT_EQ(vec[i], i+1);
}


TEST(TestASMu2D, BoundaryNodesE4)
{
  SIM2D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes.xinp"));
  sim.preprocess();

  int bcode = sim.getUniquePropertyCode("Edge4",0);
  std::vector<int> vec;
  sim.getBoundaryNodes(bcode, vec);
  ASSERT_EQ(vec.size(), 4U);
  for (int i = 0; i < 4; ++i)
    ASSERT_EQ(vec[i], 5+i);
}
