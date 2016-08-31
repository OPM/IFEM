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
#include "IFEM.h"
#include "SIM3D.h"

#include "gtest/gtest.h"

TEST(TestASMu3D, BoundaryNodesF1)
{
  SIM3D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes_3d.xinp"));
  sim.preprocess();

  int bcode = sim.getUniquePropertyCode("Face1",0);
  std::vector<int> vec;
  sim.getBoundaryNodes(bcode, vec);

  ASSERT_EQ(vec.size(), 16U);
  for (int i = 0; i < 16; ++i)
    ASSERT_EQ(vec[i], 1+4*i);
}


TEST(TestASMu3D, BoundaryNodesF2)
{
  SIM3D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes_3d.xinp"));
  sim.preprocess();

  int bcode = sim.getUniquePropertyCode("Face2",0);
  std::vector<int> vec;
  sim.getBoundaryNodes(bcode, vec);
  ASSERT_EQ(vec.size(), 16U);
  for (int i = 0; i < 16; ++i)
    ASSERT_EQ(vec[i], 2+4*i);
}


TEST(TestASMu3D, BoundaryNodesF3)
{
  SIM3D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes_3d.xinp"));
  sim.preprocess();

  int bcode = sim.getUniquePropertyCode("Face3",0);
  std::vector<int> vec;
  sim.getBoundaryNodes(bcode, vec);

  ASSERT_EQ(vec.size(), 16U);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
    ASSERT_EQ(vec[4*i+j], 16*i+j+1);
}


TEST(TestASMu3D, BoundaryNodesF4)
{
  SIM3D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes_3d.xinp"));
  sim.preprocess();

  int bcode = sim.getUniquePropertyCode("Face4",0);
  std::vector<int> vec;
  sim.getBoundaryNodes(bcode, vec);
  ASSERT_EQ(vec.size(), 16U);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
    ASSERT_EQ(vec[4*i+j], 5+16*i+j);
}


TEST(TestASMu3D, BoundaryNodesF5)
{
  SIM3D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes_3d.xinp"));
  sim.preprocess();

  int bcode = sim.getUniquePropertyCode("Face5",0);
  std::vector<int> vec;
  sim.getBoundaryNodes(bcode, vec);
  ASSERT_EQ(vec.size(), 16U);
  for (int i = 0; i < 16; ++i)
    ASSERT_EQ(vec[i], i+1);
}


TEST(TestASMu3D, BoundaryNodesF6)
{
  SIM3D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes_3d.xinp"));
  sim.preprocess();

  int bcode = sim.getUniquePropertyCode("Face6",0);
  std::vector<int> vec;
  sim.getBoundaryNodes(bcode, vec);
  ASSERT_EQ(vec.size(), 16U);
  for (int i = 0; i < 16; ++i)
    ASSERT_EQ(vec[i], 17+i);
}
