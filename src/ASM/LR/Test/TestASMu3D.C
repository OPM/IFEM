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


static void getBoundaryNodes (const char* faceName, std::vector<int>& nodes)
{
  SIM3D sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("src/ASM/LR/Test/refdata/boundary_nodes_3d.xinp"));
  ASSERT_TRUE(sim.createFEMmodel());
  sim.getBoundaryNodes(sim.getUniquePropertyCode(faceName,0),nodes);
}


TEST(TestASMu3D, BoundaryNodesF1)
{
  std::vector<int> vec;
  getBoundaryNodes("Face1",vec);
  ASSERT_EQ(vec.size(), 16U);
  for (int i = 0; i < 16; ++i)
    ASSERT_EQ(vec[i], 1+4*i);
}


TEST(TestASMu3D, BoundaryNodesF2)
{
  std::vector<int> vec;
  getBoundaryNodes("Face2",vec);
  ASSERT_EQ(vec.size(), 16U);
  for (int i = 0; i < 16; ++i)
    ASSERT_EQ(vec[i], 2+4*i);
}


TEST(TestASMu3D, BoundaryNodesF3)
{
  std::vector<int> vec;
  getBoundaryNodes("Face3",vec);
  ASSERT_EQ(vec.size(), 16U);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      ASSERT_EQ(vec[4*i+j], 16*i+j+1);
}


TEST(TestASMu3D, BoundaryNodesF4)
{
  std::vector<int> vec;
  getBoundaryNodes("Face4",vec);
  ASSERT_EQ(vec.size(), 16U);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      ASSERT_EQ(vec[4*i+j], 5+16*i+j);
}


TEST(TestASMu3D, BoundaryNodesF5)
{
  std::vector<int> vec;
  getBoundaryNodes("Face5",vec);
  ASSERT_EQ(vec.size(), 16U);
  for (int i = 0; i < 16; ++i)
    ASSERT_EQ(vec[i], i+1);
}


TEST(TestASMu3D, BoundaryNodesF6)
{
  std::vector<int> vec;
  getBoundaryNodes("Face6",vec);
  ASSERT_EQ(vec.size(), 16U);
  for (int i = 0; i < 16; ++i)
    ASSERT_EQ(vec[i], 17+i);
}
