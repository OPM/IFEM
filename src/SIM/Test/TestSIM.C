//==============================================================================
//!
//! \file TestSIM.C
//!
//! \date Apr 27 2016
//!
//! \author Arne Morten Kvarving / NTNU
//!
//! \brief Tests for various expected SIM behavior.
//!
//==============================================================================

#include "SIM2D.h"

#include "gtest/gtest.h"

TEST(TestSIM, UniqueBoundaryNodes)
{
  SIM2D sim(1);
  ASSERT_TRUE(sim.read("src/SIM/Test/refdata/boundary_nodes.xinp"));
  sim.preprocess();

  int bcode = sim.getUniquePropertyCode("dir",0);
  std::vector<int> vec;
  sim.getBoundaryNodes(bcode, vec);

  std::sort(vec.begin(), vec.end());
  ASSERT_TRUE(std::unique(vec.begin(), vec.end()) == vec.end());
}
