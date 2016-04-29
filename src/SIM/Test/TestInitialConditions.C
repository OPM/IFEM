//==============================================================================
//!
//! \file TestInitialConditions.C
//!
//! \date Apr 6 2016
//!
//! \author Timo van Opstal / NTNU
//!
//! \brief Tests initial condition handling.
//!
//==============================================================================

#include "SIM2D.h"

#include "gtest/gtest.h"

TEST(TestInitialConditions, Parse)
{
  SIM2D sim({4}, false);
  EXPECT_TRUE(sim.read("src/SIM/Test/refdata/input.xinp"));

  // Recognize both comp and component attributes and correct priority
  // Boundary conditions
  RealFunc* func;
  const Vec3 X;
  for (int i = 1; i < 5; i++) {
    func = sim.getSclFunc(i);
    ASSERT_FLOAT_EQ((float) i, (*func)(X));
  }
  // Initial conditions
  const SIMdependency::InitialCondMap& ic = sim.getICs();
  SIMdependency::InitialCondMap::const_iterator foo = ic.find("nofile");
  std::vector<SIMdependency::ICInfo> bar = foo->second;
  for (std::vector<SIMdependency::ICInfo>::iterator info = bar.begin();
       info != bar.end(); info++)
    switch (info->component) {
      case 1:
        ASSERT_EQ(info->function, "1");
        break;
      case 2:
        ASSERT_EQ(info->function, "2");
        break;
      case 3:
        ASSERT_EQ(info->function, "3");
        break;
      case 4:
        ASSERT_EQ(info->function, "4");
        break;
      default:
        ASSERT_TRUE( false );
      }
}
