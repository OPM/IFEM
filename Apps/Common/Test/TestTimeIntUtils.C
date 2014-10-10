//==============================================================================
//!
//! \file TestTimeIntUtils.C
//!
//! \date Oct 7 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for various helpers for time integration.
//!
//==============================================================================

#include "TimeIntUtils.h"

#include "gtest/gtest.h"

using namespace TimeIntegration;

TEST(TestTimeIntUtils, Order)
{
  EXPECT_EQ(Order(EULER), 1);
  EXPECT_EQ(Order(BE),    1);
  EXPECT_EQ(Order(HEUN),  2);
  EXPECT_EQ(Order(BDF2),  2);
  EXPECT_EQ(Order(RK3),   3);
  EXPECT_EQ(Order(RK4),   4);
}

TEST(TestTimeIntUtils, Steps)
{
  EXPECT_EQ(Steps(EULER), 1);
  EXPECT_EQ(Steps(BE),    1);
  EXPECT_EQ(Steps(HEUN),  1);
  EXPECT_EQ(Steps(BDF2),  2);
  EXPECT_EQ(Steps(RK3),   1);
  EXPECT_EQ(Steps(RK4),   1);
}
