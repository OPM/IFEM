//==============================================================================
//!
//! \file TestSAWallLaw.C
//!
//! \date Oct 13 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for various a SA parametrization of a turbulent
//!        boundary layer.
//!
//==============================================================================

#include "SAWallLaw.h"

#include "gtest/gtest.h"

TEST(TestSAWallLaw, ComputeYplus)
{
  SAWallLaw sa;

  double result1;
  ASSERT_TRUE(sa.computeYplus(0.01, 0.3, 0.3, result1));
  ASSERT_FLOAT_EQ(result1, 0.1000000000775523);
}

TEST(TestSAWallLaw, ComputeUstar)
{
  SAWallLaw sa;

  double result1;
  ASSERT_TRUE(sa.computeUstar(0.01, 0.3, 0.3, result1));
  ASSERT_FLOAT_EQ(result1, 3.000000002410403);
}

TEST(TestSAWallLaw, ComputeTauB)
{
  SAWallLaw sa;

  double result1;
  ASSERT_TRUE(sa.computeTauB(0.01, 0.3, 0.3, result1));
  ASSERT_FLOAT_EQ(result1, 30.00000004820805);
}

TEST(TestSAWallLaw, ComputeNu)
{
  SAWallLaw sa;

  double result1;
  ASSERT_TRUE(sa.computeNu(0.01, 0.3, 0.3, result1));
  ASSERT_FLOAT_EQ(result1, 0.01230000000988265);
}
