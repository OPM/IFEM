//==============================================================================
//!
//! \file TestSpalding.C
//!
//! \date Oct 13 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for various a spalding parametrization of a turbulent
//!        boundary layer.
//!
//==============================================================================

#include "Spalding.h"

#include "gtest/gtest.h"

TEST(TestSpalding, ComputeTauB)
{
  Spalding spald;

  double result1;
  ASSERT_TRUE(spald.computeTauB(1.0, 1.0, 0.3, 0.3, result1));
  ASSERT_FLOAT_EQ(result1, 0.300038486432162);
}
