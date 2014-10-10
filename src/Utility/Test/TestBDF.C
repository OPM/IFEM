//==============================================================================
//!
//! \file TestBDF.C
//!
//! \date Oct 10 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for helper functions for BDF based time stepping.
//!
//==============================================================================

#include "BDF.h"

#include "gtest/gtest.h"

TEST(TestBDF, BDF_1)
{
  TimeIntegration::BDF bdf(1);

  bdf.advanceStep();
  ASSERT_EQ(bdf.getOrder(), 1);
  ASSERT_EQ(bdf.getActualOrder(), 1);
  ASSERT_FLOAT_EQ(bdf[0], 1.0);
  ASSERT_FLOAT_EQ(bdf[1], -1.0);
  bdf.advanceStep();
  ASSERT_EQ(bdf.getOrder(), 1);
  ASSERT_FLOAT_EQ(bdf[0], 1.0);
  ASSERT_FLOAT_EQ(bdf[1], -1.0);
}

TEST(TestBDF, BDF_2)
{
  TimeIntegration::BDF bdf(2);

  bdf.advanceStep();
  ASSERT_EQ(bdf.getOrder(), 1);
  ASSERT_EQ(bdf.getActualOrder(), 2);
  ASSERT_FLOAT_EQ(bdf[0], 1.0);
  ASSERT_FLOAT_EQ(bdf[1], -1.0);
  bdf.advanceStep();
  ASSERT_EQ(bdf.getOrder(), 2);
  ASSERT_FLOAT_EQ(bdf[0], 1.5);
  ASSERT_FLOAT_EQ(bdf[1], -2.0);
  ASSERT_FLOAT_EQ(bdf[2], 0.5);
}

TEST(TestBDF, BDFD2_1)
{
  TimeIntegration::BDFD2 bdf(1);

  bdf.advanceStep(0.1, 0.1);
  ASSERT_EQ(bdf.getOrder(), 1);
  ASSERT_EQ(bdf.getActualOrder(), 1);
  ASSERT_FLOAT_EQ(bdf[0], 2.0);
  ASSERT_FLOAT_EQ(bdf[1], -2.0);
  bdf.advanceStep(0.1, 0.1);
  ASSERT_EQ(bdf.getOrder(), 1);
  ASSERT_FLOAT_EQ(bdf[0], 1.0);
  ASSERT_FLOAT_EQ(bdf[1], -2.0);
  ASSERT_FLOAT_EQ(bdf[2], 1.0);
}

TEST(TestBDF, BDFD2_2)
{
  TimeIntegration::BDFD2 bdf(2);

  bdf.advanceStep(0.1, 0.1);
  ASSERT_EQ(bdf.getOrder(), 1);
  ASSERT_EQ(bdf.getActualOrder(), 2);
  ASSERT_FLOAT_EQ(bdf[0], 2.0);
  ASSERT_FLOAT_EQ(bdf[1], -2.0);
  bdf.advanceStep(0.1, 0.1);
  ASSERT_EQ(bdf.getOrder(), 1);
  ASSERT_FLOAT_EQ(bdf[0], 2.5);
  ASSERT_FLOAT_EQ(bdf[1], -8.0);
  ASSERT_FLOAT_EQ(bdf[2], 5.5);
  bdf.advanceStep(0.1, 0.1);
  ASSERT_EQ(bdf.getOrder(), 2);
  ASSERT_FLOAT_EQ(bdf[0], 2.0);
  ASSERT_FLOAT_EQ(bdf[1], -5.0);
  ASSERT_FLOAT_EQ(bdf[2], 4.0);
  ASSERT_FLOAT_EQ(bdf[3], -1.0);
}
