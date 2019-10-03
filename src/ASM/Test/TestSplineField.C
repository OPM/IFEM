//==============================================================================
//!
//! \file TestSplineField.C
//!
//! \date Oct 6 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for scalar spline fields. 
//!
//==============================================================================

#include "Field.h"
#include "ItgPoint.h"
#include "ASMSquare.h"
#include "ASMCube.h"

#include "gtest/gtest.h"
#include <array>


TEST(TestSplineField, Value2D)
{
  ASMSquare patch(1);
  std::vector<double> sc = {0.0, 1.0, 1.0, 2.0}; // x + y
  Field* fscalar = Field::create(&patch,sc);
  static std::vector<std::array<double,3>> tests_scalar = {{{{0.5, 0.5, 1.0}},
                                                            {{1.0, 0.0, 1.0}},
                                                            {{0.0, 1.0, 1.0}},
                                                            {{1.0, 1.0, 2.0}}}};
  for (const auto& it : tests_scalar)
    EXPECT_FLOAT_EQ(fscalar->valueFE(ItgPoint(it[0],it[1])),it[2]);
}

TEST(TestSplineField, Grad2D)
{
  ASMSquare patch(1);
  std::vector<double> sc = {0.0, 1.0, 1.0, 2.0}; // x + y
  Field* fscalar = Field::create(&patch,sc);
  static std::vector<std::array<double,2>> tests_scalar = {{{{0.5, 0.5}},
                                                            {{1.0, 0.0}},
                                                            {{0.0, 1.0}},
                                                            {{1.0, 1.0}}}};
  for (const auto& it : tests_scalar) {
    Vector v(2);
    ASSERT_TRUE(fscalar->gradFE(ItgPoint(it[0],it[1]),v));
    EXPECT_FLOAT_EQ(v(1),1.0);
    EXPECT_FLOAT_EQ(v(2),1.0);
  }
}


TEST(TestSplineField, Value3D)
{
  ASMCube patch(1);
  std::vector<double> sc = {0.0, 1.0, 1.0, 2.0, 1.0, 2.0, 2.0, 3.0}; // x + y + z
  Field* fscalar = Field::create(&patch,sc);
  static std::vector<std::array<double,4>> tests_scalar = {{{{0.5, 0.5, 0.5, 1.5}},
                                                            {{0.0, 0.0, 0.0, 0.0}},
                                                            {{1.0, 0.0, 0.0, 1.0}},
                                                            {{0.0, 1.0, 0.0, 1.0}},
                                                            {{1.0, 1.0, 0.0, 2.0}},
                                                            {{0.0, 0.0, 1.0, 1.0}},
                                                            {{1.0, 0.0, 1.0, 2.0}},
                                                            {{0.0, 1.0, 1.0, 2.0}},
                                                            {{1.0, 1.0, 1.0, 3.0}}}};
  for (const auto& it : tests_scalar)
    EXPECT_FLOAT_EQ(fscalar->valueFE(ItgPoint(it[0],it[1],it[2])),it[3]);
}

TEST(TestSplineField, Grad3D)
{
  ASMCube patch(1);
  std::vector<double> sc = {0.0, 1.0, 1.0, 2.0, 1.0, 2.0, 2.0, 3.0}; // x + y + z
  Field* fscalar = Field::create(&patch,sc);
  static std::vector<std::array<double,3>> tests_scalar = {{{{0.5, 0.5, 0.5}},
                                                            {{0.0, 0.0, 0.0}},
                                                            {{1.0, 0.0, 0.0}},
                                                            {{0.0, 1.0, 0.0}},
                                                            {{1.0, 1.0, 0.0}},
                                                            {{0.0, 0.0, 1.0}},
                                                            {{1.0, 0.0, 1.0}},
                                                            {{0.0, 1.0, 1.0}},
                                                            {{1.0, 1.0, 1.0}}}};
  for (const auto& it : tests_scalar) {
    Vector v(3);
    ASSERT_TRUE(fscalar->gradFE(ItgPoint(it[0],it[1],it[2]),v));
    EXPECT_FLOAT_EQ(v(1),1.0);
    EXPECT_FLOAT_EQ(v(2),1.0);
    EXPECT_FLOAT_EQ(v(3),1.0);
  }
}
