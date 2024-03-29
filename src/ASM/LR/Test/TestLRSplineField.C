//==============================================================================
//!
//! \file TestLRSplineField.C
//!
//! \date Mar 9 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for scalar LR spline fields.
//!
//==============================================================================

#include "Field.h"
#include "FiniteElement.h"
#include "ASMuCube.h"
#include "ASMuSquare.h"

#include "gtest/gtest.h"

#include <memory>


TEST(TestLRSplineField, Value2D)
{
  ASMuSquare patch(1);
  EXPECT_TRUE(patch.generateFEMTopology());

  std::vector<double> sc = {0.0, 1.0, 1.0, 2.0}; // x + y
  std::unique_ptr<Field> fscalar(Field::create(&patch, sc));
  static std::vector<std::array<double,3>> tests_scalar = {{{0.5, 0.5, 1.0}},
                                                           {{1.0, 0.0, 1.0}},
                                                           {{0.0, 1.0, 1.0}},
                                                           {{1.0, 1.0, 2.0}}};
  for (const auto& it : tests_scalar) {
    FiniteElement fe;
    fe.u = it[0];
    fe.v = it[1];
    ASSERT_FLOAT_EQ(fscalar->valueFE(fe), it[2]);
  }
}


TEST(TestLRSplineField, Grad2D)
{
  ASMuSquare patch(1);
  EXPECT_TRUE(patch.generateFEMTopology());

  std::vector<double> sc = {0.0, 1.0, 1.0, 2.0}; // x + y
  std::unique_ptr<Field> fscalar(Field::create(&patch, sc));
  static std::vector<std::array<double,2>> tests_scalar = {{{0.5, 0.5}},
                                                           {{1.0, 0.0}},
                                                           {{0.0, 1.0}},
                                                           {{1.0, 1.0}}};
  for (const auto& it : tests_scalar) {
    FiniteElement fe;
    fe.u = it[0];
    fe.v = it[1];
    Vector v(2);
    fscalar->gradFE(fe, v);
    ASSERT_FLOAT_EQ(v(1), 1.0);
    ASSERT_FLOAT_EQ(v(2), 1.0);
  }
}


TEST(TestLRSplineField, Value3D)
{
  ASMuCube patch(1);
  EXPECT_TRUE(patch.generateFEMTopology());

  std::vector<double> sc = {0.0, 1.0, 1.0, 2.0, 1.0, 2.0, 2.0, 3.0}; // x + y + z
  std::unique_ptr<Field> fscalar(Field::create(&patch, sc));
  static std::vector<std::array<double,4>> tests_scalar = {{{0.5, 0.5, 0.5, 1.5}},
                                                           {{0.0, 0.0, 0.0, 0.0}},
                                                           {{1.0, 0.0, 0.0, 1.0}},
                                                           {{0.0, 1.0, 0.0, 1.0}},
                                                           {{1.0, 1.0, 0.0, 2.0}},
                                                           {{0.0, 0.0, 1.0, 1.0}},
                                                           {{1.0, 0.0, 1.0, 2.0}},
                                                           {{0.0, 1.0, 1.0, 2.0}},
                                                           {{1.0, 1.0, 1.0, 3.0}}};
  for (const auto& it : tests_scalar) {
    FiniteElement fe;
    fe.u = it[0];
    fe.v = it[1];
    fe.w = it[2];
    ASSERT_FLOAT_EQ(fscalar->valueFE(fe), it[3]);
  }
}


TEST(TestLRSplineField, Grad3D)
{
  ASMuCube patch(1);
  EXPECT_TRUE(patch.generateFEMTopology());

  std::vector<double> sc = {0.0, 1.0, 1.0, 2.0, 1.0, 2.0, 2.0, 3.0}; // x + y + z
  std::unique_ptr<Field> fscalar(Field::create(&patch, sc));
  static std::vector<std::array<double,3>> tests_scalar = {{{0.5, 0.5, 0.5}},
                                                           {{0.0, 0.0, 0.0}},
                                                           {{1.0, 0.0, 0.0}},
                                                           {{0.0, 1.0, 0.0}},
                                                           {{1.0, 1.0, 0.0}},
                                                           {{0.0, 0.0, 1.0}},
                                                           {{1.0, 0.0, 1.0}},
                                                           {{0.0, 1.0, 1.0}},
                                                           {{1.0, 1.0, 1.0}}};
  for (const auto& it : tests_scalar) {
    FiniteElement fe;
    fe.u = it[0];
    fe.v = it[1];
    fe.w = it[2];
    Vector v(3);
    fscalar->gradFE(fe, v);
    ASSERT_FLOAT_EQ(v(1), 1.0);
    ASSERT_FLOAT_EQ(v(2), 1.0);
    ASSERT_FLOAT_EQ(v(3), 1.0);
  }
}
