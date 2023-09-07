//==============================================================================
//!
//! \file TestSplineFields.C
//!
//! \date Oct 6 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for vector spline fields.
//!
//==============================================================================

#include "Field.h"
#include "Fields.h"
#include "ItgPoint.h"
#include "ASMSquare.h"
#include "ASMCube.h"

#include "gtest/gtest.h"
#include <array>
#include <memory>


TEST(TestSplineFields, Value2D)
{
  ASMSquare patch;

  // {x+y+x*y, x-y+x*y}
  std::vector<double> vc = {0.0,  0.0,
                            1.0,  1.0,
                            1.0, -1.0,
                            3.0, 1.0};
  std::unique_ptr<Fields> fvector(Fields::create(&patch,vc));
  std::unique_ptr<Field> fscalar(Field::create(&patch,vc,1,2));
  static std::vector<std::array<double,4>> tests_vector =
        {{{{0.5, 0.5, 1.25, 0.25}},
          {{1.0, 0.0, 1.0,  1.0}},
          {{0.0, 1.0, 1.0, -1.0}},
          {{1.0, 1.0, 3.0,  1.0}}}};
  for (const auto& it : tests_vector) {
    ItgPoint fe(it[0],it[1]);
    Vector v(2);
    ASSERT_TRUE(fvector->valueFE(fe,v));
    EXPECT_FLOAT_EQ(v(1), it[2]);
    EXPECT_FLOAT_EQ(v(2), it[3]);
    EXPECT_FLOAT_EQ(fscalar->valueFE(fe), it[3]);
  }
}


TEST(TestSplineFields, Value2Dmx)
{
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  ASMmxSquare patch({1,1,1});
  patch.raiseOrder(1,1);
  EXPECT_TRUE(patch.generateFEMTopology());

  // {x+y+x*y, x-y+x*y}
  std::vector<double> vc = {0.0,
                            0.5,
                            1.0,
                            1.0,
                            2.0,
                            3.0,
                            // y
                            0.0,
                            1.0,
                           -0.5,
                            1.0,
                           -1.0,
                            1.0,
                            // p
                            0.0, 1.0, 1.0, 2.0}; // x + y
  std::unique_ptr<Fields> fvector(Fields::create(&patch,vc,12));
  std::unique_ptr<Field> fscalar(Field::create(&patch,vc,3,1));
  static std::vector<std::array<double,5>> tests_vector =
                          {{{{0.5, 0.5, 1.25, 0.25, 1.0}},
                            {{1.0, 0.0, 1.0,  1.0,  1.0}},
                            {{0.0, 1.0, 1.0, -1.0,  1.0}},
                            {{1.0, 1.0, 3.0,  1.0,  2.0}}}};
  for (const auto& it : tests_vector) {
    Vector v(2);
    ItgPoint fe(it[0],it[1]);
    ASSERT_TRUE(fvector->valueFE(fe,v));
    EXPECT_FLOAT_EQ(v(1), it[2]);
    EXPECT_FLOAT_EQ(v(2), it[3]);
    EXPECT_FLOAT_EQ(fscalar->valueFE(fe), it[4]);
  }
}


TEST(TestSplineFields, Grad2D)
{
  ASMSquare patch;

  // {x+y+x*y, x-y+x*y}
  std::vector<double> vc = {0.0, 0.0, 1.0, 1.0, 1.0, -1.0, 3.0, 1.0};
  std::unique_ptr<Fields> fvector(Fields::create(&patch,vc));
  static std::vector<std::array<double,6>> tests_vector =
        {{{{0.5, 0.5, 1.5, 1.5, 1.5, -0.5}},
          {{1.0, 0.0, 1.0, 2.0, 1.0, 0.0}},
          {{0.0, 1.0, 2.0, 1.0, 2.0, -1.0}},
          {{1.0, 1.0, 2.0, 2.0, 2.0, 0.0}}}};
  for (const auto& it : tests_vector) {
    Matrix gradu(2,2);
    ASSERT_TRUE(fvector->gradFE(ItgPoint(it[0],it[1]),gradu));
    EXPECT_FLOAT_EQ(gradu(1,1), it[2]);
    EXPECT_FLOAT_EQ(gradu(1,2), it[3]);
    EXPECT_FLOAT_EQ(gradu(2,1), it[4]);
    EXPECT_FLOAT_EQ(gradu(2,2), it[5]);
  }
}


TEST(TestSplineFields, Value3D)
{
  ASMCube patch;

  // {x+y+z, x+y-z, x-y+z}
  std::vector<double> vc = {0.0,  0.0,  0.0,
                            1.0,  1.0,  1.0,
                            1.0,  1.0, -1.0,
                            2.0,  2.0,  0.0,
                            1.0, -1.0,  1.0,
                            2.0,  0.0,  2.0,
                            2.0,  0.0,  0.0,
                            3.0,  1.0,  1.0};
  std::unique_ptr<Fields> fvector(Fields::create(&patch,vc));
  std::unique_ptr<Field> fscalar(Field::create(&patch,vc,1,2));
  static std::vector<std::array<double,6>> tests_scalar =
      {{{{0.5, 0.5, 0.5, 1.5,  0.5,  0.5}},
        {{0.0, 0.0, 0.0, 0.0,  0.0,  0.0}},
        {{1.0, 0.0, 0.0, 1.0,  1.0,  1.0}},
        {{0.0, 1.0, 0.0, 1.0,  1.0, -1.0}},
        {{1.0, 1.0, 0.0, 2.0,  2.0,  0.0}},
        {{0.0, 0.0, 1.0, 1.0, -1.0,  1.0}},
        {{1.0, 0.0, 1.0, 2.0,  0.0,  2.0}},
        {{0.0, 1.0, 1.0, 2.0,  0.0,  0.0}},
        {{1.0, 1.0, 1.0, 3.0,  1.0,  1.0}}}};
  for (const auto& it : tests_scalar) {
    ItgPoint fe(it.data());
    Vector v(3);
    ASSERT_TRUE(fvector->valueFE(fe,v));
    EXPECT_FLOAT_EQ(v(1), it[3]);
    EXPECT_FLOAT_EQ(v(2), it[4]);
    EXPECT_FLOAT_EQ(v(3), it[5]);
    EXPECT_FLOAT_EQ(fscalar->valueFE(fe), it[4]);
  }
}


TEST(TestSplineFields, Grad3D)
{
  ASMCube patch;

  // {x+y+z+x*y*z, x+y-z+x*y*z, x-y+z+x*y*z}
  std::vector<double> vc = {0.0,  0.0,  0.0,
                            1.0,  1.0,  1.0,
                            1.0,  1.0, -1.0,
                            2.0,  2.0,  0.0,
                            1.0, -1.0,  1.0,
                            2.0,  0.0,  2.0,
                            2.0,  0.0,  0.0,
                            4.0,  2.0,  2.0};
  std::unique_ptr<Fields> fvector(Fields::create(&patch,vc));
  static std::vector<std::pair<std::array<double,3>,
                               std::array<double,9>>> tests_vector =
    {{{{{0.5, 0.5, 0.5}}, {{1.25,  1.25,  1.25,
                            1.25,  1.25, -0.75,
                            1.25, -0.75,  1.25}}},
      {{{0.0, 0.0, 0.0}}, {{1.0,  1.0,  1.0,
                            1.0,  1.0, -1.0,
                            1.0, -1.0,  1.0}}},
      {{{1.0, 0.0, 0.0}}, {{1.0,  1.0,  1.0,
                            1.0,  1.0, -1.0,
                            1.0, -1.0,  1.0}}},
      {{{0.0, 1.0, 0.0}}, {{1.0,  1.0,  1.0,
                            1.0,  1.0, -1.0,
                            1.0, -1.0,  1.0}}},
      {{{1.0, 1.0, 0.0}}, {{1.0,  1.0,  2.0,
                            1.0,  1.0,  0.0,
                            1.0, -1.0,  2.0}}},
      {{{0.0, 0.0, 1.0}}, {{1.0,  1.0,  1.0,
                            1.0,  1.0, -1.0,
                            1.0, -1.0,  1.0}}},
      {{{1.0, 0.0, 1.0}}, {{1.0,  2.0,  1.0,
                            1.0,  2.0, -1.0,
                            1.0,  0.0,  1.0}}},
      {{{0.0, 1.0, 1.0}}, {{2.0,  1.0,  1.0,
                            2.0,  1.0, -1.0,
                            2.0, -1.0,  1.0}}},
      {{{1.0, 1.0, 1.0}}, {{2.0,  2.0,  2.0,
                            2.0,  2.0,  0.0,
                            2.0,  0.0,  2.0}}}}};
  for (const auto& it : tests_vector) {
    Matrix gradu(3,3);
    ASSERT_TRUE(fvector->gradFE(ItgPoint(it.first.data()),gradu));
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j <3; ++j)
        EXPECT_FLOAT_EQ(gradu(i+1,j+1), it.second[i*3+j]);
  }
}


TEST(TestSplineFields, Value3Dmx)
{
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  ASMmxCube patch({1,1,1,1});
  patch.raiseOrder(1,1,1);
  EXPECT_TRUE(patch.generateFEMTopology());

  // {x+y+z+x*y*z, x+y-z+x*y*z, x-y+z+x*y*z}
  std::vector<double> vc = {0.0, 0.5, 1.0,
                            1.0, 1.5, 2.0,
                            1.0, 1.5, 2.0,
                            2.0, 3.0, 4.0,
                            // y
                            0.0, 1.0,
                            0.5, 1.5,
                            1.0, 2.0,
                           -1.0, 0.0,
                           -0.5, 1.0,
                            0.0, 2.0,
                            // z
                            0.0, 1.0,
                           -1.0, 0.0,
                            0.5, 1.5,
                           -0.5, 1.0,
                            1.0, 2.0,
                            0.0, 2.0,
                            // p
                            0.0, 1.0, 1.0, 2.0, 1.0, 2.0, 2.0, 3.0};
  std::unique_ptr<Fields> fvector(Fields::create(&patch,vc,123));
  std::unique_ptr<Field> fscalar(Field::create(&patch,vc,4,1));
  static std::vector<std::array<double,7>> tests_vector =
                         {{{{0.5, 0.5, 0.5, 1.625, 0.625, 0.625, 1.5}},
                           {{0.0, 0.0, 0.0,   0.0,   0.0,   0.0, 0.0}},
                           {{1.0, 0.0, 0.0,   1.0,   1.0,   1.0, 1.0}},
                           {{0.0, 1.0, 0.0,   1.0,   1.0,  -1.0, 1.0}},
                           {{1.0, 1.0, 0.0,   2.0,   2.0,   0.0, 2.0}},
                           {{0.0, 0.0, 1.0,   1.0,  -1.0,   1.0, 1.0}},
                           {{1.0, 0.0, 1.0,   2.0,   0.0,   2.0, 2.0}},
                           {{0.0, 1.0, 1.0,   2.0,   0.0,   0.0, 2.0}},
                           {{1.0, 1.0, 1.0,   4.0,   2.0,   2.0, 3.0}}}};
  for (const auto& it : tests_vector) {
    ItgPoint fe(it.data());
    Vector v(3);
    ASSERT_TRUE(fvector->valueFE(fe,v));
    EXPECT_FLOAT_EQ(v(1), it[3]);
    EXPECT_FLOAT_EQ(v(2), it[4]);
    EXPECT_FLOAT_EQ(v(3), it[5]);
    EXPECT_FLOAT_EQ(fscalar->valueFE(fe), it[6]);
  }
}
