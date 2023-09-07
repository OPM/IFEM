//==============================================================================
//!
//! \file TestLRSplineFields.C
//!
//! \date Mar 9 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for vector LR spline fields.
//!
//==============================================================================

#include "Field.h"
#include "Fields.h"
#include "FiniteElement.h"
#include "ASMmxBase.h"
#include "ASMuCube.h"
#include "ASMuSquare.h"

#include "gtest/gtest.h"
#include <array>
#include <memory>
#include <vector>


TEST(TestLRSplineFields, Value2D)
{
  ASMuSquare patch;
  EXPECT_TRUE(patch.generateFEMTopology());

  // {x+y+x*y, x-y+x*y}
  std::vector<double> vc = {0.0,  0.0,
                            1.0,  1.0,
                            1.0, -1.0,
                            3.0, 1.0};
  std::unique_ptr<Fields> fvector(Fields::create(&patch, vc));
  std::unique_ptr<Field> fscalar(Field::create(&patch, vc, 1, 2));
  static std::vector<std::array<double,4>> tests_vector =
        {{{0.5, 0.5, 1.25, 0.25}},
         {{1.0, 0.0, 1.0,  1.0}},
         {{0.0, 1.0, 1.0, -1.0}},
         {{1.0, 1.0, 3.0,  1.0}}};
  for (const auto& it : tests_vector) {
    FiniteElement fe;
    fe.u = it[0];
    fe.v = it[1];
    Vector v(2);
    fvector->valueFE(fe, v);
    EXPECT_FLOAT_EQ(v(1), it[2]);
    EXPECT_FLOAT_EQ(v(2), it[3]);
    EXPECT_FLOAT_EQ(fscalar->valueFE(fe), it[3]);
  }
}


TEST(TestLRSplineFields, Value2Dmx)
{
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  ASMmxuSquare patch({1,1,1});
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
  std::unique_ptr<Fields> fvector(Fields::create(&patch, vc, 12));
  std::unique_ptr<Field> fscalar(Field::create(&patch, vc, 3, 1));
  static std::vector<std::array<double,5>> tests_vector =
                          {{{0.5, 0.5, 1.25, 0.25, 1.0}},
                           {{1.0, 0.0, 1.0,  1.0,  1.0}},
                           {{0.0, 1.0, 1.0, -1.0,  1.0}},
                           {{1.0, 1.0, 3.0,  1.0,  2.0}}};
  for (const auto& it : tests_vector) {
    FiniteElement fe;
    fe.u = it[0];
    fe.v = it[1];
    Vector v(2);
    fvector->valueFE(fe, v);
    EXPECT_FLOAT_EQ(v(1), it[2]);
    EXPECT_FLOAT_EQ(v(2), it[3]);
    EXPECT_FLOAT_EQ(fscalar->valueFE(fe), it[4]);
  }
}


TEST(TestLRSplineFields, Grad2D)
{
  ASMuSquare patch;
  EXPECT_TRUE(patch.generateFEMTopology());

  // {x+y+x*y, x-y+x*y}
  std::vector<double> vc = {0.0, 0.0, 1.0, 1.0, 1.0, -1.0, 3.0, 1.0};
  std::unique_ptr<Fields> fvector(Fields::create(&patch, vc));
  static std::vector<std::array<double,6>> tests_vector =
        {{{0.5, 0.5, 1.5, 1.5, 1.5, -0.5}},
         {{1.0, 0.0, 1.0, 2.0, 1.0, 0.0}},
         {{0.0, 1.0, 2.0, 1.0, 2.0, -1.0}},
         {{1.0, 1.0, 2.0, 2.0, 2.0, 0.0}}};
  for (const auto& it : tests_vector) {
    FiniteElement fe;
    fe.u = it[0];
    fe.v = it[1];
    Matrix gradu(2,2);
    fvector->gradFE(fe, gradu);
    EXPECT_FLOAT_EQ(gradu(1,1), it[2]);
    EXPECT_FLOAT_EQ(gradu(1,2), it[3]);
    EXPECT_FLOAT_EQ(gradu(2,1), it[4]);
    EXPECT_FLOAT_EQ(gradu(2,2), it[5]);
  }
}


TEST(TestLRSplineFields, Value3D)
{
  ASMuCube patch;
  EXPECT_TRUE(patch.generateFEMTopology());

  // {x+y+z, x+y-z, x-y+z}
  std::vector<double> vc = {0.0,  0.0,  0.0,
                            1.0,  1.0,  1.0,
                            1.0,  1.0, -1.0,
                            2.0,  2.0,  0.0,
                            1.0, -1.0,  1.0,
                            2.0,  0.0,  2.0,
                            2.0,  0.0,  0.0,
                            3.0,  1.0,  1.0};
  std::unique_ptr<Fields> fvector(Fields::create(&patch, vc));
  std::unique_ptr<Field> fscalar(Field::create(&patch, vc, 1, 2));
  static std::vector<std::array<double,6>> tests_scalar =
      {{{0.5, 0.5, 0.5, 1.5,  0.5,  0.5}},
       {{0.0, 0.0, 0.0, 0.0,  0.0,  0.0}},
       {{1.0, 0.0, 0.0, 1.0,  1.0,  1.0}},
       {{0.0, 1.0, 0.0, 1.0,  1.0, -1.0}},
       {{1.0, 1.0, 0.0, 2.0,  2.0,  0.0}},
       {{0.0, 0.0, 1.0, 1.0, -1.0,  1.0}},
       {{1.0, 0.0, 1.0, 2.0,  0.0,  2.0}},
       {{0.0, 1.0, 1.0, 2.0,  0.0,  0.0}},
       {{1.0, 1.0, 1.0, 3.0,  1.0,  1.0}}};
  for (const auto& it : tests_scalar) {
    FiniteElement fe;
    fe.u = it[0];
    fe.v = it[1];
    fe.w = it[2];
    Vector v(3);
    fvector->valueFE(fe, v);
    EXPECT_FLOAT_EQ(v(1), it[3]);
    EXPECT_FLOAT_EQ(v(2), it[4]);
    EXPECT_FLOAT_EQ(v(3), it[5]);
    EXPECT_FLOAT_EQ(fscalar->valueFE(fe), it[4]);
  }
}


TEST(TestLRSplineFields, Grad3D)
{
  ASMuCube patch;
  EXPECT_TRUE(patch.generateFEMTopology());

  // {x+y+z+x*y*z, x+y-z+x*y*z, x-y+z+x*y*z}
  std::vector<double> vc = {0.0,  0.0,  0.0,
                            1.0,  1.0,  1.0,
                            1.0,  1.0, -1.0,
                            2.0,  2.0,  0.0,
                            1.0, -1.0,  1.0,
                            2.0,  0.0,  2.0,
                            2.0,  0.0,  0.0,
                            4.0,  2.0,  2.0};
  std::unique_ptr<Fields> fvector(Fields::create(&patch, vc));
  static std::vector<std::pair<std::array<double,3>,
                               std::array<double,9>>> tests_vector =
    {{{{0.5, 0.5, 0.5}}, {{1.25,  1.25,  1.25,
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
                           2.0,  0.0,  2.0}}}};
  for (const auto& it : tests_vector) {
    FiniteElement fe;
    fe.u = it.first[0];
    fe.v = it.first[1];
    fe.w = it.first[2];
    Matrix gradu(3,3);
    fvector->gradFE(fe, gradu);
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j <3; ++j)
        EXPECT_FLOAT_EQ(gradu(i+1,j+1), it.second[i*3+j]);
  }
}


TEST(TestLRSplineFields, Value3Dmx)
{
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  ASMmxuCube patch({1,1,1,1});
  patch.raiseOrder(1,1,1,false);
  ASSERT_TRUE(patch.generateFEMTopology());

  // {x+y+z+x*y*z, x+y-z+x*y*z, x-y+z+x*y*z}
  std::vector<double> vc = {0.0, 0.5, 1.0,
                            1.0, 1.5, 2.0,
                            1.0, 1.5, 2.0,
                            2.0, 3.0, 4.0,
//                             y
                            0.0, 1.0,
                            0.5, 1.5,
                            1.0, 2.0,
                           -1.0, 0.0,
                           -0.5, 1.0,
                            0.0, 2.0,
//                             z
                            0.0, 1.0,
                           -1.0, 0.0,
                            0.5, 1.5,
                           -0.5, 1.0,
                            1.0, 2.0,
                            0.0, 2.0,
//                             p
                            0.0, 1.0, 1.0, 2.0, 1.0, 2.0, 2.0, 3.0};
  std::unique_ptr<Fields> fvector(Fields::create(&patch, vc, 123));
  std::unique_ptr<Field> fscalar(Field::create(&patch, vc, 4, 1));
  static std::vector<std::array<double,7>> tests_vector =
                         {{{0.5, 0.5, 0.5, 1.625, 0.625, 0.625, 1.5}},
                          {{0.0, 0.0, 0.0,   0.0,   0.0,   0.0, 0.0}},
                          {{1.0, 0.0, 0.0,   1.0,   1.0,   1.0, 1.0}},
                          {{0.0, 1.0, 0.0,   1.0,   1.0,  -1.0, 1.0}},
                          {{1.0, 1.0, 0.0,   2.0,   2.0,   0.0, 2.0}},
                          {{0.0, 0.0, 1.0,   1.0,  -1.0,   1.0, 1.0}},
                          {{1.0, 0.0, 1.0,   2.0,   0.0,   2.0, 2.0}},
                          {{0.0, 1.0, 1.0,   2.0,   0.0,   0.0, 2.0}},
                          {{1.0, 1.0, 1.0,   4.0,   2.0,   2.0, 3.0}}};
  for (const auto& it : tests_vector) {
    FiniteElement fe;
    fe.u = it[0];
    fe.v = it[1];
    fe.w = it[2];
    Vector v(3);
    fvector->valueFE(fe, v);
    EXPECT_FLOAT_EQ(v(1), it[3]);
    EXPECT_FLOAT_EQ(v(2), it[4]);
    EXPECT_FLOAT_EQ(v(3), it[5]);
    EXPECT_FLOAT_EQ(fscalar->valueFE(fe), it[6]);
  }
}
