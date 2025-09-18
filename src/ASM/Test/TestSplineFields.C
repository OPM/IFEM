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

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <memory>


using Catch::Matchers::WithinRel;


TEST_CASE("TestSplineFields.Value2D")
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
    REQUIRE(fvector->valueFE(fe,v));
    REQUIRE_THAT(v(1), WithinRel(it[2]));
    REQUIRE_THAT(v(2), WithinRel(it[3]));
    REQUIRE_THAT(fscalar->valueFE(fe), WithinRel(it[3]));
  }
}


TEST_CASE("TestSplineFields.Value2Dmx")
{
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  ASMmxSquare patch({1,1,1});
  patch.raiseOrder(1,1);
  REQUIRE(patch.generateFEMTopology());

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
    REQUIRE(fvector->valueFE(fe,v));
    REQUIRE_THAT(v(1), WithinRel(it[2]));
    REQUIRE_THAT(v(2), WithinRel(it[3]));
    REQUIRE_THAT(fscalar->valueFE(fe), WithinRel(it[4]));
  }
}


TEST_CASE("TestSplineFields.Grad2D")
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
    REQUIRE(fvector->gradFE(ItgPoint(it[0],it[1]),gradu));
    REQUIRE_THAT(gradu(1,1), WithinRel(it[2]));
    REQUIRE_THAT(gradu(1,2), WithinRel(it[3]));
    REQUIRE_THAT(gradu(2,1), WithinRel(it[4]));
    REQUIRE_THAT(gradu(2,2), WithinRel(it[5]));
  }
}


TEST_CASE("TestSplineFields.Value3D")
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
    REQUIRE(fvector->valueFE(fe,v));
    REQUIRE_THAT(v(1), WithinRel(it[3]));
    REQUIRE_THAT(v(2), WithinRel(it[4]));
    REQUIRE_THAT(v(3), WithinRel(it[5]));
    REQUIRE_THAT(fscalar->valueFE(fe), WithinRel(it[4]));
  }
}


TEST_CASE("TestSplineFields.Grad3D")
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
    REQUIRE(fvector->gradFE(ItgPoint(it.first.data()),gradu));
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j <3; ++j)
        REQUIRE_THAT(gradu(i+1,j+1), WithinRel(it.second[i*3+j]));
  }
}


TEST_CASE("TestSplineFields.Value3Dmx")
{
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  ASMmxCube patch({1,1,1,1});
  patch.raiseOrder(1,1,1);
  REQUIRE(patch.generateFEMTopology());

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
    REQUIRE(fvector->valueFE(fe,v));
    REQUIRE_THAT(v(1), WithinRel(it[3]));
    REQUIRE_THAT(v(2), WithinRel(it[4]));
    REQUIRE_THAT(v(3), WithinRel(it[5]));
    REQUIRE_THAT(fscalar->valueFE(fe), WithinRel(it[6]));
  }
}
