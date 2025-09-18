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

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <memory>

using Catch::Matchers::WithinRel;


TEST_CASE("TestSplineField.Value2D")
{
  ASMSquare patch(1);
  std::vector<double> sc = {0.0, 1.0, 1.0, 2.0}; // x + y
  std::unique_ptr<Field>fscalar(Field::create(&patch,sc));
  static std::vector<std::array<double,3>> tests_scalar = {{{{0.5, 0.5, 1.0}},
                                                            {{1.0, 0.0, 1.0}},
                                                            {{0.0, 1.0, 1.0}},
                                                            {{1.0, 1.0, 2.0}}}};
  for (const auto& it : tests_scalar)
    REQUIRE_THAT(fscalar->valueFE(ItgPoint(it[0],it[1])), WithinRel(it[2]));
}


TEST_CASE("TestSplineField.Grad2D")
{
  ASMSquare patch(1);
  std::vector<double> sc = {0.0, 1.0, 1.0, 2.0}; // x + y
  std::unique_ptr<Field> fscalar(Field::create(&patch,sc));
  static std::vector<std::array<double,2>> tests_scalar = {{{{0.5, 0.5}},
                                                            {{1.0, 0.0}},
                                                            {{0.0, 1.0}},
                                                            {{1.0, 1.0}}}};
  for (const auto& it : tests_scalar) {
    Vector v(2);
    REQUIRE(fscalar->gradFE(ItgPoint(it[0],it[1]),v));
    REQUIRE_THAT(v(1), WithinRel(1.0));
    REQUIRE_THAT(v(2), WithinRel(1.0));
  }
}


TEST_CASE("TestSplineField.Value3D")
{
  ASMCube patch(1);
  std::vector<double> sc = {0.0, 1.0, 1.0, 2.0, 1.0, 2.0, 2.0, 3.0}; // x + y + z
  std::unique_ptr<Field> fscalar(Field::create(&patch,sc));
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
    REQUIRE_THAT(fscalar->valueFE(ItgPoint(it[0],it[1],it[2])), WithinRel(it[3]));
}


TEST_CASE("TestSplineField.Grad3D")
{
  ASMCube patch(1);
  std::vector<double> sc = {0.0, 1.0, 1.0, 2.0, 1.0, 2.0, 2.0, 3.0}; // x + y + z
  std::unique_ptr<Field> fscalar(Field::create(&patch,sc));
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
    REQUIRE(fscalar->gradFE(ItgPoint(it[0],it[1],it[2]),v));
    REQUIRE_THAT(v(1), WithinRel(1.0));
    REQUIRE_THAT(v(2), WithinRel(1.0));
    REQUIRE_THAT(v(3), WithinRel(1.0));
  }
}
