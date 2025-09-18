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

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <memory>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;


TEST_CASE("TestLRSplineField.Value2D")
{
  ASMuSquare patch(1);
  REQUIRE(patch.generateFEMTopology());

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
    REQUIRE_THAT(fscalar->valueFE(fe), WithinRel(it[2]));
  }
}


TEST_CASE("TestLRSplineField.Grad2D")
{
  ASMuSquare patch(1);
  REQUIRE(patch.generateFEMTopology());

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
    REQUIRE_THAT(v(1), WithinRel(1.0));
    REQUIRE_THAT(v(2), WithinRel(1.0));
  }
}


TEST_CASE("TestLRSplineField.Value3D")
{
  ASMuCube patch(1);
  REQUIRE(patch.generateFEMTopology());

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
    REQUIRE_THAT(fscalar->valueFE(fe), WithinRel(it[3]));
  }
}


TEST_CASE("TestLRSplineField.Grad3D")
{
  ASMuCube patch(1);
  REQUIRE(patch.generateFEMTopology());

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
    REQUIRE_THAT(v(1), WithinRel(1.0));
    REQUIRE_THAT(v(2), WithinRel(1.0));
    REQUIRE_THAT(v(3), WithinRel(1.0));
  }
}
