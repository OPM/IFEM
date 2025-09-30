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

#include "Catch2Support.h"

#include <array>
#include <memory>
#include <vector>


TEST_CASE("TestLRSplineFields.Value2D")
{
  ASMuSquare patch;
  REQUIRE(patch.generateFEMTopology());

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
    REQUIRE_THAT(v(1), WithinRel(it[2]));
    REQUIRE_THAT(v(2), WithinRel(it[3]));
    REQUIRE_THAT(fscalar->valueFE(fe), WithinRel(it[3]));
  }
}


TEST_CASE("TestLRSplineFields.Value2Dmx")
{
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  ASMmxuSquare patch({1,1,1});
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
    REQUIRE_THAT(v(1), WithinRel(it[2]));
    REQUIRE_THAT(v(2), WithinRel(it[3]));
    REQUIRE_THAT(fscalar->valueFE(fe), WithinRel(it[4]));
  }
}


TEST_CASE("TestLRSplineFields.Grad2D")
{
  ASMuSquare patch;
  REQUIRE(patch.generateFEMTopology());

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
    REQUIRE_THAT(gradu(1,1), WithinRel(it[2]));
    REQUIRE_THAT(gradu(1,2), WithinRel(it[3]));
    REQUIRE_THAT(gradu(2,1), WithinRel(it[4]));
    REQUIRE_THAT(gradu(2,2), WithinRel(it[5]));
  }
}


TEST_CASE("TestLRSplineFields.Value3D")
{
  ASMuCube patch;
  REQUIRE(patch.generateFEMTopology());

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
    REQUIRE_THAT(v(1), WithinRel(it[3]));
    REQUIRE_THAT(v(2), WithinRel(it[4]));
    REQUIRE_THAT(v(3), WithinRel(it[5]));
    REQUIRE_THAT(fscalar->valueFE(fe), WithinRel(it[4]));
  }
}


TEST_CASE("TestLRSplineFields.Grad3D")
{
  ASMuCube patch;
  REQUIRE(patch.generateFEMTopology());

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
        REQUIRE_THAT(gradu(i+1,j+1), WithinRel(it.second[i*3+j]));
  }
}


TEST_CASE("TestLRSplineFields.Value3Dmx")
{
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  ASMmxuCube patch({1,1,1,1});
  patch.raiseOrder(1,1,1,false);
  REQUIRE(patch.generateFEMTopology());

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
    REQUIRE_THAT(v(1), WithinRel(it[3]));
    REQUIRE_THAT(v(2), WithinRel(it[4]));
    REQUIRE_THAT(v(3), WithinRel(it[5]));
    REQUIRE_THAT(fscalar->valueFE(fe), WithinRel(it[6]));
  }
}
