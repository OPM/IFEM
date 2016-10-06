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
#include "FiniteElement.h"
#include "SIM2D.h"
#include "SIM3D.h"

#include "gtest/gtest.h"


TEST(TestSplineField, Value2D)
{
  SIM2D sim(1);
  sim.createDefaultModel();

  std::vector<double> sc = {0.0, 1.0, 1.0, 2.0}; // x + y
  Field* fscalar = Field::create(sim.getPatch(1), sc);
  static std::vector<std::array<double,3>> tests_scalar = {{0.5, 0.5, 1.0},
                                                           {1.0, 0.0, 1.0},
                                                           {0.0, 1.0, 1.0},
                                                           {1.0, 1.0, 2.0}};
  for (const auto& it : tests_scalar) {
    FiniteElement fe;
    fe.u = it[0];
    fe.v = it[1];
    ASSERT_FLOAT_EQ(fscalar->valueFE(fe), it[2]);
  }
}


TEST(TestSplineField, Grad2D)
{
  SIM2D sim(1);
  sim.createDefaultModel();

  std::vector<double> sc = {0.0, 1.0, 1.0, 2.0}; // x + y
  Field* fscalar = Field::create(sim.getPatch(1), sc);
  static std::vector<std::array<double,2>> tests_scalar = {{0.5, 0.5},
                                                           {1.0, 0.0},
                                                           {0.0, 1.0},
                                                           {1.0, 1.0}};
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


TEST(TestSplineField, Value3D)
{
  SIM3D sim(1);
  sim.createDefaultModel();

  std::vector<double> sc = {0.0, 1.0, 1.0, 2.0, 1.0, 2.0, 2.0, 3.0}; // x + y + z
  Field* fscalar = Field::create(sim.getPatch(1), sc);
  static std::vector<std::array<double,4>> tests_scalar = {{0.5, 0.5, 0.5, 1.5},
                                                           {0.0, 0.0, 0.0, 0.0},
                                                           {1.0, 0.0, 0.0, 1.0},
                                                           {0.0, 1.0, 0.0, 1.0},
                                                           {1.0, 1.0, 0.0, 2.0},
                                                           {0.0, 0.0, 1.0, 1.0},
                                                           {1.0, 0.0, 1.0, 2.0},
                                                           {0.0, 1.0, 1.0, 2.0},
                                                           {1.0, 1.0, 1.0, 3.0}};
  for (const auto& it : tests_scalar) {
    FiniteElement fe;
    fe.u = it[0];
    fe.v = it[1];
    fe.w = it[2];
    ASSERT_FLOAT_EQ(fscalar->valueFE(fe), it[3]);
  }
}


TEST(TestSplineField, Grad3D)
{
  SIM3D sim(1);
  sim.createDefaultModel();

  std::vector<double> sc = {0.0, 1.0, 1.0, 2.0, 1.0, 2.0, 2.0, 3.0}; // x + y + z
  Field* fscalar = Field::create(sim.getPatch(1), sc);
  static std::vector<std::array<double,3>> tests_scalar = {{0.5, 0.5, 0.5},
                                                           {0.0, 0.0, 0.0},
                                                           {1.0, 0.0, 0.0},
                                                           {0.0, 1.0, 0.0},
                                                           {1.0, 1.0, 0.0},
                                                           {0.0, 0.0, 1.0},
                                                           {1.0, 0.0, 1.0},
                                                           {0.0, 1.0, 1.0},
                                                           {1.0, 1.0, 1.0}};
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
