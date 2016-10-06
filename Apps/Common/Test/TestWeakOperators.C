//==============================================================================
//!
//! \file TestWeakoperators.C
//!
//! \date Feb 16 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for various mesh quality indicators
//!
//==============================================================================

#include "WeakOperators.h"
#include "gtest/gtest.h"
#include "FiniteElement.h"

typedef std::vector<std::vector<double>> DoubleVec;
const auto&& check_matrix_equal = [](const Matrix& A, const DoubleVec& B)
                                  {
                                    for (size_t i=1;i<=A.rows();++i)
                                      for (size_t j=1;j<=A.cols();++j)
                                        ASSERT_NEAR(A(i,j), B[i-1][j-1], 1e-13);
                                  };

static FiniteElement getFE()
{
  FiniteElement fe(2);
  fe.dNdX.resize(2,2);
  fe.dNdX(1,1) = 1.0;
  fe.dNdX(2,1) = 2.0;
  fe.dNdX(1,2) = 3.0;
  fe.dNdX(2,2) = 4.0;
  fe.N(1) = 1.0;
  fe.N(2) = 2.0;

  return fe;
}

TEST(TestWeakOperators, Advection)
{
  FiniteElement fe = getFE();

  Matrix EM_scalar(2, 2);
  Vec3 U;
  U[0] = 1.0; U[1] = 2.0;
  // First component only
  WeakOperators::Advection(EM_scalar, fe, U);
  const DoubleVec EM_scalar_ref = {{ 7.0, 10},
                                   {14.0, 20.0}};
  check_matrix_equal(EM_scalar, EM_scalar_ref);

  Matrix EM_vec(2*2, 2*2);
  U[0] = 3.0; U[1] = 4.0;
  WeakOperators::Advection(EM_vec, fe, U, 2.0);
  const DoubleVec EM_vec_ref = {{30.0,  0.0, 44.0,  0.0},
                                { 0.0, 30.0,  0.0, 44.0},
                                {60.0,  0.0, 88.0,  0.0},
                                { 0.0, 60.0,  0.0, 88.0}};
  check_matrix_equal(EM_vec, EM_vec_ref);
}


TEST(TestWeakOperators, Divergence)
{
  FiniteElement fe = getFE();

  Matrix EM_scalar(2, 2*2);
  WeakOperators::Divergence(EM_scalar, fe);
  const DoubleVec EM_scalar_ref = {{1.0, 3.0, 2.0, 4.0},
                                   {2.0, 6.0, 4.0, 8.0}};
  check_matrix_equal(EM_scalar, EM_scalar_ref);

  Vector EV_scalar(4);
  Vec3 D;
  D[0] = 1.0; D[1] = 2.0;
  WeakOperators::Divergence(EV_scalar, fe, D, 1.0);
  ASSERT_NEAR(EV_scalar(1),  7.0, 1e-13);
  ASSERT_NEAR(EV_scalar(2), 10.0, 1e-13);
  ASSERT_NEAR(EV_scalar(3),  0.0, 1e-13);
  ASSERT_NEAR(EV_scalar(4),  0.0, 1e-13);
}


TEST(TestWeakOperators, Gradient)
{
  FiniteElement fe = getFE();

  Matrix EM_scalar(2*2,2);
  // single component + pressure
  WeakOperators::Gradient(EM_scalar, fe);
  const DoubleVec EM_scalar_ref = {{-1.0,-2.0},
                                   {-3.0,-6.0},
                                   {-2.0,-4.0},
                                   {-4.0,-8.0}};
  check_matrix_equal(EM_scalar, EM_scalar_ref);

  Vector EV_scalar(4);
  WeakOperators::Gradient(EV_scalar, fe);
  ASSERT_NEAR(EV_scalar(1),  fe.dNdX(1,1), 1e-13);
  ASSERT_NEAR(EV_scalar(2),  fe.dNdX(1,2), 1e-13);
  ASSERT_NEAR(EV_scalar(3),  fe.dNdX(2,1), 1e-13);
  ASSERT_NEAR(EV_scalar(4),  fe.dNdX(2,2), 1e-13);
}


TEST(TestWeakOperators, Laplacian)
{
  FiniteElement fe = getFE();

  // single scalar block
  Matrix EM_scalar(2,2);
  WeakOperators::Laplacian(EM_scalar, fe);

  const DoubleVec EM_scalar_ref = {{10.0, 14.0},
                                   {14.0, 20.0}};

  check_matrix_equal(EM_scalar, EM_scalar_ref);

  // multiple (2) blocks in 3 component element matrix
  Matrix EM_multi(2*2,2*2);
  WeakOperators::Laplacian(EM_multi, fe, 1.0);

  const DoubleVec EM_multi_ref = {{10.0,  0.0, 14.0,  0.0},
                                  { 0.0, 10.0,  0.0, 14.0},
                                  {14.0,  0.0, 20.0,  0.0},
                                   {0.0, 14.0,  0.0, 20.0}};

  check_matrix_equal(EM_multi, EM_multi_ref);

  // stress formulation
  Matrix EM_stress(2*2,2*2);
  WeakOperators::Laplacian(EM_stress, fe, 1.0, true);
  const DoubleVec EM_stress_ref = {{11.0,  3.0, 16.0,  6.0},
                                   { 3.0, 19.0,  4.0, 26.0},
                                   {16.0,  4.0, 24.0,  8.0},
                                   { 6.0, 26.0,  8.0, 36.0}};

  check_matrix_equal(EM_stress, EM_stress_ref);

  Matrix EM_coeff(2, 2);
  WeakOperators::LaplacianCoeff(EM_coeff, fe.dNdX, fe, 2.0);
  const DoubleVec EM_coeff_ref = {{104.0, 148.0},
                                  {152.0, 216.0}};
  check_matrix_equal(EM_coeff, EM_coeff_ref);
}


TEST(TestWeakOperators, Mass)
{
  FiniteElement fe = getFE();

  Matrix EM_scalar(2,2);
  WeakOperators::Mass(EM_scalar, fe, 2.0);
  const DoubleVec EM_scalar_ref = {{2.0, 4.0},
                                   {4.0, 8.0}};
  check_matrix_equal(EM_scalar, EM_scalar_ref);

  Matrix EM_vec(2*2,2*2);
  WeakOperators::Mass(EM_vec, fe);

  const DoubleVec EM_vec_ref = {{1.0, 0.0, 2.0, 0.0},
                                {0.0, 1.0, 0.0, 2.0},
                                {2.0, 0.0, 4.0, 0.0},
                                {0.0, 2.0, 0.0, 4.0}};
  check_matrix_equal(EM_vec, EM_vec_ref);
}


TEST(TestWeakOperators, Source)
{
  FiniteElement fe = getFE();

  Vector EV_scalar(2);
  WeakOperators::Source(EV_scalar, fe, 2.0);
  ASSERT_NEAR(EV_scalar(1),  2.0, 1e-13);
  ASSERT_NEAR(EV_scalar(2),  4.0, 1e-13);

  Vector EV_vec(4);
  Vec3 f;
  f[0] = 1.0; f[1] = 2.0;
  WeakOperators::Source(EV_vec, fe, f, 1.0);
  ASSERT_NEAR(EV_vec(1),  1.0, 1e-13);
  ASSERT_NEAR(EV_vec(2),  2.0, 1e-13);
  ASSERT_NEAR(EV_vec(3),  2.0, 1e-13);
  ASSERT_NEAR(EV_vec(4),  4.0, 1e-13);
  EV_vec.fill(0.0);
  WeakOperators::Source(EV_vec, fe, 2.0, 0);
  ASSERT_NEAR(EV_vec(1),  2.0, 1e-13);
  ASSERT_NEAR(EV_vec(2),  2.0, 1e-13);
  ASSERT_NEAR(EV_vec(3),  4.0, 1e-13);
  ASSERT_NEAR(EV_vec(4),  4.0, 1e-13);
  EV_vec.fill(0.0);
  WeakOperators::Source(EV_vec, fe, 2.0, 1);
  ASSERT_NEAR(EV_vec(1),  2.0, 1e-13);
  ASSERT_NEAR(EV_vec(2),  0.0, 1e-13);
  ASSERT_NEAR(EV_vec(3),  4.0, 1e-13);
  ASSERT_NEAR(EV_vec(4),  0.0, 1e-13);
  EV_vec.fill(0.0);
  WeakOperators::Source(EV_vec, fe, 2.0, 2);
  ASSERT_NEAR(EV_vec(1),  0.0, 1e-13);
  ASSERT_NEAR(EV_vec(2),  2.0, 1e-13);
  ASSERT_NEAR(EV_vec(3),  0.0, 1e-13);
  ASSERT_NEAR(EV_vec(4),  4.0, 1e-13);
}
