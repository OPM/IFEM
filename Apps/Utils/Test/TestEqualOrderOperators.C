//==============================================================================
//!
//! \file TestEqualOrderOperators.C
//!
//! \date Feb 16 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for various discrete equal-order operators.
//!
//==============================================================================

#include "EqualOrderOperators.h"
#include "FiniteElement.h"

#include "Catch2Support.h"


namespace
{
  void check_matrix_equal (const Matrix& A, const Real2DMat& B)
  {
    for (size_t i = 1; i <= A.rows(); ++i)
      for (size_t j = 1; j <= A.cols(); ++j)
        REQUIRE_THAT(A(i,j), WithinRel(B[i-1][j-1]));
  }

  class MyElement : public FiniteElement
  {
  public:
    MyElement() : FiniteElement(2)
    {
      dNdX.resize(2,2);
      dNdX(1,1) = 1.0;
      dNdX(2,1) = 2.0;
      dNdX(1,2) = 3.0;
      dNdX(2,2) = 4.0;
      N(1) = 1.0;
      N(2) = 2.0;
    }
  };
}


TEST_CASE("TestEqualOrderOperators.Advection")
{
  MyElement fe;

  Matrix EM_scalar(2, 2);
  Vec3   U(1.0, 2.0);

  // First component only
  EqualOrderOperators::Weak::Advection(EM_scalar, fe, U);
  const Real2DMat EM_scalar_ref = {{ 7.0, 10},
                                   {14.0, 20.0}};
  check_matrix_equal(EM_scalar, EM_scalar_ref);

  Matrix EM_vec(2*2, 2*2);
  Vec3   V(3.0, 4.0);

  EqualOrderOperators::Weak::Advection(EM_vec, fe, V, 2.0);
  const Real2DMat EM_vec_ref = {{30.0,  0.0, 44.0,  0.0},
                                { 0.0, 30.0,  0.0, 44.0},
                                {60.0,  0.0, 88.0,  0.0},
                                { 0.0, 60.0,  0.0, 88.0}};
  check_matrix_equal(EM_vec, EM_vec_ref);
}


TEST_CASE("TestEqualOrderOperators.Divergence")
{
  MyElement fe;

  Matrix EM_scalar(2, 2*2);
  EqualOrderOperators::Weak::Divergence(EM_scalar, fe);
  const Real2DMat EM_scalar_ref = {{1.0, 3.0, 2.0, 4.0},
                                   {2.0, 6.0, 4.0, 8.0}};
  check_matrix_equal(EM_scalar, EM_scalar_ref);

  Vector EV_scalar(4);
  Vec3   D(1.0, 2.0);

  EqualOrderOperators::Weak::Divergence(EV_scalar, fe, D, 1.0);
  REQUIRE_THAT(EV_scalar(1), WithinRel(7.0));
  REQUIRE_THAT(EV_scalar(2), WithinRel(10.0));
  REQUIRE_THAT(EV_scalar(3), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(EV_scalar(4), WithinAbs(0.0, 1e-13));
}


TEST_CASE("TestEqualOrderOperators.Gradient")
{
  MyElement fe;

  Matrix EM_scalar(2*2,2);
  // single component + pressure
  EqualOrderOperators::Weak::Gradient(EM_scalar, fe);
  const Real2DMat EM_scalar_ref = {{-1.0,-2.0},
                                   {-3.0,-6.0},
                                   {-2.0,-4.0},
                                   {-4.0,-8.0}};
  check_matrix_equal(EM_scalar, EM_scalar_ref);

  Vector EV_scalar(4);
  EqualOrderOperators::Residual::Gradient(EV_scalar, fe);
  REQUIRE_THAT(EV_scalar(1), WithinRel(fe.dNdX(1,1)));
  REQUIRE_THAT(EV_scalar(2), WithinRel(fe.dNdX(1,2)));
  REQUIRE_THAT(EV_scalar(3), WithinRel(fe.dNdX(2,1)));
  REQUIRE_THAT(EV_scalar(4), WithinRel(fe.dNdX(2,2)));
}


TEST_CASE("TestEqualOrderOperators.Laplacian")
{
  MyElement fe;

  // single scalar block
  Matrix EM_scalar(2,2);
  EqualOrderOperators::Weak::Laplacian(EM_scalar, fe);

  const Real2DMat EM_scalar_ref = {{10.0, 14.0},
                                   {14.0, 20.0}};

  check_matrix_equal(EM_scalar, EM_scalar_ref);

  // multiple (2) blocks in 3 component element matrix
  Matrix EM_multi(2*2,2*2);
  EqualOrderOperators::Weak::Laplacian(EM_multi, fe, 1.0);

  const Real2DMat EM_multi_ref = {{10.0,  0.0, 14.0,  0.0},
                                  { 0.0, 10.0,  0.0, 14.0},
                                  {14.0,  0.0, 20.0,  0.0},
                                   {0.0, 14.0,  0.0, 20.0}};

  check_matrix_equal(EM_multi, EM_multi_ref);

  // stress formulation
  Matrix EM_stress(2*2,2*2);
  EqualOrderOperators::Weak::Laplacian(EM_stress, fe, 1.0, true);
  const Real2DMat EM_stress_ref = {{11.0,  3.0, 16.0,  6.0},
                                   { 3.0, 19.0,  4.0, 26.0},
                                   {16.0,  4.0, 24.0,  8.0},
                                   { 6.0, 26.0,  8.0, 36.0}};

  check_matrix_equal(EM_stress, EM_stress_ref);

  Matrix EM_coeff(2, 2);
  EqualOrderOperators::Weak::LaplacianCoeff(EM_coeff, fe.dNdX, fe, 2.0);
  const Real2DMat EM_coeff_ref = {{104.0, 148.0},
                                  {152.0, 216.0}};
  check_matrix_equal(EM_coeff, EM_coeff_ref);

  Matrix EM_coefv(2, 2);
  EqualOrderOperators::Weak::LaplacianCoeff(EM_coefv, Vec3(1.0, 2.0), fe);
  const Real2DMat EM_coeff_refv = {{19.0, 26.0},
                                   {26.0, 36.0}};
  check_matrix_equal(EM_coefv, EM_coeff_refv);
}


TEST_CASE("TestEqualOrderOperators.Mass")
{
  MyElement fe;

  Matrix EM_scalar(2,2);
  EqualOrderOperators::Weak::Mass(EM_scalar, fe, 2.0);
  const Real2DMat EM_scalar_ref = {{2.0, 4.0},
                                   {4.0, 8.0}};
  check_matrix_equal(EM_scalar, EM_scalar_ref);

  Matrix EM_vec(2*2,2*2);
  EqualOrderOperators::Weak::Mass(EM_vec, fe);

  const Real2DMat EM_vec_ref = {{1.0, 0.0, 2.0, 0.0},
                                {0.0, 1.0, 0.0, 2.0},
                                {2.0, 0.0, 4.0, 0.0},
                                {0.0, 2.0, 0.0, 4.0}};
  check_matrix_equal(EM_vec, EM_vec_ref);
}


TEST_CASE("TestEqualOrderOperators.Source")
{
  MyElement fe;

  Vector EV_scalar(2);
  EqualOrderOperators::Weak::Source(EV_scalar, fe, 2.0);
  REQUIRE_THAT(EV_scalar(1), WithinRel(2.0));
  REQUIRE_THAT(EV_scalar(2), WithinRel(4.0));

  Vector EV_vec(4);
  Vec3 f(1.0, 2.0);

  EqualOrderOperators::Weak::Source(EV_vec, fe, f, 1.0);
  REQUIRE_THAT(EV_vec(1), WithinRel(1.0));
  REQUIRE_THAT(EV_vec(2), WithinRel(2.0));
  REQUIRE_THAT(EV_vec(3), WithinRel(2.0));
  REQUIRE_THAT(EV_vec(4), WithinRel(4.0));
  EV_vec.fill(0.0);
  EqualOrderOperators::Weak::Source(EV_vec, fe, 2.0, 0);
  REQUIRE_THAT(EV_vec(1), WithinRel(2.0));
  REQUIRE_THAT(EV_vec(2), WithinRel(2.0));
  REQUIRE_THAT(EV_vec(3), WithinRel(4.0));
  REQUIRE_THAT(EV_vec(4), WithinRel(4.0));
  EV_vec.fill(0.0);
  EqualOrderOperators::Weak::Source(EV_vec, fe, 2.0, 1);
  REQUIRE_THAT(EV_vec(1), WithinRel(2.0));
  REQUIRE_THAT(EV_vec(2), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(EV_vec(3), WithinRel(4.0));
  REQUIRE_THAT(EV_vec(4), WithinAbs(0.0, 1e-13));
  EV_vec.fill(0.0);
  EqualOrderOperators::Weak::Source(EV_vec, fe, 2.0, 2);
  REQUIRE_THAT(EV_vec(1), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(EV_vec(2), WithinRel(2.0));
  REQUIRE_THAT(EV_vec(3), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(EV_vec(4), WithinRel(4.0));
}
