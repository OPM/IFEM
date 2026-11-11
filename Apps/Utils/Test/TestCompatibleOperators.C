//==============================================================================
//!
//! \file TestCompatibleOperators.C
//!
//! \date Oct 17 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for various discrete div-compatible operators.
//!
//==============================================================================

#include "CompatibleOperators.h"
#include "FiniteElement.h"

#include "Catch2Support.h"

#include <array>
#include <cmath>
#include <vector>


using DoubleVec = std::vector<std::vector<double>>;
const auto check_matrix_equal = [](const Matrix& A, const DoubleVec& B)
                                {
                                  for (size_t i=1;i<=A.rows();++i)
                                    for (size_t j=1;j<=A.cols();++j)
                                      REQUIRE_THAT(A(i,j), WithinRel(B[i-1][j-1]));
                                };

class MyFiniteElement : public MxFiniteElement
{
public:
  MyFiniteElement() : MxFiniteElement({6,6,4})
  {
    for (size_t i = 1; i <= 2; ++i) {
      this->grad(i).resize(6,2);
      this->basis(i).resize(6);
      for (size_t j = 1; j <= 6; ++j) {
          this->grad(i)(j,1) = 1.0+12*(i-1)+2*(j-1);
          this->grad(i)(j,2) = 2.0+12*(i-1)+2*(j-1);
          this->basis(i)(j) = j + 6*(i-1);
      }
    }

    this->grad(3).resize(4,2);
    this->basis(3).resize(4);
    for (size_t j = 1; j <= 4; ++j) {
      this->grad(3)(j,1) = 1.0 + 14.0 + 2*(j-1);
      this->grad(3)(j,2) = 2.0 + 14.0 + 2*(j-1);
      this->basis(3)(j) = 12 + j;
    }
  }
};


/*
TESTI(TestEqualOrderOperators, Advection)
{
  FiniteElement fe = getFE();

  Matrix EM_scalar(2, 2);
  Vec3 U;
  U[0] = 1.0; U[1] = 2.0;
  // First component only
  EqualOrderOperators::Weak::Advection(EM_scalar, fe, U);
  const DoubleVec EM_scalar_ref = {{ 7.0, 10},
                                   {14.0, 20.0}};
  check_matrix_equal(EM_scalar, EM_scalar_ref);

  Matrix EM_vec(2*2, 2*2);
  U[0] = 3.0; U[1] = 4.0;
  EqualOrderOperators::Weak::Advection(EM_vec, fe, U, 2.0);
  const DoubleVec EM_vec_ref = {{30.0,  0.0, 44.0,  0.0},
                                { 0.0, 30.0,  0.0, 44.0},
                                {60.0,  0.0, 88.0,  0.0},
                                { 0.0, 60.0,  0.0, 88.0}};
  check_matrix_equal(EM_vec, EM_vec_ref);
}


TESTI(TestEqualOrderOperators, Divergence)
{
  FiniteElement fe = getFE();

  Matrix EM_scalar(2, 2*2);
  EqualOrderOperators::Weak::Divergence(EM_scalar, fe);
  const DoubleVec EM_scalar_ref = {{1.0, 3.0, 2.0, 4.0},
                                   {2.0, 6.0, 4.0, 8.0}};
  check_matrix_equal(EM_scalar, EM_scalar_ref);

  Vector EV_scalar(4);
  Vec3 D;
  D[0] = 1.0; D[1] = 2.0;
  EqualOrderOperators::Weak::Divergence(EV_scalar, fe, D, 1.0);
  ASSERT_NEAR(EV_scalar(1),  7.0, 1e-13);
  ASSERT_NEAR(EV_scalar(2), 10.0, 1e-13);
  ASSERT_NEAR(EV_scalar(3),  0.0, 1e-13);
  ASSERT_NEAR(EV_scalar(4),  0.0, 1e-13);
}


TESTI(TestEqualOrderOperators, Gradient)
{
  FiniteElement fe = getFE();

  Matrix EM_scalar(2*2,2);
  // single component + pressure
  EqualOrderOperators::Weak::Gradient(EM_scalar, fe);
  const DoubleVec EM_scalar_ref = {{-1.0,-2.0},
                                   {-3.0,-6.0},
                                   {-2.0,-4.0},
                                   {-4.0,-8.0}};
  check_matrix_equal(EM_scalar, EM_scalar_ref);

  Vector EV_scalar(4);
  EqualOrderOperators::Weak::Gradient(EV_scalar, fe);
  ASSERT_NEAR(EV_scalar(1),  fe.dNdX(1,1), 1e-13);
  ASSERT_NEAR(EV_scalar(2),  fe.dNdX(1,2), 1e-13);
  ASSERT_NEAR(EV_scalar(3),  fe.dNdX(2,1), 1e-13);
  ASSERT_NEAR(EV_scalar(4),  fe.dNdX(2,2), 1e-13);
}
*/


TEST_CASE("TestCompatibleOperators.Laplacian")
{
  MyFiniteElement fe;

  std::vector<Matrix> EM(10);
  EM[1].resize(6,6);
  EM[2].resize(6,6);

  static constexpr auto idx = std::array{
      std::array{1, 4, 5},
      std::array{6, 2, 7},
      std::array{8, 9, 3},
  };

  CompatibleOperators::Weak::Laplacian(EM, fe, idx);

  const DoubleVec EM_1_ref = {{ 5.0, 11.0,  17.0,  23.0,  29.0,  35.0},
                              {11.0, 25.0,  39.0,  53.0,  67.0,  81.0},
                              {17.0, 39.0,  61.0,  83.0, 105.0, 127.0},
                              {23.0, 53.0,  83.0, 113.0, 143.0, 173.0},
                              {29.0, 67.0, 105.0, 143.0, 181.0, 219.0},
                              {35.0, 81.0, 127.0, 173.0, 219.0, 265.0}};

  check_matrix_equal(EM[1], EM_1_ref);

  const DoubleVec EM_2_ref = {{365.0, 419.0, 473.0, 527.0,  581.0,  635.0},
                              {419.0, 481.0, 543.0, 605.0,  667.0,  729.0},
                              {473.0, 543.0, 613.0, 683.0,  753.0,  823.0},
                              {527.0, 605.0, 683.0, 761.0,  839.0,  917.0},
                              {581.0, 667.0, 753.0, 839.0,  925.0, 1011.0},
                              {635.0, 729.0, 823.0, 917.0, 1011.0, 1105.0}};

  check_matrix_equal(EM[2], EM_2_ref);

  // stress formulation
  std::vector<Matrix> EM_stress(10);

  EM_stress[1].resize(6,6);
  EM_stress[2].resize(6,6);
  EM_stress[4].resize(6,6);
  EM_stress[6].resize(6,6);

  CompatibleOperators::Weak::Laplacian(EM_stress, fe, idx, 1.0, true);
  std::cout << EM_stress[1] << std::endl;

/*
  const DoubleVec EM_stress_ref = {{11.0,  3.0, 16.0,  6.0},
                                   { 3.0, 19.0,  4.0, 26.0},
                                   {16.0,  4.0, 24.0,  8.0},
                                   { 6.0, 26.0,  8.0, 36.0}};

  check_matrix_equal(EM_stress, EM_stress_ref);

  Matrix EM_coeff(2, 2);
  EqualOrderOperators::Weak::LaplacianCoeff(EM_coeff, fe.dNdX, fe, 2.0);
  const DoubleVec EM_coeff_ref = {{104.0, 148.0},
                                  {152.0, 216.0}};
  check_matrix_equal(EM_coeff, EM_coeff_ref);
*/
}


TEST_CASE("TestCompatibleOperators.Mass")
{
  MyFiniteElement fe;

  std::vector<Matrix> EM_vec(10);
  EM_vec[1].resize(6,6);
  EM_vec[2].resize(6,6);

  static constexpr auto idx = std::array{
      std::array{1, 4, 5},
      std::array{6, 2, 7},
      std::array{8, 9, 3},
  };

  CompatibleOperators::Weak::Mass(EM_vec, fe, idx);

  const DoubleVec EM_1_ref = {{1.0,  2.0,  3.0,  4.0,  5.0,  6.0},
                              {2.0,  4.0,  6.0,  8.0, 10.0, 12.0},
                              {3.0,  6.0,  9.0, 12.0, 15.0, 18.0},
                              {4.0,  8.0, 12.0, 16.0, 20.0, 24.0},
                              {5.0, 10.0, 15.0, 20.0, 25.0, 30.0},
                              {6.0, 12.0, 18.0, 24.0, 30.0, 36.0}};

  check_matrix_equal(EM_vec[1], EM_1_ref);

  const DoubleVec EM_2_ref = {{49.0, 56.0,  63.0,  70.0,  77.0,  84.0},
                              {56.0, 64.0,  72.0,  80.0,  88.0,  96.0},
                              {63.0, 72.0,  81.0,  90.0,  99.0, 108.0},
                              {70.0, 80.0,  90.0, 100.0, 110.0, 120.0},
                              {77.0, 88.0,  99.0, 110.0, 121.0, 132.0},
                              {84.0, 96.0, 108.0, 120.0, 132.0, 144.0}};

  check_matrix_equal(EM_vec[2], EM_2_ref);
}


TEST_CASE("TestCompatibleOperators.Source")
{
  MyFiniteElement fe;

  Vectors EV_scalar(4);
  EV_scalar[1].resize(6);
  EV_scalar[2].resize(6);

  CompatibleOperators::Weak::Source(EV_scalar, fe, {1,2,3}, 2.0);

  for (size_t i = 1; i <= 6; ++i) {
    REQUIRE_THAT(EV_scalar[1](i), WithinRel(2.0*fe.basis(1)(i)));
    REQUIRE_THAT(EV_scalar[2](i), WithinRel(2.0*fe.basis(2)(i)));
  }

  EV_scalar[1].fill(0.0);
  EV_scalar[2].fill(0.0);
  Vec3 f;
  f[0] = 1.0;
  f[1] = 2.0;

  CompatibleOperators::Weak::Source(EV_scalar, fe, f, {1,2,3}, 1.0);
  for (size_t i = 1; i <= 6; ++i) {
    REQUIRE_THAT(EV_scalar[1](i), WithinRel(fe.basis(1)(i)));
    REQUIRE_THAT(EV_scalar[2](i), WithinRel(2.0*fe.basis(2)(i)));
  }

  EV_scalar[1].fill(0.0);
  EV_scalar[2].fill(0.0);
  CompatibleOperators::Weak::Source(EV_scalar, fe, f, {1,2,3}, 2.0);
  for (size_t i = 1; i <= 6; ++i) {
    REQUIRE_THAT(EV_scalar[1](i), WithinRel(2.0*fe.basis(1)(i)));
    REQUIRE_THAT(EV_scalar[2](i), WithinRel(4.0*fe.basis(2)(i)));
  }
}


TEST_CASE("TestCompatibleOperators.LaplacianCoeff")
{
  MyFiniteElement fe;

  static constexpr auto idx = std::array{
      std::array{1, 4, 5},
      std::array{6, 2, 7},
      std::array{8, 9, 3},
  };

  // The identity tensor contracts the two gradients as they are, so the
  // term is the plain laplacian and is checked against it rather than
  // against numbers of its own
  const size_t nsd = fe.dNdX.cols();
  Matrix K(nsd,nsd);
  for (size_t i = 1; i <= nsd; ++i)
    K(i,i) = 1.0;

  std::vector<Matrix> EM(10), EM_ref(10);
  for (size_t b = 1; b <= 2; ++b) {
    EM[b].resize(6,6);
    EM_ref[b].resize(6,6);
  }

  CompatibleOperators::Weak::LaplacianCoeff(EM, K, fe, idx, 2.0);
  CompatibleOperators::Weak::Laplacian(EM_ref, fe, idx, 2.0);

  for (size_t b = 1; b <= 2; ++b)
    for (size_t i = 1; i <= 6; ++i)
      for (size_t j = 1; j <= 6; ++j)
        REQUIRE_THAT(EM[b](i,j), WithinAbs(EM_ref[b](i,j), 1e-13));

  // A tensor which is not a multiple of the identity has to differ from it
  K(1,2) = K(2,1) = 0.5;
  std::vector<Matrix> EM_skew(10);
  for (size_t b = 1; b <= 2; ++b)
    EM_skew[b].resize(6,6);

  CompatibleOperators::Weak::LaplacianCoeff(EM_skew, K, fe, idx, 2.0);

  double diff = 0.0;
  for (size_t b = 1; b <= 2; ++b)
    for (size_t i = 1; i <= 6; ++i)
      for (size_t j = 1; j <= 6; ++j)
        diff += std::abs(EM_skew[b](i,j) - EM_ref[b](i,j));

  REQUIRE(diff > 1.0);

  // The blocks the contributions land in are the ones idx names, so a
  // remapped diagonal has to move them and not merely relabel them. Back to
  // the identity, which is what the plain laplacian compares against.
  K(1,2) = K(2,1) = 0.0;
  static constexpr auto swapped = std::array{
      std::array{3, 4, 5},
      std::array{6, 1, 7},
      std::array{8, 9, 2},
  };

  std::vector<Matrix> EM_swap(10), EM_swap_ref(10);
  for (size_t b = 1; b <= 3; ++b) {
    EM_swap[b].resize(6,6);
    EM_swap_ref[b].resize(6,6);
  }

  CompatibleOperators::Weak::LaplacianCoeff(EM_swap, K, fe, swapped, 2.0);
  CompatibleOperators::Weak::Laplacian(EM_swap_ref, fe, swapped, 2.0);

  for (size_t b = 1; b <= 3; ++b)
    for (size_t i = 1; i <= 6; ++i)
      for (size_t j = 1; j <= 6; ++j)
        REQUIRE_THAT(EM_swap[b](i,j), WithinAbs(EM_swap_ref[b](i,j), 1e-13));
}


namespace {

constexpr size_t nsdFD = 2; //!< Dimensions the tangent is differenced over
constexpr size_t nenFD = 6; //!< Functions each velocity component carries

//! \brief Evaluates the advected field and its gradient at a nodal state.
void stateFD (const MxFiniteElement& fe, const std::vector<double>& u,
              Vec3& U, Tensor& dUdX)
{
  for (size_t k = 1; k <= nsdFD; ++k) {
    U[k-1] = 0.0;
    for (size_t i = 1; i <= nenFD; ++i)
      U[k-1] += fe.basis(k)(i) * u[(k-1)*nenFD + i-1];
    for (size_t l = 1; l <= nsdFD; ++l) {
      dUdX(k,l) = 0.0;
      for (size_t i = 1; i <= nenFD; ++i)
        dUdX(k,l) += fe.grad(k)(i,l) * u[(k-1)*nenFD + i-1];
    }
  }
}


//! \brief Returns the element residual of the convection at a nodal state.
Vectors residualFD (const MxFiniteElement& fe, const std::vector<double>& u,
                    WeakOperators::ConvectionForm form, double scale)
{
  Vec3 U;
  Tensor dUdX(nsdFD);
  stateFD(fe,u,U,dUdX);

  Vectors EV(3);
  EV[1].resize(nenFD);
  EV[2].resize(nenFD);
  static constexpr auto idx = std::array{1, 2, 0};
  CompatibleOperators::Residual::Convection(EV,fe,U,dUdX,U,idx,scale,form);

  return EV;
}

}


/*!
  The tangent of the convection is what the Newton iteration converges by, and
  nothing in the answer it converges to depends on it being right: a wrong one
  costs iterations, or the iteration itself, and says nothing while it does.
  So it is differenced against the residual it is supposed to be the tangent
  of, in every nodal degree of freedom, for each form of the term.
*/

TEST_CASE("TestCompatibleOperators.ConvectionTangent")
{
  MyFiniteElement fe;
  fe.detJxW = 1.0;

  const double scale = 0.75;

  // A state with nothing special about it, so that no term of the tangent
  // can come out right by cancelling against another
  std::vector<double> u(nsdFD*nenFD);
  for (size_t i = 0; i < u.size(); ++i)
    u[i] = 0.1 + 0.03*i - 0.002*i*i;

  static constexpr auto idx = std::array{
      std::array{1, 4, 0},
      std::array{6, 2, 0},
      std::array{0, 0, 0},
  };

  const std::array forms = {WeakOperators::CONVECTIVE,
                            WeakOperators::CONSERVATIVE,
                            WeakOperators::SKEWSYMMETRIC};

  for (const auto form : forms) {
    Vec3 U;
    Tensor dUdX(nsdFD);
    stateFD(fe,u,U,dUdX);

    std::vector<Matrix> EM(10);
    for (const int b : {1, 2, 4, 6})
      EM[b].resize(nenFD,nenFD);

    CompatibleOperators::Weak::Convection(EM,fe,U,dUdX,idx,scale,form);

    // The residual is the negative of what the tangent is taken of, so the
    // difference of the one is the negative of a column of the other. It is
    // quadratic in the state, so a central difference of it is exact at any
    // step, and the step is taken large to keep the round-off of the
    // subtraction, which grows as the step shrinks, well under the check.
    const double eps = 1.0e-3;
    for (size_t l = 1; l <= nsdFD; ++l)
      for (size_t j = 1; j <= nenFD; ++j) {
        std::vector<double> up(u), um(u);
        up[(l-1)*nenFD + j-1] += eps;
        um[(l-1)*nenFD + j-1] -= eps;
        const Vectors Rp = residualFD(fe,up,form,scale);
        const Vectors Rm = residualFD(fe,um,form,scale);
        for (size_t k = 1; k <= nsdFD; ++k)
          for (size_t i = 1; i <= nenFD; ++i) {
            const double fd = (Rp[k](i) - Rm[k](i)) / (2.0*eps);
            REQUIRE_THAT(-EM[idx[k-1][l-1]](i,j),
                         WithinAbs(fd, 1.0e-8*std::max(1.0,std::fabs(fd))));
          }
      }
  }
}
