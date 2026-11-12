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

#include <array>
#include <cmath>
#include <vector>


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


TEST_CASE("TestEqualOrderOperators.GradDiv")
{
  MyElement fe;

  // The term pairs the divergence of the two functions, so the matrix is
  // the outer product of the one vector holding every derivative of every
  // component: the entry of test (i,k) against trial (j,l) is the product
  // of dN_i/dx_k with dN_j/dx_l and nothing else.
  Matrix EM(2*2,2*2);
  EqualOrderOperators::Weak::GradDiv(EM, fe);

  const Real2DMat EM_ref = {{1.0,  3.0, 2.0,  4.0},
                            {3.0,  9.0, 6.0, 12.0},
                            {2.0,  6.0, 4.0,  8.0},
                            {4.0, 12.0, 8.0, 16.0}};

  check_matrix_equal(EM, EM_ref);

  // Which is what the stress operator is not: the two agree wherever one
  // of the indices is shared and part company where none is
  Matrix EM_stress(2*2,2*2);
  EqualOrderOperators::Weak::Stress(EM_stress, fe);
  REQUIRE_THAT(EM_stress(1,4), WithinRel(6.0));

  // A divergence free velocity is in the null space of the penalty, which
  // is the whole of what it is there to leave alone
  Vector u(2*2);
  for (size_t i = 1; i <= 2; ++i)
    for (size_t k = 1; k <= 2; ++k)
      u((i-1)*2+k) = k == 1 ? fe.dNdX(i,2) : -fe.dNdX(i,1);

  Vector Au;
  EM.multiply(u, Au);
  for (size_t i = 1; i <= Au.size(); ++i)
    REQUIRE_THAT(Au(i), WithinAbs(0.0, 1.0e-13));

  Matrix EM_scaled(2*2,2*2);
  EqualOrderOperators::Weak::GradDiv(EM_scaled, fe, 2.0);
  for (size_t i = 1; i <= 4; ++i)
    for (size_t j = 1; j <= 4; ++j)
      REQUIRE_THAT(EM_scaled(i,j), WithinRel(2.0*EM(i,j)));
}


namespace {

constexpr size_t nsdFD = 2; //!< Dimensions the tangent is differenced over
constexpr size_t nenFD = 4; //!< Functions the basis carries

//! \brief An element whose basis is of no particular shape.
class FDElement : public FiniteElement
{
public:
  //! \brief The constructor fills the basis and its derivatives.
  FDElement() : FiniteElement(nenFD)
  {
    dNdX.resize(nenFD,nsdFD);
    for (size_t i = 1; i <= nenFD; ++i) {
      N(i) = 0.3 + 0.2*i - 0.01*i*i;
      for (size_t l = 1; l <= nsdFD; ++l)
        dNdX(i,l) = 0.5 - 0.17*i + 0.31*l + 0.05*i*l;
    }
    detJxW = 1.0;
  }
};


//! \brief Evaluates the advected field and its gradient at a nodal state.
void stateFD (const FiniteElement& fe, const std::vector<double>& u,
              Vec3& U, Tensor& dUdX)
{
  for (size_t k = 1; k <= nsdFD; ++k) {
    U[k-1] = 0.0;
    for (size_t i = 1; i <= nenFD; ++i)
      U[k-1] += fe.N(i) * u[(i-1)*nsdFD + k-1];
    for (size_t l = 1; l <= nsdFD; ++l) {
      dUdX(k,l) = 0.0;
      for (size_t i = 1; i <= nenFD; ++i)
        dUdX(k,l) += fe.dNdX(i,l) * u[(i-1)*nsdFD + k-1];
    }
  }
}


//! \brief Returns the element residual of the convection at a nodal state.
Vector residualFD (const FiniteElement& fe, const std::vector<double>& u,
                   WeakOperators::ConvectionForm form, double scale)
{
  Vec3 U;
  Tensor dUdX(nsdFD);
  stateFD(fe,u,U,dUdX);

  Vector EV(nsdFD*nenFD);
  EqualOrderOperators::Residual::Convection(EV,fe,U,dUdX,U,scale,form);

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

TEST_CASE("TestEqualOrderOperators.ConvectionTangent")
{
  FDElement fe;

  const double scale = 0.75;

  // A state with nothing special about it, so that no term of the tangent
  // can come out right by cancelling against another
  std::vector<double> u(nsdFD*nenFD);
  for (size_t i = 0; i < u.size(); ++i)
    u[i] = 0.1 + 0.03*i - 0.002*i*i;

  const std::array forms = {WeakOperators::CONVECTIVE,
                            WeakOperators::CONSERVATIVE,
                            WeakOperators::SKEWSYMMETRIC};

  for (const auto form : forms) {
    Vec3 U;
    Tensor dUdX(nsdFD);
    stateFD(fe,u,U,dUdX);

    Matrix EM(nsdFD*nenFD, nsdFD*nenFD);
    EqualOrderOperators::Weak::Convection(EM,fe,U,dUdX,scale,form);

    // The residual is the negative of what the tangent is taken of, so the
    // difference of the one is the negative of a column of the other. It is
    // quadratic in the state, so a central difference of it is exact at any
    // step, and the step is taken large to keep the round-off of the
    // subtraction, which grows as the step shrinks, well under the check.
    const double eps = 1.0e-3;
    for (size_t j = 1; j <= nenFD; ++j)
      for (size_t l = 1; l <= nsdFD; ++l) {
        std::vector<double> up(u), um(u);
        up[(j-1)*nsdFD + l-1] += eps;
        um[(j-1)*nsdFD + l-1] -= eps;
        const Vector Rp = residualFD(fe,up,form,scale);
        const Vector Rm = residualFD(fe,um,form,scale);
        for (size_t i = 1; i <= nenFD; ++i)
          for (size_t k = 1; k <= nsdFD; ++k) {
            const double fd = (Rp((i-1)*nsdFD+k) - Rm((i-1)*nsdFD+k))/(2.0*eps);
            REQUIRE_THAT(-EM((i-1)*nsdFD+k,(j-1)*nsdFD+l),
                         WithinAbs(fd, 1.0e-8*std::max(1.0,std::fabs(fd))));
          }
      }
  }
}
