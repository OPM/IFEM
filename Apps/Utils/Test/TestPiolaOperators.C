//==============================================================================
//!
//! \file TestPiolaOperators.C
//!
//! \date Sep 17 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for various discrete Piola-mapped operators.
//!
//==============================================================================

#include "PiolaOperators.h"
#include "FiniteElement.h"
#include "Tensor.h"
#include "Vec3.h"

#include "Catch2Support.h"


/*!
  A finite element with a Piola mapped basis whose values and derivatives are
  arbitrary but fixed. Nothing here has to be a real mapping of a real basis:
  the tests below check that the residual operators are the weak operators
  applied to a coefficient vector, which is an algebraic identity that has to
  hold whatever the two matrices contain. The numbers are deliberately without
  structure, since a real mapping is only as general as its geometry, and an
  affine one gives a dPdX which is rank one in its component and derivative
  indices, hiding a swap of those two.
*/

class PiolaFiniteElement : public MxFiniteElement
{
public:
  PiolaFiniteElement() : MxFiniteElement({3,3,2})
  {
    detJxW = 0.75;

    // The geometry basis gradient is only read for its number of columns,
    // which is the number of space dimensions
    this->grad(1).resize(3,2);
    this->grad(2).resize(3,2);
    this->grad(3).resize(2,2);
    for (size_t b = 1; b <= 3; ++b)
      for (size_t i = 1; i <= this->basis(b).size(); ++i)
        this->basis(b)(i) = 0.5 + 0.25*i + 0.125*b;

    P.resize(2,6);
    dPdX.resize(4,6);
    for (size_t j = 1; j <= 6; ++j) {
      for (size_t k = 1; k <= 2; ++k)
        P(k,j) = 0.3 + 0.7*j - 1.1*k + 0.05*j*k;
      for (size_t r = 1; r <= 4; ++r)
        dPdX(r,j) = 1.0 - 0.3*j + 0.9*r - 0.11*j*r;
    }
  }

  //! \brief The coefficients of an arbitrary field in this basis.
  static Vector coefficients()
  {
    Vector c(6);
    for (size_t i = 1; i <= 6; ++i)
      c(i) = 0.4*i - 0.15*i*i + 1.0;
    return c;
  }

  //! \brief The gradient of that field, as velocityGradient forms it.
  Tensor gradient() const
  {
    Matrix dV;
    dV.multiplyMat(dPdX, coefficients());
    Tensor grad(2);
    grad = dV;
    return grad;
  }

  //! \brief The value of that field.
  Vec3 value() const
  {
    Vector u;
    P.multiply(coefficients(), u);
    return Vec3(u(1), u(2));
  }
};


namespace {

//! \brief Block indices, one per velocity basis, laid out as PiolaOperators
//! expects for the diagonal and the vector blocks.
const std::array<std::array<int,3>,3> matIdx {{{0,2,0},{3,1,0},{0,0,0}}};
const std::array<int,3> vecIdx {0,1,0};

//! \brief Applies the assembled matrices to a coefficient vector.
Vector applyMats(const Matrices& EM, const Vector& c)
{
  Vector out(6);
  // Basis 1 holds coefficients 1-3 and basis 2 coefficients 4-6, so the
  // diagonal blocks take their own three and the off-diagonal ones the others
  const std::array<std::array<int,2>,2> blk {{{0,2},{3,1}}};
  for (size_t bi = 0; bi < 2; ++bi)
    for (size_t bj = 0; bj < 2; ++bj) {
      const Matrix& A = EM[blk[bi][bj]];
      for (size_t i = 1; i <= 3; ++i)
        for (size_t j = 1; j <= 3; ++j)
          out(3*bi+i) += A(i,j) * c(3*bj+j);
    }
  return out;
}

//! \brief Joins the two velocity residual blocks into one vector.
Vector join(const Vectors& EV)
{
  Vector out(6);
  for (size_t b = 0; b < 2; ++b)
    for (size_t i = 1; i <= 3; ++i)
      out(3*b+i) = EV[b](i);
  return out;
}

Matrices emptyMats()
{
  Matrices EM(4);
  for (Matrix& A : EM) A.resize(3,3);
  return EM;
}

Vectors emptyVecs()
{
  Vectors EV(2);
  for (Vector& v : EV) v.resize(3);
  return EV;
}

//! \brief The convective residual for a given set of coefficients.
Vector convectionResidual(const PiolaFiniteElement& fe, const Vector& c,
                          WeakOperators::ConvectionForm form)
{
  Vector u;
  fe.P.multiply(c, u);
  const Vec3 U(u(1), u(2));

  Matrix dV;
  dV.multiplyMat(fe.dPdX, c);
  Tensor dUdX(2);
  dUdX = dV;

  Vectors EV = emptyVecs();
  PiolaOperators::Residual::Convection(EV, fe, U, dUdX, U, vecIdx, 1.0, form);
  return join(EV);
}

}


namespace {

//! \brief The three forms of the convective term.
const auto convectionForms = []()
{
  return GENERATE(WeakOperators::CONVECTIVE,
                  WeakOperators::CONSERVATIVE,
                  WeakOperators::SKEWSYMMETRIC);
};

}


TEST_CASE("TestPiolaOperators.AdvectionIsTheTangentOfConvection")
{
  // The residual of a convective term is the advection matrix, built with the
  // same advecting field, applied to the coefficients of the advected one.
  // This holds for each of the three forms, the conservative one being the
  // negated transpose of the convective one and the skew symmetric one their
  // mean, on both sides of the identity.
  const WeakOperators::ConvectionForm form = convectionForms();

  const PiolaFiniteElement fe;
  const Vector c = PiolaFiniteElement::coefficients();
  const Vec3 U = fe.value();
  const Tensor dUdX = fe.gradient();

  Matrices EM = emptyMats();
  PiolaOperators::Weak::Advection(EM, fe, U, matIdx, 1.0, form);

  Vectors EV = emptyVecs();
  PiolaOperators::Residual::Convection(EV, fe, U, dUdX, U, vecIdx, 1.0, form);

  const Vector Au = applyMats(EM, c);
  const Vector r = join(EV);
  for (size_t i = 1; i <= 6; ++i)
    REQUIRE_THAT(r(i), WithinAbs(-Au(i), 1.0e-10));
}


TEST_CASE("TestPiolaOperators.LaplacianResidualMatchesItsMatrix")
{
  // Both forms of the viscous term: the plain gradient product, and the stress
  // formulation which adds the transposed one to make up 2*eps(u):eps(v)
  const bool stress = GENERATE(false, true);

  const PiolaFiniteElement fe;
  const Vector c = PiolaFiniteElement::coefficients();
  const Tensor dUdX = fe.gradient();

  Matrices EM = emptyMats();
  PiolaOperators::Weak::Laplacian(EM, fe, matIdx, 1.0, stress);

  Vectors EV = emptyVecs();
  PiolaOperators::Residual::Laplacian(EV, fe, dUdX, vecIdx, 1.0, stress);

  const Vector Au = applyMats(EM, c);
  const Vector r = join(EV);
  for (size_t i = 1; i <= 6; ++i)
    REQUIRE_THAT(r(i), WithinAbs(-Au(i), 1.0e-10));
}


TEST_CASE("TestPiolaOperators.ConvectionIsTheTangentOfItsResidual")
{
  // The weak convection operator is the Newton tangent of the convective
  // residual, so differencing the residual has to reproduce it. The residual
  // is quadratic in the coefficients, which makes a central difference exact
  // up to round-off, hence the tight tolerance.
  const WeakOperators::ConvectionForm form = convectionForms();

  const PiolaFiniteElement fe;
  const Vector c = PiolaFiniteElement::coefficients();
  const Vec3 U = fe.value();
  const Tensor dUdX = fe.gradient();

  Matrices EM = emptyMats();
  PiolaOperators::Weak::Convection(EM, fe, U, dUdX, matIdx, 1.0, form);

  const double h = 1.0e-3;
  for (size_t j = 1; j <= 6; ++j) {
    Vector cp(c), cm(c);
    cp(j) += h;
    cm(j) -= h;
    const Vector rp = convectionResidual(fe, cp, form);
    const Vector rm = convectionResidual(fe, cm, form);

    // Column j of the tangent, which the residual carries with a minus sign
    Vector col(6);
    for (size_t i = 1; i <= 6; ++i)
      col(i) = -(rp(i) - rm(i)) / (2.0*h);

    Vector e(6);
    e(j) = 1.0;
    const Vector Ae = applyMats(EM, e);
    for (size_t i = 1; i <= 6; ++i)
      REQUIRE_THAT(Ae(i), WithinAbs(col(i), 1.0e-7));
  }
}
