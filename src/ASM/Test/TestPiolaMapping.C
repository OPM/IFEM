//==============================================================================
//!
//! \file TestPiolaMapping.C
//!
//! \date Aug 11 2023
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for utilities for Piola mapping transformations.
//!
//==============================================================================

#include "ASMs2D.h"
#include "ASMmxBase.h"
#include "ASMs2Dmx.h"
#include "ASMs3D.h"
#include "ASMs3Dmx.h"
#include "CoordinateMapping.h"
#include "FiniteElement.h"
#include "SplineUtils.h"
#include "Vec3.h"

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"

#include "Catch2Support.h"

#include <array>
#include <functional>
#include <sstream>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;


namespace utl
{
  std::vector<matrix<Real>> jacobianGradient(const matrix<Real>& dudX,
                                             const matrix3d<Real>& d2Xdu2)
  {
    std::vector<matrix<Real>> dJ;
    utl::JacobianGradient(dudX,d2Xdu2,dJ);
    return dJ;
  }

  utl::vector<Real> determinantGradient(const matrix<Real>& J,
                                        const matrix<Real>& Ji,
                                        const matrix3d<Real>& H)
  {
    utl::vector<Real> dJ;
    utl::detJacGradient(J,Ji,H,dJ);
    return dJ;
  }
}


static const char* spline2D = R"(200 1 0 0
2 0
3 3
0 0 0 1 1 1
3 3
0 0 0 1 1 1
0 0
2 -0.5
4 -1
-0.5 2.5
2.25 1.75
5 1
-1 5
2.5 4
6 3
)";

static const char* spline2D_basis = R"(200 1 0 0
2 0
4 4
0 0 0 0 1 1 1 1
4 4
0 0 0 0 1 1 1 1
0 0
0.3333333333333333 0
0.6666666666666666 0
1 0
0 0.3333333333333333
0.333333333333333 0.3333333333333331
0.6666666666666663 0.3333333333333332
0.9999999999999993 0.3333333333333333
0 0.6666666666666666
0.3333333333333332 0.6666666666666663
0.6666666666666669 0.6666666666666667
1 0.6666666666666666
0 1
0.3333333333333333 0.9999999999999993
0.6666666666666666 1
1 1
)";

static const char* spline3D = R"(700 1 0 0
3 0
4 4
0 0 0 0 1 1 1 1
4 4
0 0 0 0 1 1 1 1
4 4
0 0 0 0 1 1 1 1
0 0 0
5.551115123125783e-17 0 5.551115123125783e-17
0.3333333333333331 0 0.3333333333333331
1 0 1
0 0.6666666666666665 0
4.930380657631324e-32 0.6666666666666662 4.930380657631324e-32
0.333333333333333 0.6666666666666664 0.333333333333333
0.9999999999999993 0.6666666666666665 0.9999999999999993
0 1.333333333333333 0
5.551115123125783e-17 1.333333333333333 5.551115123125783e-17
0.3333333333333331 1.333333333333333 0.3333333333333331
1 1.333333333333333 1
0 2 0
5.551115123125783e-17 1.999999999999999 5.551115123125783e-17
0.3333333333333331 2 0.3333333333333331
1 2 1
0 1.110223024625157e-16 0.9999999999999998
4.930380657631324e-32 1.110223024625156e-16 1.000000000000001
0.333333333333333 1.110223024625157e-16 1.333333333333331
0.9999999999999993 1.110223024625157e-16 2
-1.54197642309049e-18 0.6666666666666669 0.9999999999999991
4.440892098500626e-16 0.6666666666666664 1
0.3333333333333321 0.6666666666666671 1.33333333333333
0.9999999999999997 0.6666666666666669 1.999999999999998
0.1111111111111111 1.333333333333331 0.9999999999999999
0.1111111111111104 1.33333333333333 1.000000000000001
0.4444444444444451 1.333333333333331 1.333333333333331
1.11111111111111 1.333333333333331 2
0.3333333333333333 1.999999999999999 0.9999999999999998
0.3333333333333338 1.999999999999998 1.000000000000001
0.666666666666665 1.999999999999999 1.333333333333331
1.333333333333333 1.999999999999999 2
0 -2.220446049250313e-16 2
5.551115123125783e-17 -2.220446049250312e-16 2.000000000000001
0.3333333333333331 -2.220446049250313e-16 2.333333333333332
1 -2.220446049250313e-16 3
3.39234813079909e-17 0.666666666666666 1.999999999999999
-2.590520390792031e-16 0.6666666666666655 2
0.3333333333333339 0.666666666666666 2.333333333333331
0.9999999999999986 0.666666666666666 2.999999999999999
0.2222222222222221 1.333333333333334 2
0.2222222222222228 1.333333333333333 2.000000000000001
0.5555555555555544 1.333333333333334 2.333333333333331
1.222222222222223 1.333333333333334 3.000000000000001
0.6666666666666666 2 2
0.6666666666666652 1.999999999999999 2.000000000000001
1.000000000000001 2 2.333333333333332
1.666666666666666 2 3
0 1 3
5.551115123125783e-17 0.9999999999999993 2.999999999999998
0.3333333333333331 1 3.333333333333334
1 1 4
5.551115123125783e-17 1.666666666666666 2.999999999999998
-3.33066907387547e-16 1.666666666666665 2.999999999999997
0.3333333333333334 1.666666666666665 3.333333333333331
0.9999999999999997 1.666666666666666 3.999999999999997
0.3333333333333331 2.333333333333332 3
0.3333333333333335 2.333333333333331 2.999999999999999
0.6666666666666664 2.333333333333332 3.333333333333334
1.333333333333333 2.333333333333332 4
1 3 3
0.9999999999999997 2.999999999999998 2.999999999999998
1.333333333333333 3 3.333333333333334
2 3 4
)";

static const char* spline3D_basis = R"(700 1 0 0
3 0
4 4
0 0 0 0 1 1 1 1
4 4
0 0 0 0 1 1 1 1
4 4
0 0 0 0 1 1 1 1
0 0 0
0.3333333333333333 0 0
0.6666666666666666 0 0
1 0 0
0 0.3333333333333333 0
0.333333333333333 0.3333333333333331 0
0.6666666666666663 0.3333333333333332 0
0.9999999999999993 0.3333333333333333 0
0 0.6666666666666666 0
0.3333333333333332 0.6666666666666663 0
0.6666666666666669 0.6666666666666667 0
1 0.6666666666666666 0
0 1 0
0.3333333333333333 0.9999999999999993 0
0.6666666666666666 1 0
1 1 0
0 0 0.3333333333333333
0.333333333333333 0 0.3333333333333331
0.6666666666666663 0 0.3333333333333332
0.9999999999999993 0 0.3333333333333333
0 0.333333333333333 0.3333333333333331
0.3333333333333329 0.3333333333333328 0.3333333333333329
0.6666666666666661 0.333333333333333 0.3333333333333331
0.9999999999999988 0.333333333333333 0.3333333333333331
0 0.6666666666666663 0.3333333333333332
0.333333333333333 0.666666666666666 0.333333333333333
0.6666666666666667 0.6666666666666664 0.3333333333333333
0.9999999999999993 0.6666666666666663 0.3333333333333332
0 0.9999999999999993 0.3333333333333333
0.333333333333333 0.9999999999999988 0.3333333333333331
0.6666666666666663 0.9999999999999993 0.3333333333333332
0.9999999999999993 0.9999999999999993 0.3333333333333333
0 0 0.6666666666666666
0.3333333333333332 0 0.6666666666666663
0.6666666666666669 0 0.6666666666666667
1 0 0.6666666666666666
0 0.3333333333333332 0.6666666666666663
0.333333333333333 0.333333333333333 0.666666666666666
0.6666666666666663 0.3333333333333333 0.6666666666666664
0.9999999999999993 0.3333333333333332 0.6666666666666663
0 0.6666666666666669 0.6666666666666667
0.3333333333333333 0.6666666666666664 0.6666666666666663
0.6666666666666666 0.6666666666666671 0.6666666666666666
1 0.6666666666666669 0.6666666666666667
0 1 0.6666666666666666
0.3333333333333332 0.9999999999999993 0.6666666666666663
0.6666666666666669 1 0.6666666666666667
1 1 0.6666666666666666
0 0 1
0.3333333333333333 0 0.9999999999999993
0.6666666666666666 0 1
1 0 1
0 0.3333333333333333 0.9999999999999993
0.333333333333333 0.3333333333333331 0.9999999999999988
0.6666666666666663 0.3333333333333332 0.9999999999999993
0.9999999999999993 0.3333333333333333 0.9999999999999993
0 0.6666666666666666 1
0.3333333333333332 0.6666666666666663 0.9999999999999993
0.6666666666666669 0.6666666666666667 1
1 0.6666666666666666 1
0 1 1
0.3333333333333333 0.9999999999999993 0.9999999999999993
0.6666666666666666 1 1
1 1 1
)";


TEST_CASE("TestPiolaMapping.Basis2D")
{
  ASMs2D p;
  std::stringstream g2(spline2D);

  REQUIRE(p.read(g2));
  REQUIRE(p.generateFEMTopology());

  ASMs2Dmx b(2, {1,1,1});
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  std::stringstream basis(spline2D_basis);
  REQUIRE(b.read(basis));
  REQUIRE(b.generateFEMTopology());

  constexpr double u = 0.1;
  constexpr double v = 0.1;

  Vector N;
  Matrix X, dNdu;
  Matrix J, Ji, dNdX;
  p.extractBasis(u, v, N, dNdu);
  p.getElementCoordinates(X, 1);
  double detJ = utl::Jacobian(Ji, dNdX, X, dNdu, true);
  J.multiply(X,dNdu); // J = X * dNdu

  Go::BasisPtsSf spline1, spline2;
  b.getBasis(1)->computeBasis(u,v,spline1);
  b.getBasis(2)->computeBasis(u,v,spline2);
  MxFiniteElement fe({spline1.basisValues.size(),
                      spline2.basisValues.size()});
  fe.basis(1) = spline1.basisValues;
  fe.basis(2) = spline2.basisValues;
  const auto q = std::array<std::function<double(double)>,3> {
      [](double x) { return (1.0 - x) * (1.0 - x); },
      [](double x) { return 2.0 * x * (1.0 - x); },
      [](double x) { return x * x; },
  };

  const auto c = std::array<std::function<double(double)>,4> {
    [](double x) { return (1.0 - x) * (1.0 - x) * (1.0 - x); },
    [](double x) { return 3.0 * x * (1.0 - x) * (1.0 - x); },
    [](double x) { return 3.0 * x * x * (1.0 - x); },
    [](double x) { return x * x * x; },
  };

  fe.piolaBasis(detJ, J);

  size_t k = 1;
  for (size_t j = 0; j < 3; ++j)
    for (size_t i = 0; i < 4; ++i, ++k) {
      REQUIRE_THAT(fe.P_contra.N(1,k), WithinRel((3.0*v + 4.0) * c[i](u) * q[j](v) / detJ));
      REQUIRE_THAT(fe.P_contra.N(2,k), WithinRel(   (-v - 1.0) * c[i](u) * q[j](v) / detJ));
    }

  for (size_t j = 0; j < 4; ++j)
    for (size_t i = 0; i < 3; ++i, ++k) {
      REQUIRE_THAT(fe.P_contra.N(1,k), WithinRel((3.0*v - 1.0) * q[i](u) * c[j](v) / detJ));
      REQUIRE_THAT(fe.P_contra.N(2,k),  WithinRel(  (-u + 5.0) * q[i](u) * c[j](v) / detJ));
    }
}


TEST_CASE("TestPiolaMapping.Gradient2D")
{
  ASMs2D p;
  std::stringstream g2(spline2D);

  REQUIRE(p.read(g2));
  REQUIRE(p.generateFEMTopology());

  ASMs2Dmx b(2, {1,1,1});
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  std::stringstream basis(spline2D_basis);
  REQUIRE(b.read(basis));
  REQUIRE(b.generateFEMTopology());

  constexpr double u = 0.1;
  constexpr double v = 0.1;

  BasisValues bfs(3);

  Go::BasisDerivsSf spline1, spline2;
  b.getBasis(1)->computeBasis(u,v,spline1);
  b.getBasis(2)->computeBasis(u,v,spline2);
  MxFiniteElement fe({spline1.basisValues.size(), spline2.basisValues.size()});
  SplineUtils::extractBasis(spline1, fe.basis(1), bfs[0].dNdu);
  SplineUtils::extractBasis(spline2, fe.basis(2), bfs[1].dNdu);

  Matrix X, J, Ji, dNdX;
  p.extractBasis(u, v, bfs.back().N, bfs.back().dNdu, bfs.back().d2Ndu2);
  p.getElementCoordinates(X, 1);
  double detJ = fe.detJxW = utl::Jacobian(Ji, dNdX, X, bfs.back().dNdu, true);
  J.multiply(X,bfs.back().dNdu); // J = X * dNdu
  Matrix3D H;
  H.multiply(X,bfs.back().d2Ndu2);

  const auto q = std::array<std::function<double(double)>,3> {
      [](double x) { return (1.0 - x) * (1.0 - x); },
      [](double x) { return 2.0 * x * (1.0 - x); },
      [](double x) { return x * x; },
  };

  const auto qd = std::array<std::function<double(double)>,3> {
      [](double x) { return -2.0 * (1.0 - x); },
      [](double x) { return 2.0 - 4.0 * x; },
      [](double x) { return 2.0 * x; },
  };

  const auto c = std::array<std::function<double(double)>,4> {
    [](double x) { return (1.0 - x) * (1.0 - x) * (1.0 - x); },
    [](double x) { return 3.0 * x * (1.0 - x) * (1.0 - x); },
    [](double x) { return 3.0 * x * x * (1.0 - x); },
    [](double x) { return x * x * x; },
  };

  const auto cd = std::array<std::function<double(double)>,4> {
    [](double x) { return -3.0 * (1.0 - x) * (1.0 - x); },
    [](double x) { return 9.0 * x * x - 12.0 * x + 3.0; },
    [](double x) { return 3.0 * (2.0 - 3.0 * x) * x; },
    [](double x) { return 3.0 * x * x; },
  };

  fe.piolaMapping(Ji, X, bfs);
  std::vector<Real> det = utl::determinantGradient(J, Ji, H);
  Matrices dJdX = utl::jacobianGradient(Ji, H);

  size_t k = 1;
  for (size_t j = 0; j < 3; ++j)
    for (size_t i = 0; i < 4; ++i, ++k)
      for (size_t d1 = 1; d1 <= 2; ++d1)
        for (size_t d2 = 1; d2 <= 2; ++d2)
              REQUIRE_THAT(fe.P_contra.dNdX((d1 + 2 * (d2 - 1)),k),
                           WithinRel((-det[d2-1] / (detJ * detJ) * J(d1,1) +
                           dJdX[d2-1](d1,1) / detJ) * c[i](u)*q[j](v) +
                           J(d1,1) * (cd[i](u) * q[j](v) * Ji(1,d2) +
                                      c[i](u) * qd[j](v) * Ji(2,d2)) / detJ));

  for (size_t j = 0; j < 4; ++j)
    for (size_t i = 0; i < 3; ++i, ++k)
      for (size_t d1 = 1; d1 <= 2; ++d1)
        for (size_t d2 = 1; d2 <= 2; ++d2)
              REQUIRE_THAT(fe.P_contra.dNdX(d1 + 2 * (d2 - 1),k),
                          WithinRel((-det[d2-1] / (detJ * detJ) * J(d1,2) +
                          dJdX[d2-1](d1,2) / detJ) * q[i](u) * c[j](v) +
                          J(d1,2) * (qd[i](u) * c[j](v) * Ji(1,d2) +
                                     q[i](u) * cd[j](v) * Ji(2,d2)) / detJ));
}


TEST_CASE("TestPiolaMapping.GradJ2D")
{
  ASMs2D p;
  std::stringstream g2(spline2D);

  REQUIRE(p.read(g2));
  REQUIRE(p.generateFEMTopology());

  Vec3 Xp;
  double u = 0.1;
  double v = 0.1;
  SplineUtils::point(Xp, u, v, p.getBasis());
  REQUIRE_THAT(Xp[0],  WithinRel(3.0*u*v + 4.0*u - v));
  REQUIRE_THAT(Xp[1], WithinRel(-1.0*u*v - u + 5.0 * v));

  Vector N;
  Matrix X, dNdu;
  Matrix J, Ji, dNdX;
  Matrix3D d2Ndu2, H;
  p.extractBasis(u, v, N, dNdu, d2Ndu2);
  p.getElementCoordinates(X, 1);
  double detJ = utl::Jacobian(Ji, dNdX, X, dNdu, true);
  H.multiply(X,d2Ndu2);
  J.multiply(X,dNdu); // J = X * dNdu

  REQUIRE_THAT(detJ, WithinRel(-3.0*u*v - 4.0*u + 3.0*v*v + 17*v + 19));

  REQUIRE_THAT(J(1,1), WithinRel(3.0*v + 4.0));
  REQUIRE_THAT(J(1,2), WithinRel(3.0*v - 1.0));
  REQUIRE_THAT(J(2,1), WithinRel(-v - 1.0));
  REQUIRE_THAT(J(2,2), WithinRel(-u + 5.0));

  REQUIRE_THAT(Ji(1,1), WithinRel(J(2,2) / detJ));
  REQUIRE_THAT(Ji(1,2), WithinRel(-J(1,2) / detJ));
  REQUIRE_THAT(Ji(2,1), WithinRel(-J(2,1) / detJ));
  REQUIRE_THAT(Ji(2,2), WithinRel(J(1,1) / detJ));

  REQUIRE_THAT(H(1,1,1), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(1,1,2), WithinRel(3.0));
  REQUIRE_THAT(H(1,2,1), WithinRel(3.0));
  REQUIRE_THAT(H(1,2,2), WithinAbs(0.0, 1e-13));

  REQUIRE_THAT(H(2,1,1), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(2,1,2), WithinRel(-1.0));
  REQUIRE_THAT(H(2,2,1), WithinRel(-1.0));
  REQUIRE_THAT(H(2,2,2), WithinAbs(0.0, 1e-13));

  Matrices dJdX = utl::jacobianGradient(Ji, H);

  REQUIRE_THAT(dJdX[0](1,1), WithinRel(H(1,1,1)*Ji(1,1) + H(1,1,2)*Ji(2,1)));
  REQUIRE_THAT(dJdX[0](1,2), WithinRel(H(1,1,2)*Ji(1,1) + H(1,2,2)*Ji(2,1)));
  REQUIRE_THAT(dJdX[0](2,1), WithinRel(H(2,1,1)*Ji(1,1) + H(2,1,2)*Ji(2,1)));
  REQUIRE_THAT(dJdX[0](2,2), WithinRel(H(2,1,2)*Ji(1,1) + H(2,2,2)*Ji(2,1)));

  REQUIRE_THAT(dJdX[1](1,1), WithinRel(H(1,1,1)*Ji(1,2) + H(1,1,2)*Ji(2,2)));
  REQUIRE_THAT(dJdX[1](1,2), WithinRel(H(1,1,2)*Ji(1,2) + H(1,2,2)*Ji(2,2)));
  REQUIRE_THAT(dJdX[1](2,1), WithinRel(H(2,1,1)*Ji(1,2) + H(2,1,2)*Ji(2,2)));
  REQUIRE_THAT(dJdX[1](2,2), WithinRel(H(2,1,2)*Ji(1,2) + H(2,2,2)*Ji(2,2)));
}


TEST_CASE("TestPiolaMapping.GradDetJ2D")
{
  ASMs2D p;
  std::stringstream g2(spline2D);

  REQUIRE(p.read(g2));
  REQUIRE(p.generateFEMTopology());

  double u = 0.1;
  double v = 0.1;

  Vector N;
  Matrix X, dNdu;
  Matrix J, Ji, dNdX;
  Matrix3D d2Ndu2, H;
  p.extractBasis(u, v, N, dNdu, d2Ndu2);
  p.getElementCoordinates(X, 1);
  utl::Jacobian(Ji, dNdX, X, dNdu, true);
  H.multiply(X,d2Ndu2);
  J.multiply(X,dNdu); // J = X * dNdu

  std::vector<Real> det = utl::determinantGradient(J, Ji, H);

  double Ju = -1.0;
  double Jv =  3.0 * (v - u) + 14.0;

  REQUIRE_THAT(det[0], WithinRel(Ju * Ji(1,1) + Jv * Ji(2,1)));
  REQUIRE_THAT(det[1], WithinRel(Ju * Ji(1,2) + Jv * Ji(2,2)));
}


TEST_CASE("TestPiolaMapping.Basis3D")
{
  ASMs3D p;
  std::stringstream g2(spline3D);

  REQUIRE(p.read(g2));
  REQUIRE(p.generateFEMTopology());

  ASMs3Dmx b({1,1,1});
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  std::stringstream basis(spline3D_basis);
  REQUIRE(b.read(basis));
  REQUIRE(b.generateFEMTopology());

  constexpr double u = 0.1;
  constexpr double v = 0.1;
  constexpr double w = 0.1;

  Vector N;
  Matrix X, dNdu;
  Matrix J, Ji, dNdX;
  p.extractBasis(u, v, w, N, dNdu);
  p.getElementCoordinates(X, 1);
  double detJ = utl::Jacobian(Ji, dNdX, X, dNdu, true);
  J.multiply(X,dNdu); // J = X * dNdu

  Go::BasisPts spline1, spline2, spline3;
  b.getBasis(1)->computeBasis(u,v,w,spline1);
  b.getBasis(2)->computeBasis(u,v,w,spline2);
  b.getBasis(3)->computeBasis(u,v,w,spline3);
  MxFiniteElement fe({spline1.basisValues.size(),
                      spline2.basisValues.size(),
                      spline3.basisValues.size()});
  fe.basis(1) = spline1.basisValues;
  fe.basis(2) = spline2.basisValues;
  fe.basis(3) = spline3.basisValues;
  const auto q = std::array<std::function<double(double)>,3> {
      [](double x) { return (1.0 - x) * (1.0 - x); },
      [](double x) { return 2.0 * x * (1.0 - x); },
      [](double x) { return x * x; },
  };

  const auto c = std::array<std::function<double(double)>,4> {
    [](double x) { return (1.0 - x) * (1.0 - x) * (1.0 - x); },
    [](double x) { return 3.0 * x * (1.0 - x) * (1.0 - x); },
    [](double x) { return 3.0 * x * x * (1.0 - x); },
    [](double x) { return x * x * x; },
  };

  fe.piolaBasis(detJ, J);

  size_t n = 1;
  for (size_t k = 0; k < 3; ++k)
    for (size_t j = 0; j < 3; ++j)
      for (size_t i = 0; i < 4; ++i, ++n) {
        REQUIRE_THAT(fe.P_contra.N(1,n), WithinRel(2.0 * u * c[i](u) * q[j](v) * q[k](w) / detJ));
        REQUIRE_THAT(fe.P_contra.N(2,n), WithinAbs(0.0 * c[i](u) * q[j](v) * q[k](w) / detJ, 1e-13));
        REQUIRE_THAT(fe.P_contra.N(3,n), WithinRel(2.0 * u * c[i](u) * q[j](v) * q[k](w) / detJ));
      }

  for (size_t k = 0; k < 3; ++k)
    for (size_t j = 0; j < 4; ++j)
      for (size_t i = 0; i < 3; ++i, ++n) {
        REQUIRE_THAT(fe.P_contra.N(1,n), WithinRel(2.0 * v * w * q[i](u) * c[j](v) * q[k](w) / detJ));
        REQUIRE_THAT(fe.P_contra.N(2,n), WithinRel(        2.0 * q[i](u) * c[j](v) * q[k](w) / detJ));
        REQUIRE_THAT(fe.P_contra.N(3,n), WithinAbs(        0.0 * q[i](u) * c[j](v) * q[k](w) / detJ, 1e-13));
      }

  for (size_t k = 0; k < 4; ++k)
    for (size_t j = 0; j < 3; ++j)
      for (size_t i = 0; i < 3; ++i, ++n) {
        REQUIRE_THAT(fe.P_contra.N(1,n), WithinRel(      v * v * q[i](u) * q[j](v) * c[k](w) / detJ));
        REQUIRE_THAT(fe.P_contra.N(2,n), WithinRel(3.0 * w * w * q[i](u) * q[j](v) * c[k](w) / detJ));
        REQUIRE_THAT(fe.P_contra.N(3,n), WithinRel(        3.0 * q[i](u) * q[j](v) * c[k](w) / detJ));
      }
}


TEST_CASE("TestPiolaMapping.Gradient3D")
{
  ASMs3D p;
  std::stringstream g2(spline3D);

  REQUIRE(p.read(g2));
  REQUIRE(p.generateFEMTopology());

  ASMs3Dmx b({1,1,1});
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  std::stringstream basis(spline3D_basis);
  REQUIRE(b.read(basis));
  REQUIRE(b.generateFEMTopology());

  constexpr double u = 0.1;
  constexpr double v = 0.1;
  constexpr double w = 0.1;

  BasisValues bfs(4);

  Go::BasisDerivs spline1, spline2, spline3;
  b.getBasis(1)->computeBasis(u,v,w,spline1);
  b.getBasis(2)->computeBasis(u,v,w,spline2);
  b.getBasis(3)->computeBasis(u,v,w,spline3);
  MxFiniteElement fe({spline1.basisValues.size(),
                      spline2.basisValues.size(),
                      spline3.basisValues.size()});
  SplineUtils::extractBasis(spline1, fe.basis(1), bfs[0].dNdu);
  SplineUtils::extractBasis(spline2, fe.basis(2), bfs[1].dNdu);
  SplineUtils::extractBasis(spline3, fe.basis(3), bfs[2].dNdu);

  Matrix X, J, Ji, dNdX;
  Matrix3D H;
  p.extractBasis(u, v, w, bfs.back().N, bfs.back().dNdu, bfs.back().d2Ndu2);
  p.getElementCoordinates(X, 1);
  double detJ = fe.detJxW = utl::Jacobian(Ji, dNdX, X, bfs.back().dNdu, true);
  H.multiply(X,bfs.back().d2Ndu2);
  J.multiply(X,bfs.back().dNdu); // J = X * dNdu

  const auto q = std::array<std::function<double(double)>,3> {
      [](double x) { return (1.0 - x) * (1.0 - x); },
      [](double x) { return 2.0 * x * (1.0 - x); },
      [](double x) { return x * x; },
  };

  const auto qd = std::array<std::function<double(double)>,3> {
      [](double x) { return -2.0 * (1.0 - x); },
      [](double x) { return 2.0 - 4.0 * x; },
      [](double x) { return 2.0 * x; },
  };

  const auto c = std::array<std::function<double(double)>,4> {
    [](double x) { return (1.0 - x) * (1.0 - x) * (1.0 - x); },
    [](double x) { return 3.0 * x * (1.0 - x) * (1.0 - x); },
    [](double x) { return 3.0 * x * x * (1.0 - x); },
    [](double x) { return x * x * x; },
  };

  const auto cd = std::array<std::function<double(double)>,4> {
    [](double x) { return -3.0 * (1.0 - x) * (1.0 - x); },
    [](double x) { return 9.0 * x * x - 12.0 * x + 3.0; },
    [](double x) { return 3.0 * (2.0 - 3.0 * x) * x; },
    [](double x) { return 3.0 * x * x; },
  };

  fe.piolaMapping(Ji, X, bfs);
  std::vector<Real> det = utl::determinantGradient(J, Ji, H);
  Matrices dJdX = utl::jacobianGradient(Ji, H);

  size_t n = 1;
  for (size_t k = 0; k < 3; ++k)
    for (size_t j = 0; j < 3; ++j)
      for (size_t i = 0; i < 4; ++i, ++n)
        for (size_t d1 = 1; d1 <= 3; ++d1)
          for (size_t d2 = 1; d2 <= 3; ++d2)
            REQUIRE_THAT(fe.P_contra.dNdX(d1 + 3 * (d2-1),n),
                         WithinAbs((-det[d2-1] / (detJ * detJ) * J(d1,1) +
                         dJdX[d2-1](d1,1) / detJ) * c[i](u) * q[j](v) * q[k](w) +
                         J(d1,1) * (cd[i](u) * q[j](v) * q[k](w) * Ji(1,d2) +
                                    c[i](u) * qd[j](v) * q[k](w) * Ji(2,d2) +
                                    c[i](u) * q[j](v) * qd[k](w) * Ji(3,d2)) / detJ, 1e-13));

  for (size_t k = 0; k < 3; ++k)
    for (size_t j = 0; j < 4; ++j)
      for (size_t i = 0; i < 3; ++i, ++n)
        for (size_t d1 = 1; d1 <= 3; ++d1)
          for (size_t d2 = 1; d2 <= 3; ++d2)
            REQUIRE_THAT(fe.P_contra.dNdX(d1 + 3 * (d2-1),n),
                         WithinAbs((-det[d2-1] / (detJ * detJ) * J(d1,2) +
                         dJdX[d2-1](d1,2) / detJ) * q[i](u) * c[j](v) * q[k](w) +
                         J(d1,2) * (qd[i](u) * c[j](v) * q[k](w) * Ji(1,d2) +
                                    q[i](u) * cd[j](v) * q[k](w) * Ji(2,d2) +
                                    q[i](u) * c[j](v) * qd[k](w) * Ji(3,d2)) / detJ, 1e-13));

  for (size_t k = 0; k < 4; ++k)
    for (size_t j = 0; j < 3; ++j)
      for (size_t i = 0; i < 3; ++i, ++n)
        for (size_t d1 = 1; d1 <= 3; ++d1)
          for (size_t d2 = 1; d2 <= 3; ++d2)
            REQUIRE_THAT(fe.P_contra.dNdX(d1 + 3 * (d2-1),n),
                         WithinAbs((-det[d2-1] / (detJ * detJ) * J(d1,3) +
                         dJdX[d2-1](d1,3) / detJ) * q[i](u) * q[j](v) * c[k](w) +
                         J(d1,3) * (qd[i](u) * q[j](v) * c[k](w) * Ji(1,d2) +
                                    q[i](u) * qd[j](v) * c[k](w) * Ji(2,d2) +
                                    q[i](u) * q[j](v) * cd[k](w) * Ji(3,d2)) / detJ, 1e-13));
}


TEST_CASE("TestPiolaMapping.GradJ3D")
{
  ASMs3D p;
  std::stringstream g2(spline3D);

  REQUIRE(p.read(g2));
  REQUIRE(p.generateFEMTopology());

  Vec3 Xp;
  double u = 0.1;
  double v = 0.1;
  double w = 0.1;
  SplineUtils::point(Xp, u, v, w, p.getBasis());
  REQUIRE_THAT(Xp[0], WithinRel(u * u + v * v * w));
  REQUIRE_THAT(Xp[1], WithinRel(2.0 * v + w * w * w));
  REQUIRE_THAT(Xp[2], WithinRel(3.0 * w + u * u));

  Vector N;
  Matrix X, dNdu;
  Matrix J, Ji, dNdX;
  Matrix3D d2Ndu2, H;
  p.extractBasis(u, v, w, N, dNdu, d2Ndu2);
  p.getElementCoordinates(X, 1);
  double detJ = utl::Jacobian(Ji, dNdX, X, dNdu, true);
  H.multiply(X,d2Ndu2);
  J.multiply(X,dNdu); // J = X * dNdu

  REQUIRE_THAT(J(1,1), WithinRel(2.0 * u));
  REQUIRE_THAT(J(1,2), WithinRel(2.0 * v * w));
  REQUIRE_THAT(J(1,3), WithinRel(v * v));
  REQUIRE_THAT(J(2,1), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(J(2,2), WithinRel(2.0));
  REQUIRE_THAT(J(2,3), WithinRel(3.0 * w * w));
  REQUIRE_THAT(J(3,1), WithinRel(2.0 * u));
  REQUIRE_THAT(J(3,2), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(J(3,3), WithinRel(3.0));

  auto det = [](const double a,
                const double b,
                const double c,
                const double d)
  {
      return a*d - b*c;
  };

  double detRef =   J(1,1) * det(J(2,2), J(2,3), J(3,2), J(3,3))
                  - J(1,2) * det(J(2,1), J(2,3), J(3,1), J(3,3))
                  + J(1,3) * det(J(2,1), J(2,2), J(3,1), J(3,2));

  REQUIRE_THAT(detRef, WithinRel(detJ));

  REQUIRE_THAT(Ji(1,1), WithinRel( det(J(2,2), J(1,3), J(3,2), J(3,3)) / detJ));
  REQUIRE_THAT(Ji(1,2), WithinRel(-det(J(1,2), J(1,3), J(3,2), J(3,3)) / detJ));
  REQUIRE_THAT(Ji(1,3), WithinRel( det(J(1,2), J(2,2), J(1,3), J(2,3)) / detJ));
  REQUIRE_THAT(Ji(2,1), WithinRel(-det(J(2,1), J(2,3), J(3,1), J(3,3)) / detJ));
  REQUIRE_THAT(Ji(2,2), WithinRel( det(J(1,1), J(1,3), J(3,1), J(3,3)) / detJ));
  REQUIRE_THAT(Ji(2,3), WithinRel(-det(J(1,1), J(1,3), J(2,1), J(2,3)) / detJ));
  REQUIRE_THAT(Ji(3,1), WithinRel( det(J(2,1), J(2,2), J(3,1), J(3,2)) / detJ));
  REQUIRE_THAT(Ji(3,2), WithinRel(-det(J(1,1), J(1,2), J(3,1), J(3,2)) / detJ));
  REQUIRE_THAT(Ji(3,3), WithinRel( det(J(1,1), J(1,2), J(2,1), J(2,2)) / detJ));

  REQUIRE_THAT(H(1,1,1), WithinRel(2.0));
  REQUIRE_THAT(H(1,1,2), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(1,1,3), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(1,2,1), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(1,2,2), WithinRel(2.0*w));
  REQUIRE_THAT(H(1,2,3), WithinRel(2.0*v));
  REQUIRE_THAT(H(1,3,1), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(1,3,2), WithinRel(2.0*v));
  REQUIRE_THAT(H(1,3,3), WithinAbs(0.0, 1e-13));

  REQUIRE_THAT(H(2,1,1), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(2,1,2), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(2,1,3), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(2,2,1), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(2,2,2), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(2,2,3), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(2,3,1), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(2,3,2), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(2,3,3), WithinRel(6.0*w));

  REQUIRE_THAT(H(3,1,1), WithinRel(2.0));
  REQUIRE_THAT(H(3,1,2), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(3,1,3), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(3,2,1), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(3,2,2), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(3,2,3), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(3,3,1), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(3,3,2), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(H(3,3,3), WithinAbs(0.0, 1e-13));

  Matrices dJdX = utl::jacobianGradient(Ji, H);

  REQUIRE_THAT(dJdX[0](1,1), WithinRel(H(1,1,1)*Ji(1,1) + H(1,1,2)*Ji(2,1) + H(1,1,3)*Ji(3,1)));
  REQUIRE_THAT(dJdX[0](1,2), WithinRel(H(1,1,2)*Ji(1,1) + H(1,2,2)*Ji(2,1) + H(1,2,3)*Ji(3,1)));
  REQUIRE_THAT(dJdX[0](1,3), WithinRel(H(1,1,3)*Ji(1,1) + H(1,2,3)*Ji(2,1) + H(1,3,3)*Ji(3,1)));
  REQUIRE_THAT(dJdX[0](2,1), WithinRel(H(2,1,1)*Ji(1,1) + H(2,1,2)*Ji(2,1) + H(2,1,3)*Ji(3,1)));
  REQUIRE_THAT(dJdX[0](2,2), WithinRel(H(2,1,2)*Ji(1,1) + H(2,2,2)*Ji(2,1) + H(2,2,3)*Ji(3,1)));
  REQUIRE_THAT(dJdX[0](2,3), WithinRel(H(2,1,3)*Ji(1,1) + H(2,2,3)*Ji(2,1) + H(2,3,3)*Ji(3,1)));
  REQUIRE_THAT(dJdX[0](3,1), WithinRel(H(3,1,1)*Ji(1,1) + H(3,1,2)*Ji(2,1) + H(3,1,3)*Ji(3,1)));
  REQUIRE_THAT(dJdX[0](3,2), WithinRel(H(3,1,2)*Ji(1,1) + H(3,2,2)*Ji(2,1) + H(3,2,3)*Ji(3,1)));
  REQUIRE_THAT(dJdX[0](3,3), WithinRel(H(3,1,3)*Ji(1,1) + H(3,2,3)*Ji(2,1) + H(3,3,3)*Ji(3,1)));

  REQUIRE_THAT(dJdX[1](1,1), WithinRel(H(1,1,1)*Ji(1,2) + H(1,1,2)*Ji(2,2) + H(1,1,3)*Ji(3,2)));
  REQUIRE_THAT(dJdX[1](1,2), WithinRel(H(1,1,2)*Ji(1,2) + H(1,2,2)*Ji(2,2) + H(1,2,3)*Ji(3,2)));
  REQUIRE_THAT(dJdX[1](1,3), WithinRel(H(1,1,3)*Ji(1,2) + H(1,2,3)*Ji(2,2) + H(1,3,3)*Ji(3,2)));
  REQUIRE_THAT(dJdX[1](2,1), WithinRel(H(2,1,1)*Ji(1,2) + H(2,1,2)*Ji(2,2) + H(2,1,3)*Ji(3,2)));
  REQUIRE_THAT(dJdX[1](2,2), WithinRel(H(2,1,2)*Ji(1,2) + H(2,2,2)*Ji(2,2) + H(2,2,3)*Ji(3,2)));
  REQUIRE_THAT(dJdX[1](2,3), WithinRel(H(2,1,3)*Ji(1,2) + H(2,2,3)*Ji(2,2) + H(2,3,3)*Ji(3,2)));
  REQUIRE_THAT(dJdX[1](3,1), WithinRel(H(3,1,1)*Ji(1,2) + H(3,1,2)*Ji(2,2) + H(3,1,3)*Ji(3,2)));
  REQUIRE_THAT(dJdX[1](3,2), WithinRel(H(3,1,2)*Ji(1,2) + H(3,2,2)*Ji(2,2) + H(3,2,3)*Ji(3,2)));
  REQUIRE_THAT(dJdX[1](3,3), WithinRel(H(3,1,3)*Ji(1,2) + H(3,2,3)*Ji(2,2) + H(3,3,3)*Ji(3,2)));

  REQUIRE_THAT(dJdX[2](1,1), WithinRel(H(1,1,1)*Ji(1,3) + H(1,1,2)*Ji(2,3) + H(1,1,3)*Ji(3,3)));
  REQUIRE_THAT(dJdX[2](1,2), WithinRel(H(1,1,2)*Ji(1,3) + H(1,2,2)*Ji(2,3) + H(1,2,3)*Ji(3,3)));
  REQUIRE_THAT(dJdX[2](1,3), WithinRel(H(1,1,3)*Ji(1,3) + H(1,2,3)*Ji(2,3) + H(1,3,3)*Ji(3,3)));
  REQUIRE_THAT(dJdX[2](2,1), WithinRel(H(2,1,1)*Ji(1,3) + H(2,1,2)*Ji(2,3) + H(2,1,3)*Ji(3,3)));
  REQUIRE_THAT(dJdX[2](2,2), WithinRel(H(2,1,2)*Ji(1,3) + H(2,2,2)*Ji(2,3) + H(2,2,3)*Ji(3,3)));
  REQUIRE_THAT(dJdX[2](2,3), WithinRel(H(2,1,3)*Ji(1,3) + H(2,2,3)*Ji(2,3) + H(2,3,3)*Ji(3,3)));
  REQUIRE_THAT(dJdX[2](3,1), WithinRel(H(3,1,1)*Ji(1,3) + H(3,1,2)*Ji(2,3) + H(3,1,3)*Ji(3,3)));
  REQUIRE_THAT(dJdX[2](3,2), WithinRel(H(3,1,2)*Ji(1,3) + H(3,2,2)*Ji(2,3) + H(3,2,3)*Ji(3,3)));
  REQUIRE_THAT(dJdX[2](3,3), WithinRel(H(3,1,3)*Ji(1,3) + H(3,2,3)*Ji(2,3) + H(3,3,3)*Ji(3,3)));
}


TEST_CASE("TestPiolaMapping.GradDetJ3D")
{
  ASMs3D p;
  std::stringstream g2(spline3D);

  REQUIRE(p.read(g2));
  REQUIRE(p.generateFEMTopology());

  double u = 0.1;
  double v = 0.1;
  double w = 0.1;

  Vector N;
  Matrix X, dNdu;
  Matrix J, Ji, dNdX;
  Matrix3D d2Ndu2, H;
  p.extractBasis(u, v, w, N, dNdu, d2Ndu2);
  p.getElementCoordinates(X, 1);
  utl::Jacobian(Ji, dNdX, X, dNdu, true);
  H.multiply(X,d2Ndu2);
  J.multiply(X,dNdu); // J = X * dNdu

  utl::vector<Real> ddet = utl::determinantGradient(J, Ji, H);
  utl::vector<Real> ddetu(3);

  ddetu(1) =  H(1,1,1) * (J(2,2) * J(3,3) - J(2,3) * J(3,2))
            + J(1,1) * (  H(2,2,1) * J(3,3) + J(2,2) * H(3,3,1)
                        - H(2,3,1)*J(3,2) - J(2,3) * H(3,2,1))
            - H(1,2,1) * (J(2,1)*J(3,3) - J(2,3) * J(3,1))
            - J(1,2) * (  H(2,1,1)*J(3,3) + J(2,1)*H(3,3,1)
                        - H(2,3,1)*J(3,1) - J(2,3)*H(3,1,1))
            + H(1,3,1) * (J(2,1)*J(3,2) - J(2,2)*J(3,1))
            + J(1,3) * (  H(2,1,1)*J(3,2) + J(2,1)*H(3,2,1)
                        - H(2,2,1)*J(3,1) - J(2,2)* H(3,1,1));
  ddetu(2) =  H(1,1,2) * (J(2,2) * J(3,3) - J(2,3) * J(3,2))
            + J(1,1) * (  H(2,2,2)*J(3,3) + J(2,2)*H(3,3,2)
                        - H(2,3,2)*J(3,2) - J(2,3)*H(3,2,2))
            - H(1,2,2) *  (J(2,1)*J(3,3) - J(2,3) * J(3,1))
            - J(1,2) * (  H(2,1,2)*J(3,3) + J(2,1) * H(3,3,2)
                        - H(2,3,2)*J(3,1) - J(2,3)*H(3,1,2))
            + H(1,3,2) * (J(2,1)*J(3,2) - J(2,2)*J(3,1))
            + J(1,3) * (  H(2,1,2) * J(3,2) + J(2,1) * H(3,2,2)
                        - H(2,2,2) * J(3,1) - J(2,2) * H(3,1,2));
  ddetu(3) =  H(1,1,3) * (J(2,2) * J(3,3) - J(2,3) * J(3,2))
            + J(1,1) * (  H(2,2,3)*J(3,3) + J(2,2)*H(3,3,3)
                        - H(2,3,3)*J(3,2) - J(2,3)*H(3,2,3))
            - H(1,2,3) * (J(2,1)*J(3,3) - J(2,3) * J(3,1))
            - J(1,2) * (  H(2,1,3)*J(3,3) + J(2,1) * H(3,3,3)
                        - H(2,3,3)*J(3,1) - J(2,3)*H(3,1,3))
            + H(1,3,3) * (J(2,1)*J(3,2) - J(2,2)*J(3,1))
            + J(1,3) * (  H(2,1,3)*J(3,2) + J(2,1)*H(3,2,3)
                        - H(2,2,3)*J(3,1) - J(2,2)*H(3,1,3));

  REQUIRE_THAT(ddet(1), WithinRel(ddetu(1) * Ji(1,1) + ddetu(2) * Ji(2,1) + ddetu(3) * Ji(3,1)));
  REQUIRE_THAT(ddet(2), WithinRel(ddetu(1) * Ji(1,2) + ddetu(2) * Ji(2,2) + ddetu(3) * Ji(3,2)));
  REQUIRE_THAT(ddet(3), WithinRel(ddetu(1) * Ji(1,3) + ddetu(2) * Ji(2,3) + ddetu(3) * Ji(3,3)));
}


namespace {

//! \brief Evaluates the covariant Piola mapping at a point of a 2D patch.
MxFiniteElement covariant2D (ASMs2D& p, ASMs2Dmx& b,
                             double u, double v, Matrix& J, Matrix& Nhat)
{
  BasisValues bfs(3);
  Go::BasisDerivsSf spline1, spline2;
  b.getBasis(1)->computeBasis(u,v,spline1);
  b.getBasis(2)->computeBasis(u,v,spline2);
  MxFiniteElement fe({spline1.basisValues.size(), spline2.basisValues.size()});
  SplineUtils::extractBasis(spline1, fe.basis(1), bfs[0].dNdu);
  SplineUtils::extractBasis(spline2, fe.basis(2), bfs[1].dNdu);

  Matrix X, Ji, dNdX;
  p.extractBasis(u, v, bfs.back().N, bfs.back().dNdu, bfs.back().d2Ndu2);
  p.getElementCoordinates(X, 1);
  fe.detJxW = utl::Jacobian(Ji, dNdX, X, bfs.back().dNdu, true);
  J.multiply(X,bfs.back().dNdu);
  fe.piolaMapping(Ji, X, bfs, false, true);

  Nhat.resize(2, fe.basis(1).size() + fe.basis(2).size());
  for (size_t i = 1; i <= fe.basis(1).size(); ++i)
    Nhat(1,i) = fe.basis(1)(i);
  for (size_t i = 1; i <= fe.basis(2).size(); ++i)
    Nhat(2,fe.basis(1).size()+i) = fe.basis(2)(i);

  return fe;
}


//! \brief Evaluates the covariant Piola mapping at a point of a 3D patch.
MxFiniteElement covariant3D (ASMs3D& p, ASMs3Dmx& b,
                             double u, double v, double w,
                             Matrix& J, Matrix& Nhat)
{
  BasisValues bfs(4);
  std::array<Go::BasisDerivs,3> spline;
  for (int d = 0; d < 3; ++d)
    b.getBasis(d+1)->computeBasis(u,v,w,spline[d]);
  MxFiniteElement fe({spline[0].basisValues.size(),
                      spline[1].basisValues.size(),
                      spline[2].basisValues.size()});
  for (int d = 0; d < 3; ++d)
    SplineUtils::extractBasis(spline[d], fe.basis(d+1), bfs[d].dNdu);

  Matrix X, Ji, dNdX;
  p.extractBasis(u, v, w, bfs.back().N, bfs.back().dNdu, bfs.back().d2Ndu2);
  p.getElementCoordinates(X, 1);
  fe.detJxW = utl::Jacobian(Ji, dNdX, X, bfs.back().dNdu, true);
  J.multiply(X,bfs.back().dNdu);
  fe.piolaMapping(Ji, X, bfs, false, true);

  size_t np = 0;
  for (int d = 1; d <= 3; ++d)
    np += fe.basis(d).size();
  Nhat.resize(3, np);
  size_t k = 1;
  for (int d = 1; d <= 3; ++d)
    for (size_t i = 1; i <= fe.basis(d).size(); ++i, ++k)
      Nhat(d,k) = fe.basis(d)(i);

  return fe;
}


/*!
  \brief Checks a covariant Piola mapping against its definition.
  \details The values must satisfy J^T N = Nhat, and the derivatives
  must match central differences of the values in each parametric direction.
*/

void checkCovariant (const std::function<MxFiniteElement(const RealArray&,
                                                         Matrix&, Matrix&)>& eval,
                     const RealArray& par)
{
  const size_t dim = par.size();
  Matrix J, Nhat;
  const MxFiniteElement fe = eval(par, J, Nhat);
  REQUIRE(fe.P_co.N.rows() == dim);
  REQUIRE(fe.P_co.dNdX.rows() == dim*dim);
  REQUIRE(fe.P_co.dNdX.cols() == fe.P_co.N.cols());

  Matrix JtN;
  JtN.multiply(J, fe.P_co.N, true);
  for (size_t k = 1; k <= Nhat.cols(); ++k)
    for (size_t j = 1; j <= dim; ++j)
      REQUIRE_THAT(JtN(j,k), WithinAbs(Nhat(j,k), 1e-12));

  constexpr double h = 1e-6;
  for (size_t p = 1; p <= dim; ++p) {
    RealArray parP(par), parM(par);
    parP[p-1] += h;
    parM[p-1] -= h;
    Matrix Jtmp, Ntmp;
    const MxFiniteElement feP = eval(parP, Jtmp, Ntmp);
    const MxFiniteElement feM = eval(parM, Jtmp, Ntmp);
    for (size_t k = 1; k <= Nhat.cols(); ++k)
      for (size_t j = 1; j <= dim; ++j) {
        // dN_j/du_p = sum_d dN_j/dX_d * dX_d/du_p
        double dNdu = 0.0;
        for (size_t d = 1; d <= dim; ++d)
          dNdu += fe.P_co.dNdX((d-1)*dim+j,k) * J(d,p);
        const double fd = (feP.P_co.N(j,k) - feM.P_co.N(j,k)) / (2.0*h);
        REQUIRE_THAT(dNdu, WithinAbs(fd, 1e-6));
      }
  }
}

}


TEST_CASE("TestPiolaMapping.Covariant2D")
{
  ASMs2D p;
  std::stringstream g2(spline2D);
  REQUIRE(p.read(g2));
  REQUIRE(p.generateFEMTopology());

  ASMs2Dmx b(2, {1,1,1});
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  std::stringstream basis(spline2D_basis);
  REQUIRE(b.read(basis));
  REQUIRE(b.generateFEMTopology());

  const auto eval = [&p,&b](const RealArray& par, Matrix& J, Matrix& Nhat)
  {
    return covariant2D(p, b, par[0], par[1], J, Nhat);
  };

  const RealArray par = GENERATE(RealArray{0.1, 0.1},
                                 RealArray{0.3, 0.8},
                                 RealArray{0.7, 0.4});
  checkCovariant(eval, par);
}


TEST_CASE("TestPiolaMapping.Covariant3D")
{
  ASMs3D p;
  std::stringstream g3(spline3D);
  REQUIRE(p.read(g3));
  REQUIRE(p.generateFEMTopology());

  ASMs3Dmx b({1,1,1});
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  std::stringstream basis(spline3D_basis);
  REQUIRE(b.read(basis));
  REQUIRE(b.generateFEMTopology());

  const auto eval = [&p,&b](const RealArray& par, Matrix& J, Matrix& Nhat)
  {
    return covariant3D(p, b, par[0], par[1], par[2], J, Nhat);
  };

  const RealArray par = GENERATE(RealArray{0.1, 0.1, 0.1},
                                 RealArray{0.3, 0.8, 0.5},
                                 RealArray{0.7, 0.4, 0.9});
  checkCovariant(eval, par);
}
