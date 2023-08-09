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

#include "gtest/gtest.h"

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"

#include <sstream>


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


TEST(TestPiolaMapping2D, Basis)
{
  ASMs2D p;
  std::stringstream g2(spline2D);

  ASSERT_TRUE(p.read(g2));
  ASSERT_TRUE(p.generateFEMTopology());

  ASMs2Dmx b(2, {1,1,1});
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  std::stringstream basis(spline2D_basis);
  ASSERT_TRUE(b.read(basis));
  ASSERT_TRUE(b.generateFEMTopology());

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
      EXPECT_NEAR(fe.P(1,k), (3.0*v + 4.0) * c[i](u) * q[j](v) / detJ, 1e-13);
      EXPECT_NEAR(fe.P(2,k),    (-v - 1.0) * c[i](u) * q[j](v) / detJ, 1e-13);
    }

  for (size_t j = 0; j < 4; ++j)
    for (size_t i = 0; i < 3; ++i, ++k) {
      EXPECT_NEAR(fe.P(1,k), (3.0*v - 1.0) * q[i](u) * c[j](v) / detJ, 1e-13);
      EXPECT_NEAR(fe.P(2,k),    (-u + 5.0) * q[i](u) * c[j](v) / detJ, 1e-13);
    }
}


TEST(TestPiolaMapping2D, Gradient)
{
  ASMs2D p;
  std::stringstream g2(spline2D);

  ASSERT_TRUE(p.read(g2));
  ASSERT_TRUE(p.generateFEMTopology());

  ASMs2Dmx b(2, {1,1,1});
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  std::stringstream basis(spline2D_basis);
  ASSERT_TRUE(b.read(basis));
  ASSERT_TRUE(b.generateFEMTopology());

  constexpr double u = 0.1;
  constexpr double v = 0.1;

  BasisValues bfs(3);

  Matrix X, J, Ji, dNdX;
  p.extractBasis(u, v, bfs.back().N, bfs.back().dNdu, bfs.back().d2Ndu2);
  p.getElementCoordinates(X, 1);
  Real detJ = utl::Jacobian(Ji, dNdX, X, bfs.back().dNdu, true);
  J.multiply(X,bfs.back().dNdu); // J = X * dNdu
  Matrix3D H;
  H.multiply(X,bfs.back().d2Ndu2);

  Go::BasisDerivsSf spline1, spline2;
  b.getBasis(1)->computeBasis(u,v,spline1);
  b.getBasis(2)->computeBasis(u,v,spline2);
  MxFiniteElement fe({spline1.basisValues.size(), spline2.basisValues.size()});
  SplineUtils::extractBasis(spline1, fe.basis(1), bfs[0].dNdu);
  SplineUtils::extractBasis(spline2, fe.basis(2), bfs[1].dNdu);

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

  fe.piolaMapping(detJ, Ji, X, bfs);
  std::vector<Real> det = utl::determinantGradient(J, Ji, H);
  Matrices dJdX = utl::jacobianGradient(Ji, H);

  size_t k = 1;
  for (size_t j = 0; j < 3; ++j)
    for (size_t i = 0; i < 4; ++i, ++k)
      for (size_t d1 = 1; d1 <= 2; ++d1)
        for (size_t d2 = 1; d2 <= 2; ++d2)
              EXPECT_NEAR(fe.dPdX((d1 + 2 * (d2 - 1)),k),
                          (-det[d2-1] / (detJ * detJ) * J(d1,1) +
                           dJdX[d2-1](d1,1) / detJ) * c[i](u)*q[j](v) +
                           J(d1,1) * (cd[i](u) * q[j](v) * Ji(1,d2) +
                                      c[i](u) * qd[j](v) * Ji(2,d2)) / detJ, 1e-13);

  for (size_t j = 0; j < 4; ++j)
    for (size_t i = 0; i < 3; ++i, ++k)
      for (size_t d1 = 1; d1 <= 2; ++d1)
        for (size_t d2 = 1; d2 <= 2; ++d2)
              EXPECT_NEAR(fe.dPdX(d1 + 2 * (d2 - 1),k),
                          (-det[d2-1] / (detJ * detJ) * J(d1,2) +
                          dJdX[d2-1](d1,2) / detJ) * q[i](u) * c[j](v) +
                          J(d1,2) * (qd[i](u) * c[j](v) * Ji(1,d2) +
                                     q[i](u) * cd[j](v) * Ji(2,d2)) / detJ, 1e-13);
}


TEST(TestPiolaMapping2D, GradJ)
{
  ASMs2D p;
  std::stringstream g2(spline2D);

  ASSERT_TRUE(p.read(g2));
  ASSERT_TRUE(p.generateFEMTopology());

  Vec3 Xp;
  double u = 0.1;
  double v = 0.1;
  SplineUtils::point(Xp, u, v, p.getBasis());
  EXPECT_DOUBLE_EQ(Xp[0],  3.0*u*v + 4.0*u - v);
  EXPECT_DOUBLE_EQ(Xp[1], -1.0*u*v - u + 5.0 * v);

  Vector N;
  Matrix X, dNdu;
  Matrix J, Ji, dNdX;
  Matrix3D d2Ndu2, H;
  p.extractBasis(u, v, N, dNdu, d2Ndu2);
  p.getElementCoordinates(X, 1);
  double detJ = utl::Jacobian(Ji, dNdX, X, dNdu, true);
  H.multiply(X,d2Ndu2);
  J.multiply(X,dNdu); // J = X * dNdu

  EXPECT_DOUBLE_EQ(detJ, -3.0*u*v - 4.0*u + 3.0*v*v + 17*v + 19);

  EXPECT_DOUBLE_EQ(J(1,1), 3.0*v + 4.0);
  EXPECT_DOUBLE_EQ(J(1,2), 3.0*v - 1.0);
  EXPECT_DOUBLE_EQ(J(2,1), -v - 1.0);
  EXPECT_DOUBLE_EQ(J(2,2), -u + 5.0);

  EXPECT_DOUBLE_EQ(Ji(1,1), J(2,2) / detJ);
  EXPECT_DOUBLE_EQ(Ji(1,2), -J(1,2) / detJ);
  EXPECT_DOUBLE_EQ(Ji(2,1), -J(2,1) / detJ);
  EXPECT_DOUBLE_EQ(Ji(2,2), J(1,1) / detJ);

  EXPECT_NEAR(H(1,1,1), 0.0, 1e-13);
  EXPECT_NEAR(H(1,1,2), 3.0, 1e-13);
  EXPECT_NEAR(H(1,2,1), 3.0, 1e-13);
  EXPECT_NEAR(H(1,2,2), 0.0, 1e-13);

  EXPECT_NEAR(H(2,1,1),  0.0, 1e-13);
  EXPECT_NEAR(H(2,1,2), -1.0, 1e-13);
  EXPECT_NEAR(H(2,2,1), -1.0, 1e-13);
  EXPECT_NEAR(H(2,2,2),  0.0, 1e-13);

  Matrices dJdX = utl::jacobianGradient(Ji, H);

  EXPECT_DOUBLE_EQ(dJdX[0](1,1), H(1,1,1)*Ji(1,1) + H(1,1,2)*Ji(2,1));
  EXPECT_DOUBLE_EQ(dJdX[0](1,2), H(1,1,2)*Ji(1,1) + H(1,2,2)*Ji(2,1));
  EXPECT_DOUBLE_EQ(dJdX[0](2,1), H(2,1,1)*Ji(1,1) + H(2,1,2)*Ji(2,1));
  EXPECT_DOUBLE_EQ(dJdX[0](2,2), H(2,1,2)*Ji(1,1) + H(2,2,2)*Ji(2,1));

  EXPECT_DOUBLE_EQ(dJdX[1](1,1), H(1,1,1)*Ji(1,2) + H(1,1,2)*Ji(2,2));
  EXPECT_DOUBLE_EQ(dJdX[1](1,2), H(1,1,2)*Ji(1,2) + H(1,2,2)*Ji(2,2));
  EXPECT_DOUBLE_EQ(dJdX[1](2,1), H(2,1,1)*Ji(1,2) + H(2,1,2)*Ji(2,2));
  EXPECT_DOUBLE_EQ(dJdX[1](2,2), H(2,1,2)*Ji(1,2) + H(2,2,2)*Ji(2,2));
}


TEST(TestPiolaMapping2D, GradDetJ)
{
  ASMs2D p;
  std::stringstream g2(spline2D);

  ASSERT_TRUE(p.read(g2));
  ASSERT_TRUE(p.generateFEMTopology());

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

  EXPECT_NEAR(det[0], Ju * Ji(1,1) + Jv * Ji(2,1), 1e-13);
  EXPECT_NEAR(det[1], Ju * Ji(1,2) + Jv * Ji(2,2), 1e-13);
}


TEST(TestPiolaMapping3D, Basis)
{
  ASMs3D p;
  std::stringstream g2(spline3D);

  ASSERT_TRUE(p.read(g2));
  ASSERT_TRUE(p.generateFEMTopology());

  ASMs3Dmx b({1,1,1});
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  std::stringstream basis(spline3D_basis);
  ASSERT_TRUE(b.read(basis));
  ASSERT_TRUE(b.generateFEMTopology());

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
        EXPECT_NEAR(fe.P(1,n), 2.0 * u * c[i](u) * q[j](v) * q[k](w) / detJ, 1e-13);
        EXPECT_NEAR(fe.P(2,n),     0.0 * c[i](u) * q[j](v) * q[k](w) / detJ, 1e-13);
        EXPECT_NEAR(fe.P(3,n), 2.0 * u * c[i](u) * q[j](v) * q[k](w) / detJ, 1e-13);
      }

  for (size_t k = 0; k < 3; ++k)
    for (size_t j = 0; j < 4; ++j)
      for (size_t i = 0; i < 3; ++i, ++n) {
        EXPECT_NEAR(fe.P(1,n), 2.0 * v * w * q[i](u) * c[j](v) * q[k](w) / detJ, 1e-13);
        EXPECT_NEAR(fe.P(2,n),         2.0 * q[i](u) * c[j](v) * q[k](w) / detJ, 1e-13);
        EXPECT_NEAR(fe.P(3,n),         0.0 * q[i](u) * c[j](v) * q[k](w) / detJ, 1e-13);
      }

  for (size_t k = 0; k < 4; ++k)
    for (size_t j = 0; j < 3; ++j)
      for (size_t i = 0; i < 3; ++i, ++n) {
        EXPECT_NEAR(fe.P(1,n),       v * v * q[i](u) * q[j](v) * c[k](w) / detJ, 1e-13);
        EXPECT_NEAR(fe.P(2,n), 3.0 * w * w * q[i](u) * q[j](v) * c[k](w) / detJ, 1e-13);
        EXPECT_NEAR(fe.P(3,n),         3.0 * q[i](u) * q[j](v) * c[k](w) / detJ, 1e-13);
      }
}


TEST(TestPiolaMapping3D, Gradient)
{
  ASMs3D p;
  std::stringstream g2(spline3D);

  ASSERT_TRUE(p.read(g2));
  ASSERT_TRUE(p.generateFEMTopology());

  ASMs3Dmx b({1,1,1});
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  std::stringstream basis(spline3D_basis);
  ASSERT_TRUE(b.read(basis));
  ASSERT_TRUE(b.generateFEMTopology());

  constexpr double u = 0.1;
  constexpr double v = 0.1;
  constexpr double w = 0.1;

  BasisValues bfs(4);

  Matrix X, J, Ji, dNdX;
  Matrix3D H;
  p.extractBasis(u, v, w, bfs.back().N, bfs.back().dNdu, bfs.back().d2Ndu2);
  p.getElementCoordinates(X, 1);
  Real detJ = utl::Jacobian(Ji, dNdX, X, bfs.back().dNdu, true);
  H.multiply(X,bfs.back().d2Ndu2);
  J.multiply(X,bfs.back().dNdu); // J = X * dNdu

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

  fe.piolaMapping(detJ, Ji, X, bfs);
  std::vector<Real> det = utl::determinantGradient(J, Ji, H);
  Matrices dJdX = utl::jacobianGradient(Ji, H);

  size_t n = 1;
  for (size_t k = 0; k < 3; ++k)
    for (size_t j = 0; j < 3; ++j)
      for (size_t i = 0; i < 4; ++i, ++n)
        for (size_t d1 = 1; d1 <= 3; ++d1)
          for (size_t d2 = 1; d2 <= 3; ++d2)
            EXPECT_NEAR(fe.dPdX(d1 + 3 * (d2-1),n),
                         (-det[d2-1] / (detJ * detJ) * J(d1,1) +
                         dJdX[d2-1](d1,1) / detJ) * c[i](u) * q[j](v) * q[k](w) +
                         J(d1,1) * (cd[i](u) * q[j](v) * q[k](w) * Ji(1,d2) +
                                    c[i](u) * qd[j](v) * q[k](w) * Ji(2,d2) +
                                    c[i](u) * q[j](v) * qd[k](w) * Ji(3,d2)) / detJ, 1e-13);

  for (size_t k = 0; k < 3; ++k)
    for (size_t j = 0; j < 4; ++j)
      for (size_t i = 0; i < 3; ++i, ++n)
        for (size_t d1 = 1; d1 <= 3; ++d1)
          for (size_t d2 = 1; d2 <= 3; ++d2)
            EXPECT_NEAR(fe.dPdX(d1 + 3 * (d2-1),n),
                         (-det[d2-1] / (detJ * detJ) * J(d1,2) +
                         dJdX[d2-1](d1,2) / detJ) * q[i](u) * c[j](v) * q[k](w) +
                         J(d1,2) * (qd[i](u) * c[j](v) * q[k](w) * Ji(1,d2) +
                                    q[i](u) * cd[j](v) * q[k](w) * Ji(2,d2) +
                                    q[i](u) * c[j](v) * qd[k](w) * Ji(3,d2)) / detJ, 1e-13);

  for (size_t k = 0; k < 4; ++k)
    for (size_t j = 0; j < 3; ++j)
      for (size_t i = 0; i < 3; ++i, ++n)
        for (size_t d1 = 1; d1 <= 3; ++d1)
          for (size_t d2 = 1; d2 <= 3; ++d2)
            EXPECT_NEAR(fe.dPdX(d1 + 3 * (d2-1),n),
                         (-det[d2-1] / (detJ * detJ) * J(d1,3) +
                         dJdX[d2-1](d1,3) / detJ) * q[i](u) * q[j](v) * c[k](w) +
                         J(d1,3) * (qd[i](u) * q[j](v) * c[k](w) * Ji(1,d2) +
                                    q[i](u) * qd[j](v) * c[k](w) * Ji(2,d2) +
                                    q[i](u) * q[j](v) * cd[k](w) * Ji(3,d2)) / detJ, 1e-13);
}


TEST(TestPiolaMapping3D, GradJ)
{
  ASMs3D p;
  std::stringstream g2(spline3D);

  ASSERT_TRUE(p.read(g2));
  ASSERT_TRUE(p.generateFEMTopology());

  Vec3 Xp;
  double u = 0.1;
  double v = 0.1;
  double w = 0.1;
  SplineUtils::point(Xp, u, v, w, p.getBasis());
  EXPECT_DOUBLE_EQ(Xp[0], u * u + v * v * w);
  EXPECT_DOUBLE_EQ(Xp[1], 2.0 * v + w * w * w);
  EXPECT_DOUBLE_EQ(Xp[2], 3.0 * w + u * u);

  Vector N;
  Matrix X, dNdu;
  Matrix J, Ji, dNdX;
  Matrix3D d2Ndu2, H;
  p.extractBasis(u, v, w, N, dNdu, d2Ndu2);
  p.getElementCoordinates(X, 1);
  double detJ = utl::Jacobian(Ji, dNdX, X, dNdu, true);
  H.multiply(X,d2Ndu2);
  J.multiply(X,dNdu); // J = X * dNdu

  EXPECT_NEAR(J(1,1), 2.0 * u, 1e-13);
  EXPECT_NEAR(J(1,2), 2.0 * v * w, 1e-13);
  EXPECT_NEAR(J(1,3), v * v, 1e-13);
  EXPECT_NEAR(J(2,1), 0.0, 1e-13);
  EXPECT_NEAR(J(2,2), 2.0, 1e-13);
  EXPECT_NEAR(J(2,3), 3.0 * w * w, 1e-13);
  EXPECT_NEAR(J(3,1), 2.0 * u, 1e-13);
  EXPECT_NEAR(J(3,2), 0.0, 1e-13);
  EXPECT_NEAR(J(3,3), 3.0, 1e-13);

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

  EXPECT_NEAR(detRef, detJ, 1e-13);

  EXPECT_DOUBLE_EQ(Ji(1,1),  det(J(2,2), J(1,3), J(3,2), J(3,3)) / detJ);
  EXPECT_DOUBLE_EQ(Ji(1,2), -det(J(1,2), J(1,3), J(3,2), J(3,3)) / detJ);
  EXPECT_DOUBLE_EQ(Ji(1,3),  det(J(1,2), J(2,2), J(1,3), J(2,3)) / detJ);
  EXPECT_DOUBLE_EQ(Ji(2,1), -det(J(2,1), J(2,3), J(3,1), J(3,3)) / detJ);
  EXPECT_DOUBLE_EQ(Ji(2,2),  det(J(1,1), J(1,3), J(3,1), J(3,3)) / detJ);
  EXPECT_DOUBLE_EQ(Ji(2,3), -det(J(1,1), J(1,3), J(2,1), J(2,3)) / detJ);
  EXPECT_DOUBLE_EQ(Ji(3,1),  det(J(2,1), J(2,2), J(3,1), J(3,2)) / detJ);
  EXPECT_DOUBLE_EQ(Ji(3,2), -det(J(1,1), J(1,2), J(3,1), J(3,2)) / detJ);
  EXPECT_DOUBLE_EQ(Ji(3,3),  det(J(1,1), J(1,2), J(2,1), J(2,2)) / detJ);

  EXPECT_NEAR(H(1,1,1),   2.0, 1e-13);
  EXPECT_NEAR(H(1,1,2),   0.0, 1e-13);
  EXPECT_NEAR(H(1,1,3),   0.0, 1e-13);
  EXPECT_NEAR(H(1,2,1),   0.0, 1e-13);
  EXPECT_NEAR(H(1,2,2), 2.0*w, 1e-13);
  EXPECT_NEAR(H(1,2,3), 2.0*v, 1e-13);
  EXPECT_NEAR(H(1,3,1),   0.0, 1e-13);
  EXPECT_NEAR(H(1,3,2), 2.0*v, 1e-13);
  EXPECT_NEAR(H(1,3,3),   0.0, 1e-13);

  EXPECT_NEAR(H(2,1,1),   0.0, 1e-13);
  EXPECT_NEAR(H(2,1,2),   0.0, 1e-13);
  EXPECT_NEAR(H(2,1,3),   0.0, 1e-13);
  EXPECT_NEAR(H(2,2,1),   0.0, 1e-13);
  EXPECT_NEAR(H(2,2,2),   0.0, 1e-13);
  EXPECT_NEAR(H(2,2,3),   0.0, 1e-13);
  EXPECT_NEAR(H(2,3,1),   0.0, 1e-13);
  EXPECT_NEAR(H(2,3,2),   0.0, 1e-13);
  EXPECT_NEAR(H(2,3,3), 6.0*w, 1e-13);

  EXPECT_NEAR(H(3,1,1), 2.0, 1e-13);
  EXPECT_NEAR(H(3,1,2), 0.0, 1e-13);
  EXPECT_NEAR(H(3,1,3), 0.0, 1e-13);
  EXPECT_NEAR(H(3,2,1), 0.0, 1e-13);
  EXPECT_NEAR(H(3,2,2), 0.0, 1e-13);
  EXPECT_NEAR(H(3,2,3), 0.0, 1e-13);
  EXPECT_NEAR(H(3,3,1), 0.0, 1e-13);
  EXPECT_NEAR(H(3,3,2), 0.0, 1e-13);
  EXPECT_NEAR(H(3,3,3), 0.0, 1e-13);

  Matrices dJdX = utl::jacobianGradient(Ji, H);

  EXPECT_DOUBLE_EQ(dJdX[0](1,1), H(1,1,1)*Ji(1,1) + H(1,1,2)*Ji(2,1) + H(1,1,3)*Ji(3,1));
  EXPECT_DOUBLE_EQ(dJdX[0](1,2), H(1,1,2)*Ji(1,1) + H(1,2,2)*Ji(2,1) + H(1,2,3)*Ji(3,1));
  EXPECT_DOUBLE_EQ(dJdX[0](1,3), H(1,1,3)*Ji(1,1) + H(1,2,3)*Ji(2,1) + H(1,3,3)*Ji(3,1));
  EXPECT_DOUBLE_EQ(dJdX[0](2,1), H(2,1,1)*Ji(1,1) + H(2,1,2)*Ji(2,1) + H(2,1,3)*Ji(3,1));
  EXPECT_DOUBLE_EQ(dJdX[0](2,2), H(2,1,2)*Ji(1,1) + H(2,2,2)*Ji(2,1) + H(2,2,3)*Ji(3,1));
  EXPECT_DOUBLE_EQ(dJdX[0](2,3), H(2,1,3)*Ji(1,1) + H(2,2,3)*Ji(2,1) + H(2,3,3)*Ji(3,1));
  EXPECT_DOUBLE_EQ(dJdX[0](3,1), H(3,1,1)*Ji(1,1) + H(3,1,2)*Ji(2,1) + H(3,1,3)*Ji(3,1));
  EXPECT_DOUBLE_EQ(dJdX[0](3,2), H(3,1,2)*Ji(1,1) + H(3,2,2)*Ji(2,1) + H(3,2,3)*Ji(3,1));
  EXPECT_DOUBLE_EQ(dJdX[0](3,3), H(3,1,3)*Ji(1,1) + H(3,2,3)*Ji(2,1) + H(3,3,3)*Ji(3,1));

  EXPECT_DOUBLE_EQ(dJdX[1](1,1), H(1,1,1)*Ji(1,2) + H(1,1,2)*Ji(2,2) + H(1,1,3)*Ji(3,2));
  EXPECT_DOUBLE_EQ(dJdX[1](1,2), H(1,1,2)*Ji(1,2) + H(1,2,2)*Ji(2,2) + H(1,2,3)*Ji(3,2));
  EXPECT_DOUBLE_EQ(dJdX[1](1,3), H(1,1,3)*Ji(1,2) + H(1,2,3)*Ji(2,2) + H(1,3,3)*Ji(3,2));
  EXPECT_DOUBLE_EQ(dJdX[1](2,1), H(2,1,1)*Ji(1,2) + H(2,1,2)*Ji(2,2) + H(2,1,3)*Ji(3,2));
  EXPECT_DOUBLE_EQ(dJdX[1](2,2), H(2,1,2)*Ji(1,2) + H(2,2,2)*Ji(2,2) + H(2,2,3)*Ji(3,2));
  EXPECT_DOUBLE_EQ(dJdX[1](2,3), H(2,1,3)*Ji(1,2) + H(2,2,3)*Ji(2,2) + H(2,3,3)*Ji(3,2));
  EXPECT_DOUBLE_EQ(dJdX[1](3,1), H(3,1,1)*Ji(1,2) + H(3,1,2)*Ji(2,2) + H(3,1,3)*Ji(3,2));
  EXPECT_DOUBLE_EQ(dJdX[1](3,2), H(3,1,2)*Ji(1,2) + H(3,2,2)*Ji(2,2) + H(3,2,3)*Ji(3,2));
  EXPECT_DOUBLE_EQ(dJdX[1](3,3), H(3,1,3)*Ji(1,2) + H(3,2,3)*Ji(2,2) + H(3,3,3)*Ji(3,2));

  EXPECT_DOUBLE_EQ(dJdX[2](1,1), H(1,1,1)*Ji(1,3) + H(1,1,2)*Ji(2,3) + H(1,1,3)*Ji(3,3));
  EXPECT_DOUBLE_EQ(dJdX[2](1,2), H(1,1,2)*Ji(1,3) + H(1,2,2)*Ji(2,3) + H(1,2,3)*Ji(3,3));
  EXPECT_DOUBLE_EQ(dJdX[2](1,3), H(1,1,3)*Ji(1,3) + H(1,2,3)*Ji(2,3) + H(1,3,3)*Ji(3,3));
  EXPECT_DOUBLE_EQ(dJdX[2](2,1), H(2,1,1)*Ji(1,3) + H(2,1,2)*Ji(2,3) + H(2,1,3)*Ji(3,3));
  EXPECT_DOUBLE_EQ(dJdX[2](2,2), H(2,1,2)*Ji(1,3) + H(2,2,2)*Ji(2,3) + H(2,2,3)*Ji(3,3));
  EXPECT_DOUBLE_EQ(dJdX[2](2,3), H(2,1,3)*Ji(1,3) + H(2,2,3)*Ji(2,3) + H(2,3,3)*Ji(3,3));
  EXPECT_DOUBLE_EQ(dJdX[2](3,1), H(3,1,1)*Ji(1,3) + H(3,1,2)*Ji(2,3) + H(3,1,3)*Ji(3,3));
  EXPECT_DOUBLE_EQ(dJdX[2](3,2), H(3,1,2)*Ji(1,3) + H(3,2,2)*Ji(2,3) + H(3,2,3)*Ji(3,3));
  EXPECT_DOUBLE_EQ(dJdX[2](3,3), H(3,1,3)*Ji(1,3) + H(3,2,3)*Ji(2,3) + H(3,3,3)*Ji(3,3));
}


TEST(TestPiolaMapping3D, GradDetJ)
{
  ASMs3D p;
  std::stringstream g2(spline3D);

  ASSERT_TRUE(p.read(g2));
  ASSERT_TRUE(p.generateFEMTopology());

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

  EXPECT_NEAR(ddet(1), ddetu(1) * Ji(1,1) + ddetu(2) * Ji(2,1) + ddetu(3) * Ji(3,1), 1e-13);
  EXPECT_NEAR(ddet(2), ddetu(1) * Ji(1,2) + ddetu(2) * Ji(2,2) + ddetu(3) * Ji(3,2), 1e-13);
  EXPECT_NEAR(ddet(3), ddetu(1) * Ji(1,3) + ddetu(2) * Ji(2,3) + ddetu(3) * Ji(3,3), 1e-13);
}
