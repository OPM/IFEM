//==============================================================================
//!
//! \file TestCoordinateMapping.C
//!
//! \date Oct 9 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for utilities for coordinate mapping transformations.
//!
//==============================================================================

#include "CoordinateMapping.h"
#include "ASMs1D.h"
#include "ASMs2D.h"
#include "ASMs2Dmx.h"
#include "ASMs3D.h"
#include "SplineUtils.h"

#include "GoTools/geometry/SplineSurface.h"

#include <sstream>
#include <array>

#include "gtest/gtest.h"


static const char* spline1D = "100 1 0 0\n1 0\n2 2 0 0 1 1\n0 1\n";
static const char* spline2D = "200 1 0 0\n3 0\n"
  "2 2 0 0 1 1\n"
  "2 2 0 0 1 1\n"
  "0 0 0 1 0 0 0 1 0 1 1 0\n";
static const char* spline3D = "700 1 0 0\n3 0\n"
  "2 2 0 0 1 1\n"
  "2 2 0 0 1 1\n"
  "2 2 0 0 1 1\n"
  "0 0 0 1 0 0 0 1 0 1 1 0\n"
  "0 0 1 1 0 1 0 1 1 1 1 1\n";


TEST(TestCoordinateMapping, Jacobian1D)
{
  ASMs1D p;
  std::stringstream g2(spline1D);

  ASSERT_TRUE(p.read(g2));
  ASSERT_TRUE(p.uniformRefine(3));
  ASSERT_TRUE(p.generateFEMTopology());

  Vector N;
  Matrix X, dNdU;
  p.extractBasis(0.25, N, dNdU);
  p.getElementCoordinates(X, 2);

  Matrix J, dNdX;
  double det = utl::Jacobian(J, dNdX, X, dNdU, true);

  EXPECT_FLOAT_EQ(det, 1.0);
  EXPECT_FLOAT_EQ(J(1,1), 1.0);
  EXPECT_EQ(dNdX.rows(), 2U);
  EXPECT_EQ(dNdX.cols(), 1U);
  EXPECT_FLOAT_EQ(dNdX(1,1),-4.0);
  EXPECT_FLOAT_EQ(dNdX(2,1), 4.0);
}


TEST(TestCoordinateMapping, Hessian1D)
{
  ASMs1D p;
  std::stringstream g2(spline1D);

  ASSERT_TRUE(p.read(g2));
  ASSERT_TRUE(p.raiseOrder(2));
  ASSERT_TRUE(p.uniformRefine(3));
  ASSERT_TRUE(p.generateFEMTopology());

  Vector N;
  Matrix X, dNdU;
  Matrix3D d2Ndu2;
  p.extractBasis(0.25, N, dNdU, d2Ndu2);
  p.getElementCoordinates(X, 2);

  Matrix J, dNdX;
  Matrix3D H, d2NdX2;
  double det = utl::Jacobian(J, dNdX, X, dNdU, true);
  utl::Hessian(H, d2NdX2, J, X, d2Ndu2, dNdX);

  const double dndx[] = {
    -3.0, 1.0, 2.0, 0.0
  };

  const double d2ndx2[] = {
    24.0,-40.0, 16.0, 0.0
  };

  EXPECT_FLOAT_EQ(det, 1.0);
  EXPECT_FLOAT_EQ(J(1,1), 1.0);
  EXPECT_NEAR(H(1,1,1), 0.0, 1.0e-13);

  size_t i = 0;
  EXPECT_EQ(dNdX.rows(), 4U);
  EXPECT_EQ(dNdX.cols(), 1U);
  for (double v : dNdX)
    EXPECT_FLOAT_EQ(dndx[i++], v);

  i = 0;
  EXPECT_EQ(d2NdX2.dim(1), 4U);
  EXPECT_EQ(d2NdX2.dim(2), 1U);
  EXPECT_EQ(d2NdX2.dim(3), 1U);
  for (double v : d2NdX2)
    EXPECT_FLOAT_EQ(d2ndx2[i++], v);
}


TEST(TestCoordinateMapping, Jacobian2D)
{
  ASMs2D p;
  std::stringstream g2(spline2D);

  ASSERT_TRUE(p.read(g2));
  ASSERT_TRUE(p.uniformRefine(0,3));
  ASSERT_TRUE(p.uniformRefine(1,3));
  ASSERT_TRUE(p.generateFEMTopology());

  Vector N;
  Matrix X, dNdU;
  p.extractBasis(0.25, 0.25, N, dNdU);
  p.getElementCoordinates(X, 2);

  Matrix J, dNdX;
  double det = utl::Jacobian(J, dNdX, X, dNdU, true);

  const double dndx[] = {
    -4.0, 4.0, 0.0, 0.0,
    -4.0, 0.0, 4.0, 0.0
  };

  EXPECT_FLOAT_EQ(det, 1.0);
  for (size_t i = 1; i <= J.rows(); i++)
    for (size_t j = 1; j <= J.rows(); j++)
      EXPECT_FLOAT_EQ(J(i,j), i == j ? 1.0 : 0.0);

  size_t i = 0;
  EXPECT_EQ(dNdX.rows(), 4U);
  EXPECT_EQ(dNdX.cols(), 2U);
  for (double v : dNdX)
    EXPECT_FLOAT_EQ(dndx[i++], v);
}


TEST(TestCoordinateMapping, Hessian2D)
{
  ASMs2D p;
  std::stringstream g2(spline2D);

  ASSERT_TRUE(p.read(g2));
  ASSERT_TRUE(p.raiseOrder(2,2));
  ASSERT_TRUE(p.uniformRefine(0,3));
  ASSERT_TRUE(p.uniformRefine(1,3));
  ASSERT_TRUE(p.generateFEMTopology());

  Vector N;
  Matrix X, dNdU;
  Matrix3D d2Ndu2;
  p.extractBasis(0.25, 0.25, N, dNdU, d2Ndu2);
  p.getElementCoordinates(X, 2);

  Matrix J, dNdX;
  Matrix3D H, d2NdX2;
  double det = utl::Jacobian(J, dNdX, X, dNdU, true);
  utl::Hessian(H, d2NdX2, J, X, d2Ndu2, dNdX);

  const double dndx[] = {
    -0.75, 0.25, 0.5, 0.0,-1.75, 0.5833333333333, 1.1666666666667, 0.0,
    -0.5, 0.1666666666667, 0.3333333333333, 0.0, 0.0, 0.0, 0.0, 0.0,

    -1.2857142857143, -3.0, -0.8571428571428, 0.0,
    0.4285714285714, 1.0, 0.2857142857143, 0.0,
    0.8571428571428, 2.0, 0.5714285714285, 0.0,
    0.0, 0.0, 0.0, 0.0
  };

  const double d2ndx2[] = {
    6.0, -10.0, 4.0, 0.0, 14.0, -23.333333333333, 9.333333333333, 0.0,
    4.0, -6.666666666667, 2.666666666667, 0.0, 0.0, 0.0, 0.0, 0.0,

    15.428571428571, -5.142857142857, -10.285714285714, 0.0,
    -5.142857142857, 1.714285714286, 3.428571428571, 0.0,
    -10.285714285714, 3.428571428571, 6.857142857143, 0.0,  0.0, 0.0, 0.0, 0.0,

    15.428571428571, -5.142857142857, -10.285714285714, 0.0,
    -5.142857142857, 1.714285714286, 3.428571428571, 0.0,
    -10.285714285714, 3.428571428571, 6.857142857143, 0.0,  0.0, 0.0, 0.0, 0.0,

    20.151603498542, 47.020408163265, 13.434402332362, 0.0,
    -30.227405247813, -70.530612244897, -20.151603498542, 0.0,
    10.075801749271, 23.510204081632, 6.717201166181, 0.0,  0.0, 0.0, 0.0, 0.0
  };

  EXPECT_FLOAT_EQ(det, 7.0/12.0);
  EXPECT_NEAR(J(1,1), 1.0, 1.0e-13);
  EXPECT_NEAR(J(2,2), 1.7142857142857, 1.0e-13);
  EXPECT_NEAR(J(1,2), 0.0, 1.0e-13);
  EXPECT_NEAR(J(2,1), 0.0, 1.0e-13);

  for (size_t i = 1; i <= H.dim(1); i++)
    for (size_t j = 1; j <= H.dim(2); j++)
      for (size_t k = 1; k <= H.dim(3); k++)
        EXPECT_NEAR(H(i,j,k), i==j && j==k && k==2 ? 2.0/3.0 : 0.0, 1.0e-13);

  size_t i = 0;
  EXPECT_EQ(dNdX.rows(), 16U);
  EXPECT_EQ(dNdX.cols(), 2U);
  for (double v : dNdX)
    EXPECT_NEAR(dndx[i++], v, 1.0e-13);

  i = 0;
  EXPECT_EQ(d2NdX2.dim(1), 16U);
  EXPECT_EQ(d2NdX2.dim(2), 2U);
  EXPECT_EQ(d2NdX2.dim(3), 2U);
  for (double v : d2NdX2)
    EXPECT_NEAR(d2ndx2[i++], v, 1.0e-12);
}


TEST(TestCoordinateMapping, JacobianShell)
{
  ASMs2D p(3);
  std::stringstream g2(spline2D);

  ASSERT_TRUE(p.read(g2));
  ASSERT_TRUE(p.raiseOrder(1,1));
  ASSERT_TRUE(p.uniformRefine(0,3));
  ASSERT_TRUE(p.uniformRefine(1,3));
  ASSERT_TRUE(p.generateFEMTopology());

  Vector N;
  Matrix X, dNdU;
  Matrix3D d2Ndu2;
  p.extractBasis(0.25, 0.25, N, dNdU, d2Ndu2);
  p.getElementCoordinates(X, 2);

  Matrix J, dNdX, Hs;
  Matrix3D H, d2NdX2;
  double det = utl::Jacobian(J, dNdX, X, dNdU, true);
  utl::Hessian(H, d2NdX2, J, X, d2Ndu2, dNdX);
  utl::Hessian(H, Hs);

  const double j[] = {
    1.0, 0.0, 0.0,
    0.0, 0.5, 0.0
  };

  const double h[] = {
    0.0, 0.0, 0.0,
    0.0, 2.0, 0.0,
    0.0, 0.0, 0.0
  };

  const double dndx[] = {
    -2.0, 2.0, 0.0,-2.0, 2.0, 0.0, 0.0, 0.0, 0.0,
    -2.0,-2.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0
  };

  EXPECT_FLOAT_EQ(det, 0.5);

  size_t i = 0;
  EXPECT_EQ(J.rows(), 3U);
  EXPECT_EQ(J.cols(), 2U);
  for (double v : J)
    EXPECT_FLOAT_EQ(j[i++], v);

  i = 0;
  EXPECT_EQ(dNdX.rows(), 9U);
  EXPECT_EQ(dNdX.cols(), 2U);
  for (double v : dNdX)
    EXPECT_FLOAT_EQ(dndx[i++], v);

  i = 0;
  EXPECT_EQ(Hs.rows(), 3U);
  EXPECT_EQ(Hs.cols(), 3U);
  for (double v : Hs)
    EXPECT_FLOAT_EQ(h[i++], v);
}


TEST(TestCoordinateMapping, Hessian2D_mixed)
{
  ASMs2Dmx p(2,{1,1});
  std::stringstream g2(spline2D);

  ASSERT_TRUE(p.read(g2));
  ASSERT_TRUE(p.generateFEMTopology());

  std::array<std::vector<Go::BasisDerivsSf2>, 2> splinex2;
  p.getBasis(1)->computeBasisGrid({0.5}, {0.5}, splinex2[0]);
  p.getBasis(2)->computeBasisGrid({0.5}, {0.5}, splinex2[1]);

  std::array<Vector, 2> N;
  std::array<Matrix, 2> dNxdu;
  std::array<Matrix3D, 2> d2Nxdu2;
  SplineUtils::extractBasis(splinex2[0][0], N[0], dNxdu[0], d2Nxdu2[0]);
  SplineUtils::extractBasis(splinex2[1][0], N[1], dNxdu[1], d2Nxdu2[1]);

  Matrix Jac;
  std::array<Matrix, 2> grad;

  int itgBasis = ASMmxBase::itgBasis - 1;

  Matrix Xnod;
  p.getElementCoordinates(Xnod,1);

  utl::Jacobian(Jac, grad[itgBasis], Xnod, dNxdu[itgBasis]);

  grad[1-itgBasis].multiply(dNxdu[1-itgBasis],Jac);

  const double Hess_basis1[] = {
    0.5,-1.0, 0.5, 1.0,-2.0, 1.0, 0.5,-1.0, 0.5,
    1.0, 0.0,-1.0, 0.0, 0.0, 0.0,-1.0, 0.0, 1.0,
    1.0, 0.0,-1.0, 0.0, 0.0, 0.0,-1.0, 0.0, 1.0,
    0.5, 1.0, 0.5,-1.0,-2.0,-1.0, 0.5, 1.0, 0.5
  };

  const double Hess_basis2[] = {
    0.0, 0.0, 0.0, 0.0,
    1.0,-1.0,-1.0, 1.0,
    1.0,-1.0,-1.0, 1.0,
    0.0, 0.0, 0.0, 0.0,
  };

  size_t i;
  Matrix3D Hess;
  std::array<Matrix3D, 2> hess;
  for (int b = 1; b >= 0; b--)
  {
    utl::Hessian(Hess,hess[b],Jac,Xnod,d2Nxdu2[b],grad[b],b==1);
    // geometry mapping should be the identify mapping
    const double* h = hess[b].ptr();
    i = 0;
    for (double v : d2Nxdu2[b])
      EXPECT_FLOAT_EQ(h[i++], v);
  }

  i = 0;
  for (double v : hess[0])
    EXPECT_FLOAT_EQ(Hess_basis1[i++],v);

  i = 0;
  for (double v : hess[1])
    EXPECT_FLOAT_EQ(Hess_basis2[i++],v);
}


TEST(TestCoordinateMapping, Jacobian3D)
{
  ASMs3D p;
  std::stringstream g2(spline3D);

  ASSERT_TRUE(p.read(g2));
  ASSERT_TRUE(p.uniformRefine(0,3));
  ASSERT_TRUE(p.uniformRefine(1,3));
  ASSERT_TRUE(p.uniformRefine(2,3));
  ASSERT_TRUE(p.generateFEMTopology());

  Vector N;
  Matrix X, dNdU;
  p.extractBasis(0.25, 0.25, 0.25, N, dNdU);
  p.getElementCoordinates(X, 2);

  Matrix J, dNdX;
  double det = utl::Jacobian(J, dNdX, X, dNdU, true);

  const double dndx[] = {
    -4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -4.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -4.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0
  };

  EXPECT_FLOAT_EQ(det, 1.0);
  for (size_t i = 1; i <= J.rows(); i++)
    for (size_t j = 1; j <= J.rows(); j++)
      EXPECT_FLOAT_EQ(J(i,j), i == j ? 1.0 : 0.0);

  size_t i = 0;
  EXPECT_EQ(dNdX.rows(), 8U);
  EXPECT_EQ(dNdX.cols(), 3U);
  for (double v : dNdX)
    EXPECT_FLOAT_EQ(dndx[i++], v);
}


TEST(TestCoordinateMapping, Hessian3D)
{
  ASMs3D p;
  std::stringstream g2(spline3D);

  ASSERT_TRUE(p.read(g2));
  ASSERT_TRUE(p.raiseOrder(2,2,2));
  ASSERT_TRUE(p.uniformRefine(0,3));
  ASSERT_TRUE(p.uniformRefine(1,3));
  ASSERT_TRUE(p.uniformRefine(2,3));
  ASSERT_TRUE(p.generateFEMTopology());

  Vector N;
  Matrix X, dNdU;
  Matrix3D d2Ndu2;
  p.extractBasis(0.25, 0.25, 0.25, N, dNdU, d2Ndu2);
  p.getElementCoordinates(X, 2);

  Matrix J, dNdX;
  Matrix3D H, d2NdX2;
  double det = utl::Jacobian(J, dNdX, X, dNdU, true);
  utl::Hessian(H, d2NdX2, J, X, d2Ndu2, dNdX);

  size_t i, j, k;
  EXPECT_FLOAT_EQ(det, 0.34027779);
  for (i = 1; i <= J.rows(); i++)
    for (j = 1; j <= J.rows(); j++)
      if (i == 1 && j == 1)
        EXPECT_FLOAT_EQ(J(i,j), 1.0);
      else
        EXPECT_NEAR(J(i,j), i == j ? 1.7142857142857 : 0.0, 1.0e-13);

  const double dndx[] = {
    -0.1875, 0.0625, 0.125, 0, -0.4375, 0.1458333333333,
    0.2916666666667, 0, -0.125, 0.04166666666667, 0.08333333333333, 0,
    0, 0, 0, 0, -0.4375, 0.1458333333333,
    0.2916666666667, 0, -1.0208333333333, 0.3402777777778, 0.6805555555556, 0,
    -0.2916666666667, 0.09722222222222, 0.1944444444444, 0, 0, 0,
    0, 0, -0.125, 0.04166666666667, 0.08333333333333, 0,
    -0.2916666666667, 0.09722222222222, 0.1944444444444, 0, -0.08333333333333, 0.02777777777778,
    0.05555555555556, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0,

    -0.3214285714286, -0.75, -0.2142857142857, 0, 0.1071428571429, 0.25,
    0.07142857142857, 0, 0.2142857142857, 0.5, 0.1428571428571, 0,
    0, 0, 0, 0, -0.75, -1.75,
    -0.5, 0, 0.25, 0.5833333333333, 0.1666666666667, 0,
    0.5, 1.1666666666667, 0.3333333333333, 0, 0, 0,
    0, 0, -0.2142857142857, -0.5, -0.1428571428571, 0,
    0.07142857142857, 0.1666666666667, 0.04761904761905, 0, 0.1428571428571, 0.3333333333333,
    0.0952380952381, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0,

    -0.3214285714286, -0.75, -0.2142857142857, 0, -0.75, -1.75,
    -0.5, 0, -0.2142857142857, -0.5, -0.1428571428571, 0,
    0, 0, 0, 0, 0.1071428571429, 0.25,
    0.07142857142857, 0, 0.25, 0.5833333333333, 0.1666666666667, 0,
    0.07142857142857, 0.1666666666667, 0.04761904761905, 0, 0, 0,
    0, 0, 0.2142857142857, 0.5, 0.1428571428571, 0,
    0.5, 1.1666666666667, 0.3333333333333, 0, 0.1428571428571, 0.3333333333333,
    0.0952380952381, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0.0
  };

  const double d2ndx2[] = {
    1.5, -2.5, 1, 0, 3.5, -5.833333333333,
    2.333333333333, 0, 1, -1.666666666667, 0.6666666666667, 0,
    0, 0, 0, 0, 3.5, -5.833333333333,
    2.333333333333, 0, 8.166666666667, -13.611111111111, 5.444444444444, 0,
    2.333333333333, -3.888888888889, 1.555555555556, 0, 0, 0,
    0, 0, 1, -1.666666666667, 0.6666666666667, 0,
    2.333333333333, -3.888888888889, 1.555555555556, 0, 0.6666666666667, -1.111111111111,
    0.4444444444444, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0,

    3.857142857143, -1.285714285714, -2.571428571429, 0, -1.285714285714, 0.4285714285714,
    0.8571428571429, 0, -2.571428571429, 0.8571428571429, 1.714285714286, 0,
    0, 0, 0, 0, 9, -3,
    -6, 0, -3, 1, 2, 0,
    -6, 2, 4, 0, 0, 0,
    0, 0, 2.571428571429, -0.8571428571429, -1.714285714286, 0,
    -0.8571428571429, 0.2857142857143, 0.5714285714286, 0, -1.714285714286, 0.5714285714286,
    1.142857142857, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0,

    3.857142857143, -1.285714285714, -2.571428571429, 0, 9, -3,
    -6, 0, 2.571428571429, -0.8571428571429, -1.714285714286, 0,
    0, 0, 0, 0, -1.285714285714, 0.4285714285714,
    0.8571428571429, 0, -3, 1, 2, 0,
    -0.8571428571429, 0.2857142857143, 0.5714285714286, 0, 0, 0,
    0, 0, -2.571428571429, 0.8571428571429, 1.714285714286, 0,
    -6, 2, 4, 0, -1.714285714286, 0.5714285714286,
    1.142857142857, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0,

    3.857142857143, -1.285714285714, -2.571428571429, 0, -1.285714285714, 0.4285714285714,
    0.8571428571429, 0, -2.571428571429, 0.8571428571429, 1.714285714286, 0,
    0, 0, 0, 0, 9, -3,
    -6, 0, -3, 1, 2, 0,
    -6, 2, 4, 0, 0, 0,
    0, 0, 2.571428571429, -0.8571428571429, -1.714285714286, 0,
    -0.8571428571429, 0.2857142857143, 0.5714285714286, 0, -1.714285714286, 0.5714285714286,
    1.142857142857, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0,

    5.037900874636, 11.755102040816, 3.35860058309, 0, -7.556851311953, -17.632653061224,
    -5.037900874636, 0, 2.518950437318, 5.877551020408, 1.679300291545, 0,
    0, 0, 0, 0, 11.755102040816, 27.428571428571,
    7.836734693878, 0, -17.632653061224, -41.142857142857, -11.755102040816, 0,
    5.877551020408, 13.714285714286, 3.918367346939, 0, 0, 0,
    0, 0, 3.35860058309, 7.836734693878, 2.239067055394, 0,
    -5.037900874636, -11.755102040816, -3.35860058309, 0, 1.679300291545, 3.918367346939,
    1.119533527697, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0,

    6.612244897959, 15.428571428571, 4.408163265306, 0, -2.204081632653, -5.142857142857,
    -1.469387755102, 0, -4.408163265306, -10.285714285714, -2.938775510204, 0,
    0, 0, 0, 0, -2.204081632653, -5.142857142857,
    -1.469387755102, 0, 0.734693877551, 1.714285714286, 0.4897959183673, 0,
    1.469387755102, 3.428571428571, 0.9795918367347, 0, 0, 0,
    0, 0, -4.408163265306, -10.285714285714, -2.938775510204, 0,
    1.469387755102, 3.428571428571, 0.9795918367347, 0, 2.938775510204, 6.857142857143,
    1.959183673469, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0,

    3.857142857143, -1.285714285714, -2.571428571429, 0, 9, -3,
    -6, 0, 2.571428571429, -0.8571428571429, -1.714285714286, 0,
    0, 0, 0, 0, -1.285714285714, 0.4285714285714,
    0.8571428571429, 0, -3, 1, 2, 0,
    -0.8571428571429, 0.2857142857143, 0.5714285714286, 0, 0, 0,
    0, 0, -2.571428571429, 0.8571428571429, 1.714285714286, 0,
    -6, 2, 4, 0, -1.714285714286, 0.5714285714286,
    1.142857142857, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0,

    6.612244897959, 15.428571428571, 4.408163265306, 0, -2.204081632653, -5.142857142857,
    -1.469387755102, 0, -4.408163265306, -10.285714285714, -2.938775510204, 0,
    0, 0, 0, 0, -2.204081632653, -5.142857142857,
    -1.469387755102, 0, 0.734693877551, 1.714285714286, 0.4897959183673, 0,
    1.469387755102, 3.428571428571, 0.9795918367347, 0, 0, 0,
    0, 0, -4.408163265306, -10.285714285714, -2.938775510204, 0,
    1.469387755102, 3.428571428571, 0.9795918367347, 0, 2.938775510204, 6.857142857143,
    1.959183673469, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0,

    5.037900874636, 11.755102040816, 3.35860058309, 0, 11.755102040816, 27.428571428571,
    7.836734693878, 0, 3.35860058309, 7.836734693878, 2.239067055394, 0,
    0, 0, 0, 0, -7.556851311953, -17.632653061224,
    -5.037900874636, 0, -17.632653061224, -41.142857142857, -11.755102040816, 0,
    -5.037900874636, -11.755102040816, -3.35860058309, 0, 0, 0,
    0, 0, 2.518950437318, 5.877551020408, 1.679300291545, 0,
    5.877551020408, 13.714285714286, 3.918367346939, 0, 1.679300291545, 3.918367346939,
    1.119533527697, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0
  };

  i = 0;
  EXPECT_EQ(dNdX.rows(), 64U);
  EXPECT_EQ(dNdX.cols(), 3U);
  for (double v : dNdX)
    EXPECT_NEAR(dndx[i++], v, 1.0e-13);

  i = 0;
  EXPECT_EQ(d2NdX2.dim(1), 64U);
  EXPECT_EQ(d2NdX2.dim(2), 3U);
  EXPECT_EQ(d2NdX2.dim(3), 3U);
  for (double v : d2NdX2)
    EXPECT_NEAR(d2ndx2[i++], v, 1.0e-12);

  for (i = 1; i <= H.dim(1); i++)
    for (j = 1; j <= H.dim(2); j++)
      for (k = 1; k <= H.dim(3); k++)
        if (k > 1 && j > 1 && i == j && i == k)
          EXPECT_NEAR(H(i,j,k), 2.0/3.0, 1.0e-13);
        else
          EXPECT_NEAR(H(i,j,k), 0.0, 1.0e-13);
}
