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
#include "SIM1D.h"
#include "ASMs2D.h"
#include "ASMs2Dmx.h"
#include "SIM2D.h"
#include "ASMs3D.h"
#include "SIM3D.h"
#include "SplineUtils.h"

#include "GoTools/geometry/SplineSurface.h"

#include "gtest/gtest.h"
#include <fstream>
#include <array>

static Matrix readMatrix(size_t r, size_t c, const std::string& file)
{
  Matrix result(r,c);
  std::ifstream f(file);
  for (size_t i=1;i<=r;++i)
    for (size_t j=1;j<=c;++j)
      f >> result(i,j);

  return result;
}

static Matrix3D readMatrices(size_t r, size_t c, size_t k, const std::string& file)
{
  Matrix3D result(r,c,k);
  std::ifstream f(file);
  for (size_t n=1;n<=k;++n)
    for (size_t i=1;i<=r;++i)
      for (size_t j=1;j<=c;++j)
        f >> result(i,j,n);

  return result;
}

auto&& CHECK_MATRICES_EQUAL = [](const Matrix& A, const std::string& path)
{
  Matrix B = readMatrix(A.rows(), A.cols(), path);
  for (size_t i=1;i<=A.rows();++i)
    for (size_t j=1;j<=A.cols();++j)
      ASSERT_NEAR(A(i,j), B(i,j), 1e-13);
};


auto&& CHECK_MATRICES3D_EQUAL = [](const Matrix3D& A, const std::string& path)
{
  Matrix3D B = readMatrices(A.dim(1), A.dim(2), A.dim(3), path);
  for (size_t i=1;i<=A.dim(1);++i)
    for (size_t j=1;j<=A.dim(2);++j)
      for (size_t k=1;k<=A.dim(3);++k)
        ASSERT_NEAR(A(i,j,k), B(i,j,k), 1e-13);
};


TEST(TestCoordinateMapping, Jacobian1D)
{
  SIM1D sim(1);
  sim.createDefaultModel();
  ASMs1D& p = static_cast<ASMs1D&>(*sim.getPatch(1));
  p.uniformRefine(3);
  sim.preprocess();

  Vector N;
  utl::matrix<Real> dNdU;
  p.extractBasis(0.25, N, dNdU);
  utl::matrix<Real> X;
  p.getElementCoordinates(X, 2);

  utl::matrix<Real> J;
  utl::matrix<Real> dNdX;
  Real det = utl::Jacobian(J, dNdX, X, dNdU, true);
  ASSERT_FLOAT_EQ(det, 1.0);
  ASSERT_FLOAT_EQ(J(1,1), 1.0);
  CHECK_MATRICES_EQUAL(dNdX, "src/Utility/Test/refdata/Jacobian1D_dNdX.asc");
}


TEST(TestCoordinateMapping, Hessian1D)
{
  SIM1D sim(1);
  sim.createDefaultModel();
  ASMs1D& p = static_cast<ASMs1D&>(*sim.getPatch(1));
  p.raiseOrder(2);
  p.uniformRefine(3);
  sim.preprocess();

  Vector N;
  utl::matrix<Real> dNdU;
  Matrix3D d2Ndu2;
  p.extractBasis(0.25, N, dNdU, d2Ndu2);
  utl::matrix<Real> X;
  p.getElementCoordinates(X, 2);

  utl::matrix<Real> J;
  utl::matrix<Real> dNdX;
  Matrix3D d2NdX2;
  Matrix3D H;
  Real det = utl::Jacobian(J, dNdX, X, dNdU, true);
  utl::Hessian(H, d2NdX2, J, X, d2Ndu2, dNdX);

  ASSERT_FLOAT_EQ(det, 1.0);
  ASSERT_FLOAT_EQ(J(1,1), 1.0);
  CHECK_MATRICES_EQUAL(dNdX, "src/Utility/Test/refdata/Hessian1D_dNdX.asc");
  CHECK_MATRICES3D_EQUAL(d2NdX2, "src/Utility/Test/refdata/Hessian1D_d2NdX2.asc");
  ASSERT_NEAR(H(1,1,1), 0.0, 1e-12);
}


TEST(TestCoordinateMapping, Jacobian2D)
{
  SIM2D sim(1);
  sim.createDefaultModel();
  ASMs2D& p = static_cast<ASMs2D&>(*sim.getPatch(1));
  p.uniformRefine(1,3);
  p.uniformRefine(2,3);
  sim.preprocess();

  Vector N;
  utl::matrix<Real> dNdU;
  p.extractBasis(0.25, 0.25, N, dNdU);
  utl::matrix<Real> X;
  p.getElementCoordinates(X, 2);

  utl::matrix<Real> J;
  utl::matrix<Real> dNdX;
  Real det = utl::Jacobian(J, dNdX, X, dNdU, true);
  ASSERT_FLOAT_EQ(det, 1.0);
  CHECK_MATRICES_EQUAL(J, "src/Utility/Test/refdata/Jacobian2D_J.asc");
  CHECK_MATRICES_EQUAL(dNdX, "src/Utility/Test/refdata/Jacobian2D_dNdX.asc");
}


TEST(TestCoordinateMapping, Hessian2D)
{
  SIM2D sim(1);
  sim.createDefaultModel();
  ASMs2D& p = static_cast<ASMs2D&>(*sim.getPatch(1));
  p.raiseOrder(2,2);
  p.uniformRefine(1,3);
  p.uniformRefine(2,3);
  sim.preprocess();

  Vector N;
  utl::matrix<Real> dNdU;
  Matrix3D d2Ndu2;
  p.extractBasis(0.25, 0.25, N, dNdU, d2Ndu2);
  utl::matrix<Real> X;
  p.getElementCoordinates(X, 2);

  utl::matrix<Real> J;
  utl::matrix<Real> dNdX;
  Matrix3D d2NdX2;
  Matrix3D H;
  Real det = utl::Jacobian(J, dNdX, X, dNdU, true);
  utl::Hessian(H, d2NdX2, J, X, d2Ndu2, dNdX);
  ASSERT_FLOAT_EQ(det, 1.0);
  CHECK_MATRICES_EQUAL(J, "src/Utility/Test/refdata/Hessian2D_J.asc");
  CHECK_MATRICES3D_EQUAL(H, "src/Utility/Test/refdata/Hessian2D_H.asc");
}

TEST(TestCoordinateMapping, Hessian2D_mixed)
{
  SIM2D sim({1,1}, false);
  sim.createDefaultModel();
  ASMs2Dmx& p = static_cast<ASMs2Dmx&>(*sim.getPatch(1));
  sim.preprocess();

  std::vector<std::vector<Go::BasisDerivsSf2>> splinex2(2);
  p.getBasis(1)->computeBasisGrid({0.5}, {0.5}, splinex2[0]);
  p.getBasis(2)->computeBasisGrid({0.5}, {0.5}, splinex2[1]);

  std::array<Vector, 2> N;
  std::array<Matrix, 2> dNxdu;
  std::array<Matrix3D, 2> d2Nxdu2;
  SplineUtils::extractBasis(splinex2[0][0], N[0], dNxdu[0], d2Nxdu2[0]);
  SplineUtils::extractBasis(splinex2[1][0], N[1], dNxdu[1], d2Nxdu2[1]);

  Matrix Jac;
  std::array<Matrix, 2> grad;

  int geoBasis = ASMmxBase::geoBasis;

  Matrix Xnod;
  p.getElementCoordinates(Xnod,1);

  utl::Jacobian(Jac, grad[geoBasis-1], Xnod, dNxdu[geoBasis-1]);

  grad[1-(geoBasis-1)].multiply(dNxdu[1-(geoBasis-1)],Jac);

  Matrix3D Hess;
  std::array<Matrix3D, 2> hess;
  utl::Hessian(Hess, hess[1],Jac,Xnod,d2Nxdu2[1],grad[1],true);
  utl::Hessian(Hess, hess[0],Jac,Xnod,d2Nxdu2[0],grad[0],false);

  // geometry mapping should be the identify mapping
  for (size_t d = 0; d < 2; ++ d)
    for (size_t i=0; i < d2Nxdu2[d].size(); ++i)
      ASSERT_FLOAT_EQ(d2Nxdu2[d].ptr()[i], hess[d].ptr()[i]);

  CHECK_MATRICES3D_EQUAL(hess[0], "src/Utility/Test/refdata/Hessian2D_mixed_b1.asc");
  CHECK_MATRICES3D_EQUAL(hess[1], "src/Utility/Test/refdata/Hessian2D_mixed_b2.asc");
}


TEST(TestCoordinateMapping, Jacobian3D)
{
  SIM3D sim(1);
  sim.createDefaultModel();
  ASMs3D& p = static_cast<ASMs3D&>(*sim.getPatch(1));
  p.uniformRefine(0, 3);
  p.uniformRefine(1, 3);
  p.uniformRefine(2, 3);
  sim.preprocess();

  Vector N;
  utl::matrix<Real> dNdU;
  p.extractBasis(0.25, 0.25, 0.25, N, dNdU);
  utl::matrix<Real> X;
  p.getElementCoordinates(X, 2);

  utl::matrix<Real> J;
  utl::matrix<Real> dNdX;
  Real det = utl::Jacobian(J, dNdX, X, dNdU, true);
  ASSERT_FLOAT_EQ(det, 1.0);
  CHECK_MATRICES_EQUAL(J, "src/Utility/Test/refdata/Jacobian3D_J.asc");
  CHECK_MATRICES_EQUAL(dNdX, "src/Utility/Test/refdata/Jacobian3D_dNdX.asc");
}


TEST(TestCoordinateMapping, Hessian3D)
{
  SIM3D sim(1);
  sim.createDefaultModel();
  ASMs3D& p = static_cast<ASMs3D&>(*sim.getPatch(1));
  p.raiseOrder(2,2,2);
  p.uniformRefine(0,3);
  p.uniformRefine(1,3);
  p.uniformRefine(2,3);
  sim.preprocess();

  Vector N;
  utl::matrix<Real> dNdU;
  Matrix3D d2Ndu2;
  p.extractBasis(0.25, 0.25, 0.25, N, dNdU, d2Ndu2);
  utl::matrix<Real> X;
  p.getElementCoordinates(X, 2);

  utl::matrix<Real> J;
  utl::matrix<Real> dNdX;
  Matrix3D d2NdX2;
  Matrix3D H;
  Real det = utl::Jacobian(J, dNdX, X, dNdU, true);
  utl::Hessian(H, d2NdX2, J, X, d2Ndu2, dNdX);
  ASSERT_FLOAT_EQ(det, 0.34027779);
  CHECK_MATRICES_EQUAL(J, "src/Utility/Test/refdata/Hessian3D_J.asc");
  CHECK_MATRICES_EQUAL(dNdX, "src/Utility/Test/refdata/Hessian3D_dNdX.asc");
  CHECK_MATRICES3D_EQUAL(d2NdX2, "src/Utility/Test/refdata/Hessian3D_d2NdX2.asc");
  CHECK_MATRICES3D_EQUAL(H, "src/Utility/Test/refdata/Hessian3D_H.asc");
}
