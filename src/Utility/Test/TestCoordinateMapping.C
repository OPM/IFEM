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
#include "SIM2D.h"
#include "ASMs3D.h"
#include "SIM3D.h"

#include "gtest/gtest.h"
#include <fstream>

Matrix readMatrix(size_t r, size_t c, const std::string& file)
{
  Matrix result(r,c);
  std::ifstream f(file);
  for (size_t i=1;i<=r;++i)
    for (size_t j=1;j<=c;++j)
      f >> result(i,j);

  return result;
}

Matrix3D readMatrices(size_t r, size_t c, size_t k, const std::string& file)
{
  Matrix3D result(r,c,k);
  std::ifstream f(file);
  for (size_t n=1;n<=k;++n)
    for (size_t i=1;i<=r;++i)
      for (size_t j=1;j<=c;++j)
        f >> result(i,j,n);

  return result;
}

#define CHECK_MATRICES_EQUAL(A,path) \
  do { \
  Matrix B = readMatrix(A.rows(), A.cols(), path); \
  for (size_t i=1;i<=A.rows();++i) \
    for (size_t j=1;j<=A.cols();++j) \
      ASSERT_NEAR(A(i,j), B(i,j), 1e-13); \
  } while(0);

#define CHECK_MATRICES3D_EQUAL(A,path) \
  do { \
    Matrix3D B = readMatrices(A.dim(1), A.dim(2), A.dim(3), path); \
    for (size_t i=1;i<=A.dim(1);++i) \
      for (size_t j=1;j<=A.dim(2);++j) \
        for (size_t k=1;k<=A.dim(3);++k) \
          ASSERT_NEAR(A(i,j,k), B(i,j,k), 1e-13); \
  } while(0);

#define printMatrix(A) \
  std::cout << "matrix is of size " << A.rows() << "x" << A.cols() << std::endl; \
  for (size_t i=1;i<=A.rows();++i) { \
    for (size_t j=1;j<=A.cols();++j) \
      std::cout << A(i,j) << " ";    \
    std::cout << std::endl; \
  }


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
  utl::zero_print_tol = 1e-16;
  std::cout.precision(17);
  std::cout << H << std::endl;
  CHECK_MATRICES3D_EQUAL(d2NdX2, "src/Utility/Test/refdata/Hessian3D_d2NdX2.asc");
  CHECK_MATRICES3D_EQUAL(H, "src/Utility/Test/refdata/Hessian3D_H.asc");
}
