//==============================================================================
//!
//! \file TestLegendre.C
//!
//! \date Oct 10 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for various utility methods for Spectral elements.
//!
//==============================================================================

#include "Legendre.h"
#include <fstream>

#include "gtest/gtest.h"

static Matrix readMatrix(size_t r, size_t c, const std::string& file)
{
  Matrix result(r,c);
  std::ifstream f(file);
  for (size_t i=1;i<=r;++i)
    for (size_t j=1;j<=c;++j)
      f >> result(i,j);

  return result;
}

#define CHECK_MATRICES_EQUAL(A,path) \
  do { \
  Matrix B = readMatrix(A.rows(), A.cols(), path); \
  for (size_t i=1;i<=A.rows();++i) \
    for (size_t j=1;j<=A.cols();++j) \
      ASSERT_NEAR(A(i,j), B(i,j), 1e-13); \
  } while(0);

TEST(TestLegendre, GL)
{
  RealArray w, p;
  ASSERT_TRUE(Legendre::GL(w, p, 3));
  ASSERT_FLOAT_EQ(w[0], 0.5555555555555561);
  ASSERT_FLOAT_EQ(w[1], 0.8888888888888888);
  ASSERT_FLOAT_EQ(w[2], 0.5555555555555561);
  ASSERT_FLOAT_EQ(p[0], -0.7745966692414832);
  ASSERT_NEAR(p[1], 0.0, 1e-14);
  ASSERT_FLOAT_EQ(p[2], 0.7745966692414832);

  ASSERT_TRUE(Legendre::GL(w, p, 4));
  ASSERT_FLOAT_EQ(w[0], 0.3478548451374518);
  ASSERT_FLOAT_EQ(w[1], 0.6521451548625461);
  ASSERT_FLOAT_EQ(w[2], 0.6521451548625461);
  ASSERT_FLOAT_EQ(w[3], 0.3478548451374518);
  ASSERT_FLOAT_EQ(p[0],-0.8611363115940535);
  ASSERT_FLOAT_EQ(p[1],-0.3399810435848561);
  ASSERT_FLOAT_EQ(p[2], 0.3399810435848561);
  ASSERT_FLOAT_EQ(p[3], 0.8611363115940535);
}

TEST(TestLegendre, GLL)
{
  Vector w, p;
  ASSERT_TRUE(Legendre::GLL(w, p, 3));
  ASSERT_FLOAT_EQ(w[0], 0.3333333333333333);
  ASSERT_FLOAT_EQ(w[1], 1.333333333333333);
  ASSERT_FLOAT_EQ(w[2], 0.3333333333333333);
  ASSERT_FLOAT_EQ(p[0], -1.0);
  ASSERT_FLOAT_EQ(p[1], 0.0);
  ASSERT_FLOAT_EQ(p[2], 1.0);

  ASSERT_TRUE(Legendre::GLL(w, p, 4));
  ASSERT_FLOAT_EQ(w[0], 0.1666666666666667);
  ASSERT_FLOAT_EQ(w[1], 0.8333333333333335);
  ASSERT_FLOAT_EQ(w[2], 0.8333333333333335);
  ASSERT_FLOAT_EQ(w[3], 0.1666666666666667);
  ASSERT_FLOAT_EQ(p[0],-1.0);
  ASSERT_FLOAT_EQ(p[1],-0.447213595499958);
  ASSERT_FLOAT_EQ(p[2], 0.447213595499958);
  ASSERT_FLOAT_EQ(p[3], 1.0);
}

TEST(TestLegendre, Eval)
{
  Real result1, result2;

  Legendre::LegendreEval(3, 0.3, result1);
  Legendre::LegendreEval(5, 0.7, result2);
  ASSERT_FLOAT_EQ(result1, -0.3825);
  ASSERT_FLOAT_EQ(result2, -0.3651987499999999);
}

TEST(TestLegendre, Derivative)
{
  Real result1, result2;

  Legendre::LegendreDerEval(3, 0.3, result1);
  Legendre::LegendreDerEval(5, 0.7, result2);
  ASSERT_FLOAT_EQ(result1, -0.8250000000000001);
  ASSERT_FLOAT_EQ(result2, -1.5335625);
}

TEST(TestLegendre, BasisDerivatives)
{
  Matrix result1, result2;

  Legendre::basisDerivatives(3, result1);
  CHECK_MATRICES_EQUAL(result1, "src/Utility/Test/refdata/Legendre_BasisDerivatives_3.asc");
  Legendre::basisDerivatives(10, result2);
  CHECK_MATRICES_EQUAL(result2, "src/Utility/Test/refdata/Legendre_BasisDerivatives_10.asc");
}
