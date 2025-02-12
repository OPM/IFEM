//==============================================================================
//!
//! \file TestMatVec.C
//!
//! \date Jan 14 2025
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Tests for some matrix-vector operations.
//!
//==============================================================================

#include "MatVec.h"

#include <gtest/gtest.h>

#include <numeric>


TEST(TestMatVec, add)
{
  Vector d(10);
  std::iota(d.begin(),d.end(),Real(0));

  Vector v = d + d;
  for (size_t i = 0; i < v.size(); i++)
    EXPECT_FLOAT_EQ(v[i], 2*Real(i));
}


TEST(TestMatVec, multiply)
{
  Matrix A(3,5);
  Vector u(5);

  std::iota(A.begin(),A.end(),Real(1));
  std::iota(u.begin(),u.end(),Real(1));

  Vector v = A * u;
  EXPECT_FLOAT_EQ(v(1),Real(135));
  EXPECT_FLOAT_EQ(v(2),Real(150));
  EXPECT_FLOAT_EQ(v(3),Real(165));
}


TEST(TestMatVec, transform)
{
  Matrix A(12,12);
  Vector v(12);
  for (size_t i = 1; i <= A.rows(); i++)
  {
    A(i,i) = v(i) = i;
    for (size_t j = i+1; j <= A.cols(); j++)
      A(i,j) = A(j,i) = i + A.rows()*(j-i);
  }

  const Real alpha = Real(0.5); // approx. 29 degrees

  Matrix T(3,3), B(A);
  Vector w(v);
  T(1,1) =  Real(1);
  T(2,2) =  T(3,3) = cos(alpha);
  T(2,3) =  sin(alpha);
  T(3,2) = -T(2,3);
  std::cout <<"T:"<< T <<"A:"<< A <<"v:"<< v;
  ASSERT_TRUE(utl::transform(A,T));
  std::cout <<"A transformed:"<< A;
  ASSERT_TRUE(utl::transform(v,T));
  std::cout <<"v transformed:"<< v;

  Matrix P(12,12);
  for (size_t k = 0; k < A.rows(); k += 3)
    for (size_t i = 1; i <= 3; i++)
      for (size_t j = 1; j <= 3; j++)
        P(k+i,k+j) = T(i,j);
  std::cout <<"P:"<< P;

  Matrix C;
  Vector u;
  B.multiply(C.multiply(P,B),P,false,true);
  std::cout <<"B:"<< B;
  ASSERT_TRUE(P.multiply(w,u));
  std::cout <<"u:"<< u;

  for (size_t i = 1; i <= A.rows(); i++)
    for (size_t j = 1; j <= A.cols(); j++)
      EXPECT_FLOAT_EQ(A(i,j), B(i,j));

  for (size_t i = 1; i <= v.size(); i++)
    EXPECT_FLOAT_EQ(v(i), u(i));
}
