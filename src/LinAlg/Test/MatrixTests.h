//==============================================================================
//!
//! \file MatrixTests.h
//!
//! \date Apr 11 2016
//!
//! \author Eivind Fonn / SINTEF
//!
//! \brief Implementation of some unit tests for matrix.
//!
//==============================================================================

#include "MatVec.h"

#include "gtest/gtest.h"

#include <limits>
#include <numeric>

namespace {

template<class Scalar>
void vectorAddTest()
{
  constexpr size_t size = 10;

  utl::vector<Scalar> d(size);
  std::iota(d.begin(), d.end(), 0.0);

  utl::vector<Scalar> d2(d);
  d2.add(d, 1.5);

  utl::vector<Scalar> d3(size / 2);
  d3.add(d, 0.5, 1, 2);

  utl::vector<Scalar> d4(size);
  d4.add(d3, 1.0, 0, 1, 0, 2);
  d4.add(d3, 1.0, 0, 1, 1, 2);
  std::cout << d4 << std::endl;

  for (size_t i = 0; i < size; ++i) {
    EXPECT_FLOAT_EQ(d2[i], 2.5*i);
    EXPECT_FLOAT_EQ(d4[i], d3[i/2]);
  }
  for (size_t i = 0; i < size / 2; ++i)
    EXPECT_FLOAT_EQ(d3[i], i + 0.5);

  if constexpr (std::is_same_v<Scalar,Real>) {
    utl::vector<Scalar> v = d + d;
    for (size_t i = 0; i < size; ++i)
      EXPECT_FLOAT_EQ(v[i], 2.0*i);
  }
}


template<class Scalar>
void vectorDotTest()
{
  constexpr size_t size = 10;
  constexpr Scalar max = Scalar(size-1) * Scalar(size) / 2.0;

  utl::vector<Scalar> d(size);
  for (size_t i = 0; i < size; ++i)
    d[i] = std::sqrt(i);

  EXPECT_FLOAT_EQ(d.dot(d), max);
}


template<class Scalar>
void vectorMultiplyTest()
{
  constexpr size_t size = 10;

  utl::vector<Scalar> d(size);
  std::iota(d.begin(), d.end(), 0.0);
  d *= 2.0;

  for (size_t i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(d[i], 2.0 * i);
}


template<class Scalar>
void vectorNormTest()
{
  constexpr size_t size = 10;
  constexpr Scalar max = Scalar(size-1) * Scalar(size) / 2.0;

  utl::vector<Scalar> d(size);
  std::iota(d.begin(), d.end(), 0.0);

  EXPECT_FLOAT_EQ(d.normInf(), size - 1.0);
  EXPECT_FLOAT_EQ(d.norm2(), max);
}


template<class Scalar>
void multiplyTest()
{
  utl::vector<Scalar> u(14), v(9), x, y;
  utl::matrix<Scalar> A(3,5), B, B2(3,5), B3;
  utl::matrix<Scalar> I(5,5), I2(3,3);

  std::iota(u.begin(),u.end(),1.0);
  std::iota(v.begin(),v.end(),1.0);
  std::iota(A.begin(),A.end(),1.0);
  std::iota(B2.begin(),B2.end(),1.0);
  for (size_t i = 1; i <= 5; ++i)
    I(i,i) = 1.0;
  for (size_t i = 1; i <= 3; ++i)
    I2(i,i) = 1.0;

  ASSERT_TRUE(A.multiply(u,x,1.0,0.0,false,3,4,1,2));
  ASSERT_TRUE(A.multiply(v,y,1.0,0.0,true,4,2));

  EXPECT_FLOAT_EQ(x(3),370.0);
  EXPECT_FLOAT_EQ(x(7),410.0);
  EXPECT_FLOAT_EQ(x(11),450.0);
  EXPECT_FLOAT_EQ(y(1),38.0);
  EXPECT_FLOAT_EQ(y(3),83.0);
  EXPECT_FLOAT_EQ(y(5),128.0);
  EXPECT_FLOAT_EQ(y(7),173.0);
  EXPECT_FLOAT_EQ(y(9),218.0);

  ASSERT_TRUE(A.multiply(u,x,1.0,-1.0,false,3,4,1,2));
  ASSERT_TRUE(A.multiply(v,y,1.0,-1.0,true,4,2));

  EXPECT_FLOAT_EQ(x.sum(),0.0);
  EXPECT_FLOAT_EQ(y.sum(),0.0);

  u.std::template vector<Scalar>::resize(5);
  ASSERT_TRUE(A.multiply(u,v));
  v *= 0.5;

  EXPECT_FLOAT_EQ(v(1),67.5);
  EXPECT_FLOAT_EQ(v(2),75.0);
  EXPECT_FLOAT_EQ(v(3),82.5);

  B.multiplyMat(A, static_cast<const std::vector<Scalar>&>(I));
  B2.multiplyMat(A, static_cast<const std::vector<Scalar>&>(I), false, true);
  B3.multiplyMat(A, static_cast<const std::vector<Scalar>&>(I2), true, false);
  EXPECT_EQ(B.rows(), A.rows());
  EXPECT_EQ(B.cols(), A.cols());
  for (size_t i = 1; i <= 3; ++i)
    for (size_t j = 1; j <= 5; ++j) {
      EXPECT_FLOAT_EQ(B(i,j), A(i,j));
      EXPECT_FLOAT_EQ(B2(i,j), 2.0*A(i,j));
      EXPECT_FLOAT_EQ(B3(j,i), A(i,j));
    }

  if constexpr (std::is_same_v<Scalar,Real>) {
    utl::vector<Scalar> v2 = A * u;
    EXPECT_FLOAT_EQ(v2(1),135.0);
    EXPECT_FLOAT_EQ(v2(2),150.0);
    EXPECT_FLOAT_EQ(v2(3),165.0);
  }
}


template<class Scalar>
void normTest()
{
  utl::matrix<Scalar> a(4,5);
  std::iota(a.begin(),a.end(),1.0);
  std::cout <<"A:"<< a;

  EXPECT_FLOAT_EQ(a.sum(),210.0);
  EXPECT_FLOAT_EQ(a.sum(5),34.0);
  EXPECT_FLOAT_EQ(a.asum(5),34.0);
  EXPECT_FLOAT_EQ(a.trace(),34.0);
  EXPECT_NEAR(a.norm2(5),sqrt(414.0),std::numeric_limits<Scalar>::epsilon()*10.0);
  EXPECT_FLOAT_EQ(a.normInf(),60.0);
}


template<class Scalar>
void outerProductTest()
{
  constexpr size_t size1 = 3;
  constexpr size_t size2 = 4;

  utl::vector<Scalar> u(size1);
  utl::vector<Scalar> v(size2);
  std::iota(u.begin(), u.end(), 1.0);
  std::iota(v.begin(), v.end(), 2.0);

  utl::matrix<Scalar> A;
  A.outer_product(u, v);

  utl::matrix<Scalar> A2(A);
  A2.outer_product(u, v, true, 2.0);

  EXPECT_EQ(A.rows(), size1);
  EXPECT_EQ(A.cols(), size2);

  for (size_t i = 1; i <= size1; ++i)
    for (size_t j = 1; j <= size2; ++j) {
      EXPECT_FLOAT_EQ(A(i,j), u(i)*v(j));
      EXPECT_FLOAT_EQ(A2(i,j), 3.0*u(i)*v(j));
    }
}


template<class Scalar>
void matrix3DMultiplyTest()
{
  std::vector<Scalar> a(10);
  utl::matrix<Scalar> A(2,5);
  utl::matrix3d<Scalar> B(5,4,3), C, D;

  std::iota(a.begin(),a.end(),1.0);
  std::iota(A.begin(),A.end(),1.0);
  std::iota(B.begin(),B.end(),1.0);

  C.multiply(A,B);
  ASSERT_TRUE(D.multiplyMat(a,B));

  typename std::vector<Scalar>::const_iterator c = C.begin();
  for (const Scalar d : D)
    EXPECT_FLOAT_EQ(d,*(c++));
}

}
