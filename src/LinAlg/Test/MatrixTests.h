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

#include "matrixnd.h"

#include "Catch2Support.h"

#include <limits>
#include <numeric>


namespace {

template<class Scalar>
void vectorAddTest()
{
  utl::vector<Scalar> d(12);
  std::iota(d.begin(), d.end(), 1.0);

  utl::vector<Scalar> d2(d);
  d2.add(d, 1.5);
  std::cout <<"d2:"<< d2 << std::endl;

  utl::vector<Scalar> d3(d.size() / 2);
  d3.add(d, 0.5, 1, 2);
  std::cout <<"d3:"<< d3 << std::endl;

  utl::vector<Scalar> d4(d.size());
  d4.add(d3, 1.0, 0, 1, 0, 2);
  d4.add(d3, 1.0, 0, 1, 1, 2);
  std::cout <<"d4:"<< d4 << std::endl;

  for (size_t i = 0; i < d2.size(); ++i) {
    REQUIRE_THAT(d2[i], WithinRel(2.5*(i+1)));
    REQUIRE_THAT(d4[i], WithinRel(d3[i/2]));
  }

  for (size_t i = 0; i < d3.size(); ++i)
    REQUIRE_THAT(d3[i], WithinRel(i + 1.0));

  d2.fill(0.0);
  d2.add(d, 1.0, 0, 3, 0, 3);
  d2.add(d, 1.0, 1, 3, 1, 3);
  d2.add(d, 1.0, 2, 3, 2, 3);
  std::cout <<"d5:"<< d2 << std::endl;
  for (size_t i = 0; i < d.size(); ++i)
    REQUIRE_THAT(d[i], WithinRel(d2[i]));
  d2.add(d, 1.0, 0, 1, 7, 3);
  std::cout <<"d6:"<< d2 << std::endl;
  for (size_t i = 7; i < d.size(); i += 3)
    REQUIRE_THAT(d2[i], WithinRel(d[i]+(i-4)/3));
}


template<class Scalar>
void vectorDotTest()
{
  constexpr size_t size = 10;
  constexpr Scalar max = Scalar((size-1)*size/2);

  utl::vector<Scalar> d(size);
  for (size_t i = 0; i < size; ++i)
    d[i] = std::sqrt(i);

  REQUIRE_THAT(d.dot(d), WithinRel(max));
}


template<class Scalar>
void vectorMultiplyTest()
{
  constexpr size_t size = 10;

  utl::vector<Scalar> d(size);
  std::iota(d.begin(), d.end(), 0.0);
  d *= 2.0;

  for (size_t i = 0; i < size; ++i)
    REQUIRE_THAT(d[i], WithinRel(2.0 * i));
}


template<class Scalar>
void vectorNormTest()
{
  constexpr size_t size = 10;
  const Scalar max = sqrt(size*(size+1)*(2*size+1) / 6);

  utl::vector<Scalar> d(size);
  std::iota(d.begin(), d.end(), 1.0);

  REQUIRE_THAT(d.normInf(), WithinRel(static_cast<Scalar>(size)));
  REQUIRE_THAT(d.norm2(), WithinRel(max));
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

  REQUIRE(A.multiply(u,x,1.0,0.0,false,3,4,1,2));
  REQUIRE(A.multiply(v,y,1.0,0.0,true,4,2));

  REQUIRE_THAT(x(3),  WithinRel(370.0));
  REQUIRE_THAT(x(7),  WithinRel(410.0));
  REQUIRE_THAT(x(11), WithinRel(450.0));
  REQUIRE_THAT(y(1),  WithinRel(38.0));
  REQUIRE_THAT(y(3),  WithinRel(83.0));
  REQUIRE_THAT(y(5),  WithinRel(128.0));
  REQUIRE_THAT(y(7),  WithinRel(173.0));
  REQUIRE_THAT(y(9),  WithinRel(218.0));

  REQUIRE(A.multiply(u,x,1.0,-1.0,false,3,4,1,2));
  REQUIRE(A.multiply(v,y,1.0,-1.0,true,4,2));

  REQUIRE_THAT(x.sum(), WithinAbs(0.0, 1e-14));
  REQUIRE_THAT(y.sum(), WithinAbs(0.0, 1e-14));

  u.resize(5,utl::RETAIN);
  REQUIRE(A.multiply(u,v));
  v *= 0.5;

  REQUIRE_THAT(v(1), WithinRel(67.5));
  REQUIRE_THAT(v(2), WithinRel(75.0));
  REQUIRE_THAT(v(3), WithinRel(82.5));

  B.multiplyMat(A, static_cast<const std::vector<Scalar>&>(I));
  B2.multiplyMat(A, static_cast<const std::vector<Scalar>&>(I), false, true);
  B3.multiplyMat(A, static_cast<const std::vector<Scalar>&>(I2), true, false);
  REQUIRE(B.rows() == A.rows());
  REQUIRE(B.cols() == A.cols());
  for (size_t i = 1; i <= 3; ++i)
    for (size_t j = 1; j <= 5; ++j) {
      REQUIRE_THAT(B(i,j), WithinRel(A(i,j)));
      REQUIRE_THAT(B2(i,j), WithinRel(2.0*A(i,j)));
      REQUIRE_THAT(B3(j,i), WithinRel(A(i,j)));
    }

  REQUIRE(A.multiply(u,v));
  REQUIRE_THAT(v(1), WithinRel(135.0));
  REQUIRE_THAT(v(2), WithinRel(150.0));
  REQUIRE_THAT(v(3), WithinRel(165.0));
}


template<class Scalar>
void normTest()
{
  constexpr Scalar eps = std::numeric_limits<Scalar>::epsilon()*10;

  utl::matrix<Scalar> a(4,5);
  std::iota(a.begin(),a.end(),1.0);
  std::cout <<"A:"<< a;

  REQUIRE_THAT(a.sum(),     WithinRel(210.0));
  REQUIRE_THAT(a.sum(5),    WithinRel(34.0));
  REQUIRE_THAT(a.asum(5),   WithinRel(34.0));
  REQUIRE_THAT(a.trace(),   WithinRel(34.0));
  REQUIRE_THAT(a.norm2(5),  WithinAbs(sqrt(414.0), eps));
  REQUIRE_THAT(a.normInf(), WithinRel(60.0));
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

  REQUIRE(A.rows() == size1);
  REQUIRE(A.cols() == size2);

  for (size_t i = 1; i <= size1; ++i)
    for (size_t j = 1; j <= size2; ++j) {
      REQUIRE_THAT(A(i,j), WithinRel(u(i)*v(j)));
      REQUIRE_THAT(A2(i,j), WithinRel(3.0*u(i)*v(j)));
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
  REQUIRE(D.multiplyMat(a,B));

  typename std::vector<Scalar>::const_iterator c = C.begin();
  for (const Scalar d : D)
    REQUIRE_THAT(d, WithinRel(*(c++)));
}

}
