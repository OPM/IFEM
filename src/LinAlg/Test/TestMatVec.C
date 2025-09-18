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

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <numeric>

using Catch::Matchers::WithinRel;


TEST_CASE("TestMatVec.Add")
{
  Vector d(10);
  std::iota(d.begin(),d.end(),Real(0));

  Vector v = d + d;
  for (size_t i = 0; i < v.size(); i++)
    REQUIRE_THAT(v[i], WithinRel(2*Real(i)));
}


TEST_CASE("TestMatVec.Multiply")
{
  Matrix A(3,5);
  Vector u(5);

  std::iota(A.begin(),A.end(),Real(1));
  std::iota(u.begin(),u.end(),Real(1));

  Vector v = A * u;
  REQUIRE_THAT(v(1), WithinRel(Real(135)));
  REQUIRE_THAT(v(2), WithinRel(Real(150)));
  REQUIRE_THAT(v(3), WithinRel(Real(165)));
}


TEST_CASE("TestMatVec.Transform")
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
  REQUIRE(utl::transform(A,T));
  std::cout <<"A transformed:"<< A;
  REQUIRE(utl::transform(v,T));
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
  REQUIRE(P.multiply(w,u));
  std::cout <<"u:"<< u;

  for (size_t i = 1; i <= A.rows(); i++)
    for (size_t j = 1; j <= A.cols(); j++)
      REQUIRE_THAT(A(i,j), WithinRel(B(i,j)));

  for (size_t i = 1; i <= v.size(); i++)
    REQUIRE_THAT(v(i), WithinRel(u(i)));
}
