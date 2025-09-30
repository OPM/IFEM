//==============================================================================
//!
//! \file TestVec3Oper.C
//!
//! \date Oct 7 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for global algebraic operations involving Vec3 objects.
//!
//==============================================================================

#include "Vec3.h"
#include "Vec3Oper.h"

#include "Catch2Support.h"

#include <sstream>


TEST_CASE("TestVec3Oper.GetAndSet")
{
  Vec3 a(1.0,2.0,3.0);
  REQUIRE_THAT(a[0], WithinRel(1.0));
  REQUIRE_THAT(a[1], WithinRel(2.0));
  REQUIRE_THAT(a[2], WithinRel(3.0));
  REQUIRE_THAT(a(1), WithinRel(1.0));
  REQUIRE_THAT(a(2), WithinRel(2.0));
  REQUIRE_THAT(a(3), WithinRel(3.0));
  REQUIRE_THAT(a.x,  WithinRel(1.0));
  REQUIRE_THAT(a.y,  WithinRel(2.0));
  REQUIRE_THAT(a.z,  WithinRel(3.0));
  a[0] = 4.0;
  a(2) = 5.0;
  a.z  = 6.0;
  REQUIRE_THAT(a[0], WithinRel(4.0));
  REQUIRE_THAT(a[1], WithinRel(5.0));
  REQUIRE_THAT(a[2], WithinRel(6.0));
}


TEST_CASE("TestVec3Oper.MxV")
{
  const double data[9] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
  utl::matrix<Real> A(3,3); A.fill(data);
  std::vector<Real> x(data, data+3);
  Vec3 result = A * x;

  REQUIRE_THAT(result.x, WithinRel(30.0));
  REQUIRE_THAT(result.y, WithinRel(36.0));
  REQUIRE_THAT(result.z, WithinRel(42.0));
}


TEST_CASE("TestVec3Oper.MultiplyScalar")
{
  Vec3 a(1.0,2.0,3.0);
  Vec3 result = a*10.0;
  Vec3 result2 = 10.0*a;

  REQUIRE(result == result2);
  REQUIRE_THAT(result.x, WithinRel(10.0));
  REQUIRE_THAT(result.y, WithinRel(20.0));
  REQUIRE_THAT(result.z, WithinRel(30.0));
}


TEST_CASE("TestVec3Oper.DivideScalar")
{
  Vec3 a(1.0,2.0,3.0);
  Vec3 result = a/10.0;

  REQUIRE_THAT(result.x, WithinRel(0.1));
  REQUIRE_THAT(result.y, WithinRel(0.2));
  REQUIRE_THAT(result.z, WithinRel(0.3));
}


TEST_CASE("TestVec3Oper.Dot")
{
  Vec3 a(1.0,2.0,3.0);
  Vec3 b(4.0,5.0,6.0);
  Real result = a * b;

  REQUIRE_THAT(result, WithinRel(32.0));
}


TEST_CASE("TestVec3Oper.Addition")
{
  Vec3 a(1.0,2.0,3.0);
  Vec3 b(4.0,5.0,6.0);
  Vec3 result = a + b;

  REQUIRE_THAT(result.x, WithinRel(5.0));
  REQUIRE_THAT(result.y, WithinRel(7.0));
  REQUIRE_THAT(result.z, WithinRel(9.0));
}


TEST_CASE("TestVec3Oper.Subtraction")
{
  Vec3 a(1.0,2.0,3.0);
  Vec3 b(4.0,5.0,6.0);
  Vec3 result = a - b;

  REQUIRE_THAT(result.x, WithinRel(-3.0));
  REQUIRE_THAT(result.y, WithinRel(-3.0));
  REQUIRE_THAT(result.z, WithinRel(-3.0));
}


TEST_CASE("TestVec3Oper.Length")
{
  Vec3 a(1.0,2.0,3.0);
  double alen = sqrt(14.0);

  REQUIRE_THAT(a.length(), WithinRel(alen));
}


TEST_CASE("TestVec3Oper.Equality")
{
  Vec3 a(1.0,2.0,3.0);
  Vec3 b(4.0,5.0,6.0);
  Vec3 c(1.0,2.0,3.0);

  REQUIRE(a == c);
  REQUIRE(!(a == b));
}


TEST_CASE("TestVec3Oper.InEquality")
{
  Vec3 a(1.0,2.0,3.0);
  Vec3 b(4.0,5.0,6.0);
  Vec3 c(1.0,2.0,3.0);

  REQUIRE(a != b);
  REQUIRE(!(a != c));
}


TEST_CASE("TestVec3Oper.Less")
{
  Vec3 a(1.0,2.0,3.0);
  Vec3 b(4.0,5.0,6.0);

  REQUIRE(a < b);
  REQUIRE(!(b < a));
}


TEST_CASE("TestVec3Oper.StreamOut")
{
  Vec3 a(1.0,2.0,3.0);
  std::stringstream str;
  str << a;

  REQUIRE(str.str() == "1 2 3");
}


TEST_CASE("TestVec3Oper.StreamIn")
{
  std::stringstream str;
  str << "1.0 2.0 3.0";
  Vec3 result;
  str >> result;

  REQUIRE_THAT(result.x, WithinRel(1.0));
  REQUIRE_THAT(result.y, WithinRel(2.0));
  REQUIRE_THAT(result.z, WithinRel(3.0));
}
