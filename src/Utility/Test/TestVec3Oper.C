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

#include "gtest/gtest.h"

TEST(TestVec3Oper, MxV)
{
  utl::matrix<Real> A(3,3);
  A.fill(1.0);
  std::vector<Real> x(3, 1.0);
  Vec3 result = A*x;
  EXPECT_TRUE(result[0] == 3.0);
  EXPECT_TRUE(result[1] == 3.0);
  EXPECT_TRUE(result[2] == 3.0);
}

TEST(TestVec3Oper, MultiplyScalar)
{
  Vec3 a;
  a[0] = a[1] = a[2] = 1.0;
  Vec3 result = a*10.0;

  Vec3 result2 = 10.0*a;

  EXPECT_TRUE(result == result2);
  EXPECT_TRUE(result[0] == 10.0);
  EXPECT_TRUE(result[1] == 10.0);
  EXPECT_TRUE(result[2] == 10.0);
}

TEST(TestVec3Oper, DivideScalar)
{
  Vec3 a;
  a[0] = a[1] = a[2] = 10.0;
  Vec3 result = a/10.0;

  EXPECT_TRUE(result[0] == 1.0);
  EXPECT_TRUE(result[1] == 1.0);
  EXPECT_TRUE(result[2] == 1.0);
}

TEST(TestVec3Oper, Dot)
{
  Vec3 a;
  a[0] = a[1] = a[2] = 1.0;
  Vec3 b;
  b[0] = b[1] = b[2] = 1.0;

  Real result = a*b;
  EXPECT_TRUE(result == 3.0);
}

TEST(TestVec3Oper, Addition)
{
  Vec3 a;
  a[0] = a[1] = a[2] = 1.0;
  Vec3 b;
  b[0] = b[1] = b[2] = 1.0;

  Vec3 result = a+b;
  EXPECT_TRUE(result[0] == 2.0);
  EXPECT_TRUE(result[1] == 2.0);
  EXPECT_TRUE(result[2] == 2.0);
}

TEST(TestVec3Oper, Subtraction)
{
  Vec3 a;
  a[0] = a[1] = a[2] = 1.0;
  Vec3 b;
  b[0] = b[1] = b[2] = 1.0;

  Vec3 result = a-b;
  EXPECT_TRUE(result[0] == 0.0);
  EXPECT_TRUE(result[1] == 0.0);
  EXPECT_TRUE(result[2] == 0.0);
}

TEST(TestVec3Oper, Equality)
{
  Vec3 a;
  a[0] = a[1] = a[2] = 1.0;
  Vec3 b;
  b[0] = b[1] = b[2] = 1.0;

  EXPECT_TRUE(a == b);
}

TEST(TestVec3Oper, InEquality)
{
  Vec3 a;
  a[0] = a[1] = a[2] = 1.0;
  Vec3 b;
  b[0] = b[1] = b[2] = 2.0;

  EXPECT_TRUE(a != b);
}

TEST(TestVec3Oper, Less)
{
  Vec3 a;
  a[0] = a[1] = a[2] = 1.0;
  Vec3 b;
  b[0] = b[1] = b[2] = 2.0;

  EXPECT_TRUE(a < b);
}

TEST(TestVec3Oper, StreamOut)
{
  Vec3 a;
  a[0] = a[1] = a[2] = 1.0;
  std::stringstream str;

  str << a;

  EXPECT_STREQ(str.str().c_str(), "1 1 1");
}

TEST(TestVec3Oper, StreamIn)
{
  std::stringstream str;
  str << "1.0 1.0 1.0";
  Vec3 result;

  str >> result;

  EXPECT_TRUE(result[0] == 1.0);
  EXPECT_TRUE(result[1] == 1.0);
  EXPECT_TRUE(result[2] == 1.0);
}
