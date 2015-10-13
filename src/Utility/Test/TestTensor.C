//==============================================================================
//!
//! \file TestTensor.C
//!
//! \date Oct 14 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for second-order tensors with some basic operations.
//!
//==============================================================================

#include "Tensor.h"
#include "Vec3.h"

#include "gtest/gtest.h"

TEST(TestTensor, Shift)
{
  const double data[9] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };

  Tensor T2(std::vector<double>(data,data+4));
  Tensor T3(std::vector<double>(data,data+9));

  Tensor t2(T2);
  t2.shift(1);
  ASSERT_FLOAT_EQ(t2(1,1), 3.0);
  ASSERT_FLOAT_EQ(t2(1,2), 1.0);
  ASSERT_FLOAT_EQ(t2(2,1), 4.0);
  ASSERT_FLOAT_EQ(t2(2,2), 2.0);
  t2 = T2;
  t2.shift(-1);
  ASSERT_FLOAT_EQ(t2(1,1), 3.0);
  ASSERT_FLOAT_EQ(t2(1,2), 1.0);
  ASSERT_FLOAT_EQ(t2(2,1), 4.0);
  ASSERT_FLOAT_EQ(t2(2,2), 2.0);
  t2 = T2;
  t2.shift(2);
  ASSERT_FLOAT_EQ(t2(1,1), 1.0);
  ASSERT_FLOAT_EQ(t2(1,2), 3.0);
  ASSERT_FLOAT_EQ(t2(2,1), 2.0);
  ASSERT_FLOAT_EQ(t2(2,2), 4.0);
  t2 = T2;
  t2.shift(-2);
  ASSERT_FLOAT_EQ(t2(1,1), 1.0);
  ASSERT_FLOAT_EQ(t2(1,2), 3.0);
  ASSERT_FLOAT_EQ(t2(2,1), 2.0);
  ASSERT_FLOAT_EQ(t2(2,2), 4.0);

  Tensor t3(T3);
  t3.shift(1);
  ASSERT_FLOAT_EQ(t3(1,1), 7.0);
  ASSERT_FLOAT_EQ(t3(1,2), 1.0);
  ASSERT_FLOAT_EQ(t3(1,3), 4.0);
  ASSERT_FLOAT_EQ(t3(2,1), 8.0);
  ASSERT_FLOAT_EQ(t3(2,2), 2.0);
  ASSERT_FLOAT_EQ(t3(2,3), 5.0);
  ASSERT_FLOAT_EQ(t3(3,1), 9.0);
  ASSERT_FLOAT_EQ(t3(3,2), 3.0);
  ASSERT_FLOAT_EQ(t3(3,3), 6.0);
  t3 = T3;
  t3.shift(-1);
  ASSERT_FLOAT_EQ(t3(1,1), 4.0);
  ASSERT_FLOAT_EQ(t3(1,2), 7.0);
  ASSERT_FLOAT_EQ(t3(1,3), 1.0);
  ASSERT_FLOAT_EQ(t3(2,1), 5.0);
  ASSERT_FLOAT_EQ(t3(2,2), 8.0);
  ASSERT_FLOAT_EQ(t3(2,3), 2.0);
  ASSERT_FLOAT_EQ(t3(3,1), 6.0);
  ASSERT_FLOAT_EQ(t3(3,2), 9.0);
  ASSERT_FLOAT_EQ(t3(3,3), 3.0);
  t3 = T3;
  t3.shift(2);
  ASSERT_FLOAT_EQ(t3(1,1), 4.0);
  ASSERT_FLOAT_EQ(t3(1,2), 7.0);
  ASSERT_FLOAT_EQ(t3(1,3), 1.0);
  ASSERT_FLOAT_EQ(t3(2,1), 5.0);
  ASSERT_FLOAT_EQ(t3(2,2), 8.0);
  ASSERT_FLOAT_EQ(t3(2,3), 2.0);
  ASSERT_FLOAT_EQ(t3(3,1), 6.0);
  ASSERT_FLOAT_EQ(t3(3,2), 9.0);
  ASSERT_FLOAT_EQ(t3(3,3), 3.0);
  t3 = T3;
  t3.shift(-2);
  ASSERT_FLOAT_EQ(t3(1,1), 7.0);
  ASSERT_FLOAT_EQ(t3(1,2), 1.0);
  ASSERT_FLOAT_EQ(t3(1,3), 4.0);
  ASSERT_FLOAT_EQ(t3(2,1), 8.0);
  ASSERT_FLOAT_EQ(t3(2,2), 2.0);
  ASSERT_FLOAT_EQ(t3(2,3), 5.0);
  ASSERT_FLOAT_EQ(t3(3,1), 9.0);
  ASSERT_FLOAT_EQ(t3(3,2), 3.0);
  ASSERT_FLOAT_EQ(t3(3,3), 6.0);
}


TEST(TestTensor, Principal)
{
  SymmTensor T1(2,true);
  SymmTensor T2(2);
  SymmTensor T3(3);

  T1(1,1) = 1.0; T1(2,2) = 2.0; T1(3,3) = 3.0;
  T2(1,1) = 1.0; T2(2,2) = 2.0;
  T3(1,1) = 1.0; T3(2,2) = 2.0; T3(3,3) = 3.0;

  double angle = 25.0*M_PI/180.0;
  Tensor Trans2(2);
  Trans2(1,1) = Trans2(2,2) = cos(angle);
  Trans2(1,2) = -sin(angle);
  Trans2(2,1) =  sin(angle);
  Tensor Trans3(Vec3(1.0,3.0,2.0));

  T1.transform(Trans2);
  T2.transform(Trans2);
  T3.transform(Trans3);

  Vec3 p;
  T1.principal(p);
  ASSERT_FLOAT_EQ(p.x, 3.0);
  ASSERT_FLOAT_EQ(p.y, 2.0);
  ASSERT_FLOAT_EQ(p.z, 1.0);
  T2.principal(p);
  ASSERT_FLOAT_EQ(p.x, 2.0);
  ASSERT_FLOAT_EQ(p.y, 1.0);
  T3.principal(p);
  ASSERT_FLOAT_EQ(p.x, 3.0);
  ASSERT_FLOAT_EQ(p.y, 2.0);
  ASSERT_FLOAT_EQ(p.z, 1.0);

  Vec3Vec dir(3);
  T1.principal(p,dir.data());
  ASSERT_FLOAT_EQ(p.x, 3.0);
  ASSERT_FLOAT_EQ(p.y, 2.0);
  ASSERT_FLOAT_EQ(p.z, 1.0);
  ASSERT_FLOAT_EQ(dir[0].x, 0.0);
  ASSERT_FLOAT_EQ(dir[0].y, 0.0);
  ASSERT_FLOAT_EQ(dir[0].z, 1.0);
  ASSERT_FLOAT_EQ(dir[1].x, Trans2(1,2));
  ASSERT_FLOAT_EQ(dir[1].y, Trans2(2,2));
  ASSERT_FLOAT_EQ(dir[1].z, 0.0);
  T2.principal(p,dir.data());
  ASSERT_FLOAT_EQ(p.x, 2.0);
  ASSERT_FLOAT_EQ(p.y, 1.0);
  ASSERT_FLOAT_EQ(dir[0].x, Trans2(1,2));
  ASSERT_FLOAT_EQ(dir[0].y, Trans2(2,2));
  ASSERT_FLOAT_EQ(dir[0].z, 0.0);
  ASSERT_FLOAT_EQ(dir[1].x,-Trans2(1,1));
  ASSERT_FLOAT_EQ(dir[1].y,-Trans2(2,1));
  ASSERT_FLOAT_EQ(dir[1].z, 0.0);
  T3.principal(p,dir.data());
  ASSERT_FLOAT_EQ(p.x, 3.0);
  ASSERT_FLOAT_EQ(p.y, 2.0);
  ASSERT_FLOAT_EQ(p.z, 1.0);
  ASSERT_NEAR(dir[0].x, Trans3(1,3), 1.0e-15);
  ASSERT_NEAR(dir[0].y, Trans3(2,3), 1.0e-15);
  ASSERT_NEAR(dir[0].z, Trans3(3,3), 1.0e-15);
  ASSERT_NEAR(dir[1].x, Trans3(1,2), 1.0e-15);
  ASSERT_NEAR(dir[1].y, Trans3(2,2), 1.0e-15);
  ASSERT_NEAR(dir[1].z, Trans3(3,2), 1.0e-15);
  ASSERT_NEAR(dir[2].x, Trans3(1,1), 1.0e-15);
  ASSERT_NEAR(dir[2].y, Trans3(2,1), 1.0e-15);
  ASSERT_NEAR(dir[2].z, Trans3(3,1), 1.0e-15);
  T3.principal(p,dir.data(),2);
  ASSERT_FLOAT_EQ(p.x, 3.0);
  ASSERT_FLOAT_EQ(p.y, 2.0);
  ASSERT_FLOAT_EQ(p.z, 1.0);
  ASSERT_NEAR(dir[0].x, Trans3(1,3), 1.0e-15);
  ASSERT_NEAR(dir[0].y, Trans3(2,3), 1.0e-15);
  ASSERT_NEAR(dir[0].z, Trans3(3,3), 1.0e-15);
  ASSERT_NEAR(dir[1].x, Trans3(1,1), 1.0e-15);
  ASSERT_NEAR(dir[1].y, Trans3(2,1), 1.0e-15);
  ASSERT_NEAR(dir[1].z, Trans3(3,1), 1.0e-15);
}
