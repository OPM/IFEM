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
#include <array>

#include "gtest/gtest.h"


TEST(TestTensor, Multiply)
{
  const double data[9] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };

  Tensor T2 = 2.0*Tensor(std::vector<double>(data,data+4));
  ASSERT_EQ(T2.dim(), 2U);
  EXPECT_FLOAT_EQ(T2(1,1), 2.0);
  EXPECT_FLOAT_EQ(T2(2,1), 4.0);
  EXPECT_FLOAT_EQ(T2(1,2), 6.0);
  EXPECT_FLOAT_EQ(T2(2,2), 8.0);

  Tensor T3 = 0.5*Tensor(std::vector<double>(data,data+9));
  ASSERT_EQ(T3.dim(), 3U);
  EXPECT_FLOAT_EQ(T3(1,1), 0.5);
  EXPECT_FLOAT_EQ(T3(2,1), 1.0);
  EXPECT_FLOAT_EQ(T3(3,1), 1.5);
  EXPECT_FLOAT_EQ(T3(1,2), 2.0);
  EXPECT_FLOAT_EQ(T3(2,2), 2.5);
  EXPECT_FLOAT_EQ(T3(3,2), 3.0);
  EXPECT_FLOAT_EQ(T3(1,3), 3.5);
  EXPECT_FLOAT_EQ(T3(2,3), 4.0);
  EXPECT_FLOAT_EQ(T3(3,3), 4.5);

  Tensor T = T2*T3;
  ASSERT_EQ(T.dim(), 2U);
  EXPECT_FLOAT_EQ(T(1,1),  7.0);
  EXPECT_FLOAT_EQ(T(2,1), 10.0);
  EXPECT_FLOAT_EQ(T(1,2), 19.0);
  EXPECT_FLOAT_EQ(T(2,2), 28.0);
}


TEST(TestTensor, Shift)
{
  const double data[9] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };

  Tensor T2(std::vector<double>(data,data+4));
  Tensor T3(std::vector<double>(data,data+9));

  Tensor t2(T2);
  t2.shift(1);
  EXPECT_FLOAT_EQ(t2(1,1), 3.0);
  EXPECT_FLOAT_EQ(t2(1,2), 1.0);
  EXPECT_FLOAT_EQ(t2(2,1), 4.0);
  EXPECT_FLOAT_EQ(t2(2,2), 2.0);
  t2 = T2;
  t2.shift(-1);
  EXPECT_FLOAT_EQ(t2(1,1), 3.0);
  EXPECT_FLOAT_EQ(t2(1,2), 1.0);
  EXPECT_FLOAT_EQ(t2(2,1), 4.0);
  EXPECT_FLOAT_EQ(t2(2,2), 2.0);
  t2 = T2;
  t2.shift(2);
  EXPECT_FLOAT_EQ(t2(1,1), 1.0);
  EXPECT_FLOAT_EQ(t2(1,2), 3.0);
  EXPECT_FLOAT_EQ(t2(2,1), 2.0);
  EXPECT_FLOAT_EQ(t2(2,2), 4.0);
  t2 = T2;
  t2.shift(-2);
  EXPECT_FLOAT_EQ(t2(1,1), 1.0);
  EXPECT_FLOAT_EQ(t2(1,2), 3.0);
  EXPECT_FLOAT_EQ(t2(2,1), 2.0);
  EXPECT_FLOAT_EQ(t2(2,2), 4.0);

  Tensor t3(T3);
  t3.shift(1);
  EXPECT_FLOAT_EQ(t3(1,1), 7.0);
  EXPECT_FLOAT_EQ(t3(1,2), 1.0);
  EXPECT_FLOAT_EQ(t3(1,3), 4.0);
  EXPECT_FLOAT_EQ(t3(2,1), 8.0);
  EXPECT_FLOAT_EQ(t3(2,2), 2.0);
  EXPECT_FLOAT_EQ(t3(2,3), 5.0);
  EXPECT_FLOAT_EQ(t3(3,1), 9.0);
  EXPECT_FLOAT_EQ(t3(3,2), 3.0);
  EXPECT_FLOAT_EQ(t3(3,3), 6.0);
  t3 = T3;
  t3.shift(-1);
  EXPECT_FLOAT_EQ(t3(1,1), 4.0);
  EXPECT_FLOAT_EQ(t3(1,2), 7.0);
  EXPECT_FLOAT_EQ(t3(1,3), 1.0);
  EXPECT_FLOAT_EQ(t3(2,1), 5.0);
  EXPECT_FLOAT_EQ(t3(2,2), 8.0);
  EXPECT_FLOAT_EQ(t3(2,3), 2.0);
  EXPECT_FLOAT_EQ(t3(3,1), 6.0);
  EXPECT_FLOAT_EQ(t3(3,2), 9.0);
  EXPECT_FLOAT_EQ(t3(3,3), 3.0);
  t3 = T3;
  t3.shift(2);
  EXPECT_FLOAT_EQ(t3(1,1), 4.0);
  EXPECT_FLOAT_EQ(t3(1,2), 7.0);
  EXPECT_FLOAT_EQ(t3(1,3), 1.0);
  EXPECT_FLOAT_EQ(t3(2,1), 5.0);
  EXPECT_FLOAT_EQ(t3(2,2), 8.0);
  EXPECT_FLOAT_EQ(t3(2,3), 2.0);
  EXPECT_FLOAT_EQ(t3(3,1), 6.0);
  EXPECT_FLOAT_EQ(t3(3,2), 9.0);
  EXPECT_FLOAT_EQ(t3(3,3), 3.0);
  t3 = T3;
  t3.shift(-2);
  EXPECT_FLOAT_EQ(t3(1,1), 7.0);
  EXPECT_FLOAT_EQ(t3(1,2), 1.0);
  EXPECT_FLOAT_EQ(t3(1,3), 4.0);
  EXPECT_FLOAT_EQ(t3(2,1), 8.0);
  EXPECT_FLOAT_EQ(t3(2,2), 2.0);
  EXPECT_FLOAT_EQ(t3(2,3), 5.0);
  EXPECT_FLOAT_EQ(t3(3,1), 9.0);
  EXPECT_FLOAT_EQ(t3(3,2), 3.0);
  EXPECT_FLOAT_EQ(t3(3,3), 6.0);
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
  EXPECT_FLOAT_EQ(p.x, 3.0);
  EXPECT_FLOAT_EQ(p.y, 2.0);
  EXPECT_FLOAT_EQ(p.z, 1.0);
  T2.principal(p);
  EXPECT_FLOAT_EQ(p.x, 2.0);
  EXPECT_FLOAT_EQ(p.y, 1.0);
  T3.principal(p);
  EXPECT_FLOAT_EQ(p.x, 3.0);
  EXPECT_FLOAT_EQ(p.y, 2.0);
  EXPECT_FLOAT_EQ(p.z, 1.0);

  std::array<Vec3,3> dir;
  T1.principal(p,dir.data());
  EXPECT_FLOAT_EQ(p.x, 3.0);
  EXPECT_FLOAT_EQ(p.y, 2.0);
  EXPECT_FLOAT_EQ(p.z, 1.0);
  EXPECT_FLOAT_EQ(dir[0].x, 0.0);
  EXPECT_FLOAT_EQ(dir[0].y, 0.0);
  EXPECT_FLOAT_EQ(dir[0].z, 1.0);
  EXPECT_FLOAT_EQ(dir[1].x, Trans2(1,2));
  EXPECT_FLOAT_EQ(dir[1].y, Trans2(2,2));
  EXPECT_FLOAT_EQ(dir[1].z, 0.0);
  T2.principal(p,dir.data());
  EXPECT_FLOAT_EQ(p.x, 2.0);
  EXPECT_FLOAT_EQ(p.y, 1.0);
  EXPECT_FLOAT_EQ(dir[0].x, Trans2(1,2));
  EXPECT_FLOAT_EQ(dir[0].y, Trans2(2,2));
  EXPECT_FLOAT_EQ(dir[0].z, 0.0);
  EXPECT_FLOAT_EQ(dir[1].x,-Trans2(1,1));
  EXPECT_FLOAT_EQ(dir[1].y,-Trans2(2,1));
  EXPECT_FLOAT_EQ(dir[1].z, 0.0);
  T3.principal(p,dir.data());
  EXPECT_FLOAT_EQ(p.x, 3.0);
  EXPECT_FLOAT_EQ(p.y, 2.0);
  EXPECT_FLOAT_EQ(p.z, 1.0);
  EXPECT_NEAR(dir[0].x, Trans3(1,3), 1.0e-15);
  EXPECT_NEAR(dir[0].y, Trans3(2,3), 1.0e-15);
  EXPECT_NEAR(dir[0].z, Trans3(3,3), 1.0e-15);
  EXPECT_NEAR(dir[1].x, Trans3(1,2), 1.0e-15);
  EXPECT_NEAR(dir[1].y, Trans3(2,2), 1.0e-15);
  EXPECT_NEAR(dir[1].z, Trans3(3,2), 1.0e-15);
  EXPECT_NEAR(dir[2].x, Trans3(1,1), 1.0e-15);
  EXPECT_NEAR(dir[2].y, Trans3(2,1), 1.0e-15);
  EXPECT_NEAR(dir[2].z, Trans3(3,1), 1.0e-15);
  T3.principal(p,dir.data(),2);
  EXPECT_FLOAT_EQ(p.x, 3.0);
  EXPECT_FLOAT_EQ(p.y, 2.0);
  EXPECT_FLOAT_EQ(p.z, 1.0);
  EXPECT_NEAR(dir[0].x, Trans3(1,3), 1.0e-15);
  EXPECT_NEAR(dir[0].y, Trans3(2,3), 1.0e-15);
  EXPECT_NEAR(dir[0].z, Trans3(3,3), 1.0e-15);
  EXPECT_NEAR(dir[1].x, Trans3(1,1), 1.0e-15);
  EXPECT_NEAR(dir[1].y, Trans3(2,1), 1.0e-15);
  EXPECT_NEAR(dir[1].z, Trans3(3,1), 1.0e-15);

  EXPECT_FLOAT_EQ(T1.vonMises(), sqrt(3.0));
  EXPECT_FLOAT_EQ(T2.vonMises(), sqrt(3.0));
  EXPECT_FLOAT_EQ(T3.vonMises(), sqrt(3.0));
}


TEST(TestTensor, Transform)
{
  const double data[9] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };

  SymmTensor T1(std::vector<double>(data,data+4));
  SymmTensor T2(std::vector<double>(data,data+3));
  SymmTensor T3(std::vector<double>(data,data+6));
  SymmTensor S1(T1), S2(T2), S3(T3);

  double angle = 25.0*M_PI/180.0;
  Tensor Trans2(2), Trans3(3);
  Trans2(1,1) = Trans2(2,2) = cos(angle);
  Trans2(1,2) = -sin(angle);
  Trans2(2,1) =  sin(angle);
  Trans3 = Trans2; Trans3(3,3) = 1.0;

  T1.transform(Trans2);
  T2.transform(Trans2);
  T3.transform(Trans2);
  S1.transform(Trans3);
  S2.transform(Trans3);
  S3.transform(Trans3);

  size_t i, j;
  for (i = 1; i <= 3; i++)
    for (j = i; j <= 3; j++)
      EXPECT_FLOAT_EQ(T1(i,j),S1(i,j));
  for (i = 1; i <= 2; i++)
    for (j = i; j <= 2; j++)
      EXPECT_FLOAT_EQ(T2(i,j),S2(i,j));
  for (i = 1; i <= 3; i++)
    for (j = i; j <= 3; j++)
      EXPECT_FLOAT_EQ(T3(i,j),S3(i,j));
}


TEST(TestTensor, RotationAngles)
{
  Tensor Rot(0.4,1.6,1.1);
  Vec3 angles = Rot.rotVec();

  EXPECT_FLOAT_EQ(angles.x, 0.4);
  EXPECT_FLOAT_EQ(angles.y, 1.6);
  EXPECT_FLOAT_EQ(angles.z, 1.1);
}


TEST(TestTensor, Determinant)
{
  const double data[9] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };

  Tensor R1(std::vector<double>(data,data+1));
  Tensor R2(std::vector<double>(data,data+4));
  Tensor R3(std::vector<double>(data,data+9));

  EXPECT_EQ(R1.dim(),1);
  EXPECT_EQ(R2.dim(),2);
  EXPECT_EQ(R3.dim(),3);

  EXPECT_FLOAT_EQ(R1.det(), 1.0);
  EXPECT_FLOAT_EQ(R2.det(),-2.0);
  EXPECT_FLOAT_EQ(R3.det(), 0.0);

  SymmTensor S0(std::vector<double>(data,data+1));
  SymmTensor S1(std::vector<double>(data,data+3));
  SymmTensor S2(std::vector<double>(data,data+4));
  SymmTensor S3(std::vector<double>(data,data+6));

  EXPECT_EQ(S0.dim(),1);
  EXPECT_EQ(S1.dim(),2);
  EXPECT_EQ(S2.dim(),2);
  EXPECT_EQ(S3.dim(),3);

  EXPECT_FLOAT_EQ(S0.det(),  1.0);
  EXPECT_FLOAT_EQ(S1.det(), -7.0);
  EXPECT_FLOAT_EQ(S2.det(),-42.0);
  EXPECT_FLOAT_EQ(S3.det(),101.0);

  Tensor T0(1), T1(2), T2(3), T3(3);
  T0 = S0; T1 = S1; T2 = S2; T3 = S3;

  EXPECT_FLOAT_EQ(S0.det(),T0.det());
  EXPECT_FLOAT_EQ(S1.det(),T1.det());
  EXPECT_FLOAT_EQ(S2.det(),T2.det());
  EXPECT_FLOAT_EQ(S3.det(),T3.det());
}
