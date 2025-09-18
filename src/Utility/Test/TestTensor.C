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

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;


TEST_CASE("TestTensor.Multiply")
{
  const double data[9] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };

  Tensor T2 = 2.0*Tensor(std::vector<double>(data,data+4));
  REQUIRE(T2.dim() == 2);
  REQUIRE_THAT(T2(1,1), WithinRel(2.0));
  REQUIRE_THAT(T2(2,1), WithinRel(4.0));
  REQUIRE_THAT(T2(1,2), WithinRel(6.0));
  REQUIRE_THAT(T2(2,2), WithinRel(8.0));

  Tensor T3 = 0.5*Tensor(std::vector<double>(data,data+9));
  REQUIRE(T3.dim() == 3);
  REQUIRE_THAT(T3(1,1), WithinRel(0.5));
  REQUIRE_THAT(T3(2,1), WithinRel(1.0));
  REQUIRE_THAT(T3(3,1), WithinRel(1.5));
  REQUIRE_THAT(T3(1,2), WithinRel(2.0));
  REQUIRE_THAT(T3(2,2), WithinRel(2.5));
  REQUIRE_THAT(T3(3,2), WithinRel(3.0));
  REQUIRE_THAT(T3(1,3), WithinRel(3.5));
  REQUIRE_THAT(T3(2,3), WithinRel(4.0));
  REQUIRE_THAT(T3(3,3), WithinRel(4.5));

  Tensor T = T2*T3;
  REQUIRE(T.dim() == 2);
  REQUIRE_THAT(T(1,1), WithinRel(7.0));
  REQUIRE_THAT(T(2,1), WithinRel(10.0));
  REQUIRE_THAT(T(1,2), WithinRel(19.0));
  REQUIRE_THAT(T(2,2), WithinRel(28.0));
}


TEST_CASE("TestTensor.Shift")
{
  const double data[9] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };

  Tensor T2(std::vector<double>(data,data+4));
  Tensor T3(std::vector<double>(data,data+9));

  Tensor t2(T2);
  t2.shift(1);
  REQUIRE_THAT(t2(1,1), WithinRel(3.0));
  REQUIRE_THAT(t2(1,2), WithinRel(1.0));
  REQUIRE_THAT(t2(2,1), WithinRel(4.0));
  REQUIRE_THAT(t2(2,2), WithinRel(2.0));
  t2 = T2;
  t2.shift(-1);
  REQUIRE_THAT(t2(1,1), WithinRel(3.0));
  REQUIRE_THAT(t2(1,2), WithinRel(1.0));
  REQUIRE_THAT(t2(2,1), WithinRel(4.0));
  REQUIRE_THAT(t2(2,2), WithinRel(2.0));
  t2 = T2;
  t2.shift(2);
  REQUIRE_THAT(t2(1,1), WithinRel(1.0));
  REQUIRE_THAT(t2(1,2), WithinRel(3.0));
  REQUIRE_THAT(t2(2,1), WithinRel(2.0));
  REQUIRE_THAT(t2(2,2), WithinRel(4.0));
  t2 = T2;
  t2.shift(-2);
  REQUIRE_THAT(t2(1,1), WithinRel(1.0));
  REQUIRE_THAT(t2(1,2), WithinRel(3.0));
  REQUIRE_THAT(t2(2,1), WithinRel(2.0));
  REQUIRE_THAT(t2(2,2), WithinRel(4.0));

  Tensor t3(T3);
  t3.shift(1);
  REQUIRE_THAT(t3(1,1), WithinRel(7.0));
  REQUIRE_THAT(t3(1,2), WithinRel(1.0));
  REQUIRE_THAT(t3(1,3), WithinRel(4.0));
  REQUIRE_THAT(t3(2,1), WithinRel(8.0));
  REQUIRE_THAT(t3(2,2), WithinRel(2.0));
  REQUIRE_THAT(t3(2,3), WithinRel(5.0));
  REQUIRE_THAT(t3(3,1), WithinRel(9.0));
  REQUIRE_THAT(t3(3,2), WithinRel(3.0));
  REQUIRE_THAT(t3(3,3), WithinRel(6.0));
  t3 = T3;
  t3.shift(-1);
  REQUIRE_THAT(t3(1,1), WithinRel(4.0));
  REQUIRE_THAT(t3(1,2), WithinRel(7.0));
  REQUIRE_THAT(t3(1,3), WithinRel(1.0));
  REQUIRE_THAT(t3(2,1), WithinRel(5.0));
  REQUIRE_THAT(t3(2,2), WithinRel(8.0));
  REQUIRE_THAT(t3(2,3), WithinRel(2.0));
  REQUIRE_THAT(t3(3,1), WithinRel(6.0));
  REQUIRE_THAT(t3(3,2), WithinRel(9.0));
  REQUIRE_THAT(t3(3,3), WithinRel(3.0));
  t3 = T3;
  t3.shift(2);
  REQUIRE_THAT(t3(1,1), WithinRel(4.0));
  REQUIRE_THAT(t3(1,2), WithinRel(7.0));
  REQUIRE_THAT(t3(1,3), WithinRel(1.0));
  REQUIRE_THAT(t3(2,1), WithinRel(5.0));
  REQUIRE_THAT(t3(2,2), WithinRel(8.0));
  REQUIRE_THAT(t3(2,3), WithinRel(2.0));
  REQUIRE_THAT(t3(3,1), WithinRel(6.0));
  REQUIRE_THAT(t3(3,2), WithinRel(9.0));
  REQUIRE_THAT(t3(3,3), WithinRel(3.0));
  t3 = T3;
  t3.shift(-2);
  REQUIRE_THAT(t3(1,1), WithinRel(7.0));
  REQUIRE_THAT(t3(1,2), WithinRel(1.0));
  REQUIRE_THAT(t3(1,3), WithinRel(4.0));
  REQUIRE_THAT(t3(2,1), WithinRel(8.0));
  REQUIRE_THAT(t3(2,2), WithinRel(2.0));
  REQUIRE_THAT(t3(2,3), WithinRel(5.0));
  REQUIRE_THAT(t3(3,1), WithinRel(9.0));
  REQUIRE_THAT(t3(3,2), WithinRel(3.0));
  REQUIRE_THAT(t3(3,3), WithinRel(6.0));
}


TEST_CASE("TestTensor.Principal")
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
  REQUIRE_THAT(p.x, WithinRel(3.0));
  REQUIRE_THAT(p.y, WithinRel(2.0));
  REQUIRE_THAT(p.z, WithinRel(1.0));
  T2.principal(p);
  REQUIRE_THAT(p.x, WithinRel(2.0));
  REQUIRE_THAT(p.y, WithinRel(1.0));
  T3.principal(p);
  REQUIRE_THAT(p.x, WithinRel(3.0));
  REQUIRE_THAT(p.y, WithinRel(2.0));
  REQUIRE_THAT(p.z, WithinRel(1.0));

  std::array<Vec3,3> dir;
  T1.principal(p,dir.data());
  REQUIRE_THAT(p.x, WithinRel(3.0));
  REQUIRE_THAT(p.y, WithinRel(2.0));
  REQUIRE_THAT(p.z, WithinRel(1.0));
  REQUIRE_THAT(dir[0].x, WithinAbs(0.0, 1e-14));
  REQUIRE_THAT(dir[0].y, WithinAbs(0.0, 1e-14));
  REQUIRE_THAT(dir[0].z, WithinRel(1.0));
  REQUIRE_THAT(dir[1].x, WithinRel(Trans2(1,2)));
  REQUIRE_THAT(dir[1].y, WithinRel(Trans2(2,2)));
  REQUIRE_THAT(dir[1].z, WithinAbs(0.0, 1e-14));
  T2.principal(p,dir.data());
  REQUIRE_THAT(p.x, WithinRel(2.0));
  REQUIRE_THAT(p.y, WithinRel(1.0));
  REQUIRE_THAT(dir[0].x, WithinRel(Trans2(1,2)));
  REQUIRE_THAT(dir[0].y, WithinRel(Trans2(2,2)));
  REQUIRE_THAT(dir[0].z, WithinAbs(0.0, 1e-14));
  REQUIRE_THAT(dir[1].x, WithinRel(-Trans2(1,1)));
  REQUIRE_THAT(dir[1].y, WithinRel(-Trans2(2,1)));
  REQUIRE_THAT(dir[1].z, WithinAbs(0.0, 1e-14));
  T3.principal(p,dir.data());
  REQUIRE_THAT(p.x, WithinRel(3.0));
  REQUIRE_THAT(p.y, WithinRel(2.0));
  REQUIRE_THAT(p.z, WithinRel(1.0));
  REQUIRE_THAT(dir[0].x, WithinRel(Trans3(1,3), 1.0e-15));
  REQUIRE_THAT(dir[0].y, WithinRel(Trans3(2,3), 1.0e-15));
  REQUIRE_THAT(dir[0].z, WithinRel(Trans3(3,3), 1.0e-15));
  REQUIRE_THAT(dir[1].x, WithinAbs(Trans3(1,2), 1.0e-15));
  REQUIRE_THAT(dir[1].y, WithinRel(Trans3(2,2), 1.0e-15));
  REQUIRE_THAT(dir[1].z, WithinRel(Trans3(3,2), 1.0e-15));
  REQUIRE_THAT(dir[2].x, WithinRel(Trans3(1,1), 1.0e-15));
  REQUIRE_THAT(dir[2].y, WithinRel(Trans3(2,1), 1.0e-15));
  REQUIRE_THAT(dir[2].z, WithinRel(Trans3(3,1), 1.0e-15));
  T3.principal(p,dir.data(),2);
  REQUIRE_THAT(p.x, WithinRel(3.0));
  REQUIRE_THAT(p.y, WithinRel(2.0));
  REQUIRE_THAT(p.z, WithinRel(1.0));
  REQUIRE_THAT(dir[0].x, WithinRel(Trans3(1,3), 1.0e-15));
  REQUIRE_THAT(dir[0].y, WithinRel(Trans3(2,3), 1.0e-15));
  REQUIRE_THAT(dir[0].z, WithinRel(Trans3(3,3), 1.0e-15));
  REQUIRE_THAT(dir[1].x, WithinRel(Trans3(1,1), 1.0e-15));
  REQUIRE_THAT(dir[1].y, WithinRel(Trans3(2,1), 1.0e-15));
  REQUIRE_THAT(dir[1].z, WithinRel(Trans3(3,1), 1.0e-15));

  REQUIRE(T1.vonMises() == sqrt(3.0));
  REQUIRE(T2.vonMises() == sqrt(3.0));
  REQUIRE(T3.vonMises() == sqrt(3.0));
}


TEST_CASE("TestTensor.Transform")
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
      REQUIRE_THAT(T1(i,j), WithinRel(S1(i,j)));
  for (i = 1; i <= 2; i++)
    for (j = i; j <= 2; j++)
      REQUIRE_THAT(T2(i,j), WithinRel(S2(i,j)));
  for (i = 1; i <= 3; i++)
    for (j = i; j <= 3; j++)
      REQUIRE_THAT(T3(i,j), WithinRel(S3(i,j)));
}


TEST_CASE("TestTensor.RotationAngles")
{
  Tensor Rot(0.4,1.6,1.1);
  Vec3 angles = Rot.rotVec();

  REQUIRE_THAT(angles.x, WithinRel(0.4));
  REQUIRE_THAT(angles.y, WithinRel(1.6));
  REQUIRE_THAT(angles.z, WithinRel(1.1));
}


TEST_CASE("TestTensor.Determinant")
{
  const double data[9] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };

  Tensor R1(std::vector<double>(data,data+1));
  Tensor R2(std::vector<double>(data,data+4));
  Tensor R3(std::vector<double>(data,data+9));

  REQUIRE(R1.dim() == 1);
  REQUIRE(R2.dim() == 2);
  REQUIRE(R3.dim() == 3);

  REQUIRE_THAT(R1.det(), WithinRel(1.0));
  REQUIRE_THAT(R2.det(), WithinRel(-2.0));
  REQUIRE_THAT(R3.det(), WithinAbs(0.0, 1e-14));

  SymmTensor S0(std::vector<double>(data,data+1));
  SymmTensor S1(std::vector<double>(data,data+3));
  SymmTensor S2(std::vector<double>(data,data+4));
  SymmTensor S3(std::vector<double>(data,data+6));

  REQUIRE(S0.dim() == 1);
  REQUIRE(S1.dim() == 2);
  REQUIRE(S2.dim() == 2);
  REQUIRE(S3.dim() == 3);

  REQUIRE_THAT(S0.det(), WithinRel(1.0));
  REQUIRE_THAT(S1.det(), WithinRel(-7.0));
  REQUIRE_THAT(S2.det(), WithinRel(-42.0));
  REQUIRE_THAT(S3.det(), WithinRel(101.0));

  Tensor T0(1), T1(2), T2(3), T3(3);
  T0 = S0; T1 = S1; T2 = S2; T3 = S3;

  REQUIRE_THAT(S0.det(), WithinRel(T0.det()));
  REQUIRE_THAT(S1.det(), WithinRel(T1.det()));
  REQUIRE_THAT(S2.det(), WithinRel(T2.det()));
  REQUIRE_THAT(S3.det(), WithinRel(T3.det()));

  T3 = R2;
  REQUIRE_THAT(T3(1,1), WithinRel(1.0));
  REQUIRE_THAT(T3(2,1), WithinRel(2.0));
  REQUIRE_THAT(T3(1,2), WithinRel(3.0));
  REQUIRE_THAT(T3(2,2), WithinRel(4.0));
  for (int i = 1; i <= 3; i++)
  {
    REQUIRE_THAT(T3(i,3), WithinAbs(0.0, 1e-14));
    REQUIRE_THAT(T3(3,i), WithinAbs(0.0, 1e-14));
  }
}


TEST_CASE("TestTensor.Rotate")
{
  const double alpha = 0.5; // approx. 29 degrees

  const double c = cos(alpha);
  const double s = sin(alpha);
  const double z = 0.0;

  Vec3 a2( z, c, s);
  Vec3 a3( z,-s, c);
  Vec3 b1( c, z, s);
  Vec3 b3(-s, z, c);
  Vec3 c1( c, s, z);
  Vec3 c2(-s, c, z);

  Vec3 i1(1.0, z, z);
  Vec3 i2(z, 1.0, z);
  Vec3 i3(z, z, 1.0);

  Tensor A(3,true), B(3,true);
  A.rotate(alpha,1); B.postMult(Tensor(i1,a2,a3));
  std::cout <<"A:\n"<< A <<"B:\n"<< B;
  REQUIRE(A.equal(B,0.0));
  A.rotate(alpha,2); B.postMult(Tensor(b1,i2,b3));
  std::cout <<"A:\n"<< A <<"B:\n"<< B;
  REQUIRE(A.equal(B,0.0));
  A.rotate(alpha,3); B.postMult(Tensor(c1,c2,i3));
  std::cout <<"A:\n"<< A <<"B:\n"<< B;
  REQUIRE(A.equal(B,0.0));

  Tensor C(2,true), D(2,true);
  C.rotate(alpha,3); D.postMult(Tensor(c1,c2,i3));
  std::cout <<"C:\n"<< C <<"D:\n"<< D;
  REQUIRE(C.equal(D,0.0));
}
