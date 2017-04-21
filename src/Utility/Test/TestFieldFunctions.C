//==============================================================================
//!
//! \file TestFieldFunctions.C
//!
//! \date Nov 7 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for parsing of field functions.
//!
//==============================================================================

#include "FieldFunctions.h"

#include "gtest/gtest.h"


TEST(TestFieldFunctions, 2D1P)
{
  FieldFunction f2D_scalar("src/Utility/Test/refdata/Field2D-1P",
                           "Stokes-2", "v", 0);

  VecFieldFunction f2D_vec("src/Utility/Test/refdata/Field2D-1P",
                           "Stokes-1", "u", 0);

  TensorFieldFunction f2D_ten("src/Utility/Test/refdata/Field2D-1P",
                              "Stokes-1", "u_x,x|u_x,y|u_y,x|u_y,y", 0);

  STensorFieldFunction f2D_sten("src/Utility/Test/refdata/Field2D-1P",
                                "Stokes-1", "u_x,x|u_y,y|u_x,y", 0);

  double param[3] = {0.5, 0.5, 0.0};
  Vec4 X(param);
  double scal = f2D_scalar(X);
  double x = 1.0, y = 1.0;
  EXPECT_NEAR(scal, x, 1e-14);
  Vec3 vec = f2D_vec(X);
  EXPECT_NEAR(vec[0],  x, 1e-14);
  EXPECT_NEAR(vec[1], -y, 1e-14);
  Tensor ten = f2D_ten(X); // see above
  EXPECT_NEAR(ten(1,1), 1.0, 1e-14);
  EXPECT_NEAR(ten(1,2), 0.0, 1e-14);
  EXPECT_NEAR(ten(2,1), 0.0, 1e-14);
  EXPECT_NEAR(ten(2,2), -1.0, 1e-14);
  SymmTensor sten = f2D_sten(X);
  EXPECT_NEAR(sten(1,1), 1.0, 1e-14);
  EXPECT_NEAR(sten(1,2), 0.0, 1e-14);
  EXPECT_NEAR(sten(2,1), 0.0, 1e-14);
  EXPECT_NEAR(sten(2,2), -1.0, 1e-14);

  param[0] = 1.0;
  param[1] = 0.25;
  scal = f2D_scalar(X);
  x = 2.0, y = 0.5;
  vec = f2D_vec(X);
  EXPECT_NEAR(scal, x, 1e-14);
  EXPECT_NEAR(vec[0], x, 1e-14);
  EXPECT_NEAR(vec[1], -y, 1e-14);
}


TEST(TestFieldFunctions, 2D2P)
{
  FieldFunction f2D_scalar("src/Utility/Test/refdata/Field2D-2P",
                           "Stokes-2", "v", 0);

  VecFieldFunction f2D_vec("src/Utility/Test/refdata/Field2D-2P",
                           "Stokes-1", "u", 0);

  double param[3] = {0.5, 0.5, 0.0};
  Vec4 X(param);
  for (size_t pIdx = 0; pIdx < 2; ++pIdx) {
    f2D_scalar.initPatch(pIdx);
    double scal = f2D_scalar(X);
    EXPECT_NEAR(scal, 0.5 + pIdx, 1e-14);
    f2D_vec.initPatch(pIdx);
    Vec3 vec = f2D_vec(X);
    EXPECT_NEAR(vec[0], 0.5 + pIdx, 1e-14);
    EXPECT_NEAR(vec[1], -1.0, 1e-14);
  }
}


TEST(TestFieldFunctions, 3D1P)
{
  FieldFunction f3D_scalar("src/Utility/Test/refdata/Field3D-1P",
                           "Stokes-2", "v", 0);

  VecFieldFunction f3D_vec("src/Utility/Test/refdata/Field3D-1P",
                           "Stokes-1", "u", 0);

  TensorFieldFunction f3D_ten("src/Utility/Test/refdata/Field3D-1P",
                              "Stokes-1", "u_x,x|u_y,x|u_z,x|u_x,y|u_y,y|u_z,y|u_x,z|u_x,y|u_z,z", 0);

  STensorFieldFunction f3D_sten("src/Utility/Test/refdata/Field3D-1P",
                                "Stokes-1", "u_x,x|u_y,y|u_z,z|u_x,y|u_z,x|u_x,y", 0);

  double param[3] = {0.25, 0.25, 0.25};
  Vec4 X(param);
  double scal = f3D_scalar(X);
  double x = 0.5;
  double y = 0.5;
  double z = 0.5;
  EXPECT_NEAR(scal, 0.0, 1e-14);
  Vec3 vec = f3D_vec(X);
  EXPECT_NEAR(vec[0], y*(1.0-y), 1e-14); // y*(1-y)
  EXPECT_NEAR(vec[1], z*(1.0-z), 1e-14); // z*(1-z)
  EXPECT_NEAR(vec[2], x*(1.0-x), 1e-14); // x*(1-x)
  param[0] = param[1] = param[2] = 0.5;
  Tensor ten = f3D_ten(X);
  x = y = z = 1.0;
  EXPECT_NEAR(ten(1,1),  0.0, 1e-14);
  EXPECT_NEAR(ten(1,2), 1.0-2*y, 1e-14);
  EXPECT_NEAR(ten(1,3),  0.0, 1e-14);
  EXPECT_NEAR(ten(2,1),  0.0, 1e-14);
  EXPECT_NEAR(ten(2,2),  0.0, 1e-14);
  EXPECT_NEAR(ten(2,3), 1.0-2*z, 1e-14);
  EXPECT_NEAR(ten(3,1), 1.0-2*x, 1e-14);
  EXPECT_NEAR(ten(3,2),  0.0, 1e-14);
  EXPECT_NEAR(ten(3,3),  0.0, 1e-14);

  SymmTensor sten = f3D_sten(X);
  EXPECT_NEAR(sten(1,1),  0.0, 1e-14);
  EXPECT_NEAR(sten(1,2), 1.0-2*x, 1e-14);
  EXPECT_NEAR(sten(1,3), 1.0-2*y, 1e-14);
  EXPECT_NEAR(sten(2,1), 1.0-2*x, 1e-14);
  EXPECT_NEAR(sten(2,2),  0.0, 1e-14);
  EXPECT_NEAR(sten(2,3), 1.0-2*x, 1e-14);
  EXPECT_NEAR(sten(3,1), 1.0-2*y, 1e-14);
  EXPECT_NEAR(sten(3,2), 1.0-2*x, 1e-14);
  EXPECT_NEAR(sten(3,3),  0.0, 1e-14);
}


TEST(TestFieldFunctions, 3D2P)
{
  FieldFunction f3D_scalar("src/Utility/Test/refdata/Field3D-2P",
                           "Stokes-2", "v", 0);

  VecFieldFunction f3D_vec("src/Utility/Test/refdata/Field3D-2P",
                           "Stokes-1", "u", 0);

  TensorFieldFunction f3D_ten("src/Utility/Test/refdata/Field3D-2P",
                              "Stokes-1", "u_x,x|u_y,x|u_z,x|u_x,y|u_y,y|u_z,y|u_x,z|u_x,y|u_z,z", 0);

  STensorFieldFunction f3D_sten("src/Utility/Test/refdata/Field3D-2P",
                                "Stokes-1", "u_x,x|u_y,y|u_z,z|u_x,y|u_y,z|u_z,x", 0);

  double param[3] = {0.25, 0.25, 0.25};
  Vec4 X(param);
  for (size_t pIdx = 0; pIdx < 2; ++pIdx) {
    f3D_scalar.initPatch(pIdx);
    double scal = f3D_scalar(X);
    double x = 0.25 + pIdx;
    double y = 0.5;
    double z = 0.5;
    EXPECT_NEAR(scal, 1.0-2*x, 1e-14);
    f3D_vec.initPatch(pIdx);
    Vec3 vec = f3D_vec(X);
    EXPECT_NEAR(vec[0], y*(1.0-y), 1e-14); // y*(1-y)
    EXPECT_NEAR(vec[1], z*(1.0-z), 1e-14); // z*(1-z)
    EXPECT_NEAR(vec[2], x*(1.0-x), 1e-14); // x*(1-x)
    f3D_ten.initPatch(pIdx);
    Tensor ten = f3D_ten(X);
    EXPECT_NEAR(ten(1,1),  0.0, 1e-14);
    EXPECT_NEAR(ten(1,2), 1.0-2*y, 1e-14);
    EXPECT_NEAR(ten(1,3),  0.0, 1e-14);
    EXPECT_NEAR(ten(2,1),  0.0, 1e-14);
    EXPECT_NEAR(ten(2,2),  0.0, 1e-14);
    EXPECT_NEAR(ten(2,3), 1.0-2*z, 1e-14);
    EXPECT_NEAR(ten(3,1), 1.0-2*x, 1e-14);
    EXPECT_NEAR(ten(3,2),  0.0, 1e-14);
    EXPECT_NEAR(ten(3,3),  0.0, 1e-14);

    f3D_sten.initPatch(pIdx);
    SymmTensor sten = f3D_sten(X);
    EXPECT_NEAR(sten(1,1),  0.0, 1e-14);
    EXPECT_NEAR(sten(1,2), 1.0-2*y, 1e-14);
    EXPECT_NEAR(sten(1,3), 1.0-2*x, 1e-14);
    EXPECT_NEAR(sten(2,1), 1.0-2*y, 1e-14);
    EXPECT_NEAR(sten(2,2),  0.0, 1e-14);
    EXPECT_NEAR(sten(2,3), 1.0-2*z, 1e-14);
    EXPECT_NEAR(sten(3,1), 1.0-2*x, 1e-14);
    EXPECT_NEAR(sten(3,2), 1.0-2*z, 1e-14);
    EXPECT_NEAR(sten(3,3),  0.0, 1e-14);
  }
}
