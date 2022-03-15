//==============================================================================
//!
//! \file TestFieldFunctionsLR.C
//!
//! \date Nov 7 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for parsing of LR field functions.
//!
//==============================================================================

#include "FieldFunctions.h"

#include "gtest/gtest.h"


TEST(TestFieldFunctions, 2D1PLR)
{
  // Fields contains: u_x = x^3*y^2 / 3
  //                  u_y = -x^2*y^3 / 3
  //                    p = x^2*y^2
  FieldFunction f2D_scalar("src/Utility/Test/refdata/Field2D-1PLR",
                           "Stokes-2", "p", 0);

  ScalarGradFieldFunction f2Dg_scalar("src/Utility/Test/refdata/Field2D-1PLR",
                                      "Stokes-2", "p", 0);

  VecFieldFunction f2D_vec("src/Utility/Test/refdata/Field2D-1PLR",
                           "Stokes-1", "u", 0);

  ScalarLaplacianFieldFunction f2Dl_scalar("src/Utility/Test/refdata/Field2D-1PLR",
                                           "Stokes-2", "p", 0);

  TensorFieldFunction f2D_ten("src/Utility/Test/refdata/Field2D-1PLR",
                              "Stokes-1", "u_x,x|u_x,y|u_y,x|u_y,y", 0);

  VecGradFieldFunction f2Dg_vec("src/Utility/Test/refdata/Field2D-1PLR",
                               "Stokes-1", "u", 0);

  VecLaplacianFieldFunction f2Dl_vec("src/Utility/Test/refdata/Field2D-1PLR",
                                     "Stokes-1", "u", 0);

  STensorFieldFunction f2D_sten("src/Utility/Test/refdata/Field2D-1PLR",
                                "Stokes-1", "u_x,x|u_y,y|u_x,y", 0);

  double param[3] = {0.5, 0.5, 0.0};
  Vec4 X(param);
  double scal = f2D_scalar(X);
  double x = 0.5, y = 0.5;
  EXPECT_NEAR(scal, x*x*y*y, 1e-14);

  Vec3 sg = f2Dg_scalar(X);
  EXPECT_NEAR(sg[0], 2*x*y*y, 1e-14);
  EXPECT_NEAR(sg[1], x*x*2*y, 1e-14);

  Vec3 vec = f2D_vec(X);
  EXPECT_NEAR(vec[0],  x*x*x*y*y/3.0, 1e-14);
  EXPECT_NEAR(vec[1], -x*x*y*y*y/3, 1e-14);

  Vec3 slap = f2Dl_scalar(X);
  EXPECT_NEAR(slap[0], 2*y*y, 1e-13);
  EXPECT_NEAR(slap[1], 2*x*x, 1e-13);

  Tensor vecg = f2Dg_vec(X);
  EXPECT_NEAR(vecg(1,1), x*x*y*y, 1e-14);
  EXPECT_NEAR(vecg(1,2), x*x*x*2.0*y/3.0, 1e-14);
  EXPECT_NEAR(vecg(2,1), -2.0*x*y*y*y/3.0, 1e-14);
  EXPECT_NEAR(vecg(2,2), -x*x*y*y, 1e-14);

  Tensor ten = f2D_ten(X);
  EXPECT_NEAR(vecg(1,1), x*x*y*y, 1e-14);
  EXPECT_NEAR(vecg(1,2), x*x*x*2.0*y/3.0, 1e-14);
  EXPECT_NEAR(vecg(2,1), -2.0*x*y*y*y/3.0, 1e-14);
  EXPECT_NEAR(vecg(2,2), -x*x*y*y, 1e-14);

  Tensor lap = f2Dl_vec(X);
  EXPECT_NEAR(lap(1,1), 2.0*x*y*y, 1e-13);
  EXPECT_NEAR(lap(1,2), x*x*x*2.0/3.0, 1e-13);
  EXPECT_NEAR(lap(2,1), -2.0*x*x*x/3.0, 1e-14);
  EXPECT_NEAR(lap(2,2), -x*x*2.0*y, 1e-14);

  SymmTensor sten = f2D_sten(X);
  EXPECT_NEAR(sten(1,1), x*x*y*y, 1e-14);
  EXPECT_NEAR(sten(1,2), x*x*x*2.0*y/3.0, 1e-14);
  EXPECT_NEAR(sten(2,1), x*x*x*2.0*y/3.0, 1e-14);
  EXPECT_NEAR(sten(2,2), -x*x*y*y, 1e-14);
}


TEST(TestFieldFunctions, 2D1PLRmx)
{
  // Fields contains: u_x = x^3*y^2 / 3
  //                  u_y = -x^2*y^3 / 3
  //                    p = x^2*y^2
  FieldFunction f2D_scalar("src/Utility/Test/refdata/Field2D-1PLRmx",
                           "Stokes-3", "p", 0);

  ScalarGradFieldFunction f2Dg_scalar("src/Utility/Test/refdata/Field2D-1PLRmx",
                                      "Stokes-3", "p", 0);

  VecFieldFunction f2D_vec("src/Utility/Test/refdata/Field2D-1PLRmx",
                           "Stokes", "u_x|u_y", 0);

  VecGradFieldFunction f2Dg_vec("src/Utility/Test/refdata/Field2D-1PLRmx",
                               "Stokes", "u_x|u_y", 0);

  VecLaplacianFieldFunction f2Dl_vec("src/Utility/Test/refdata/Field2D-1PLRmx",
                                     "Stokes", "u_x|u_y", 0);

  double param[3] = {0.5, 0.5, 0.0};
  Vec4 X(param);
  double scal = f2D_scalar(X);
  double x = 0.5, y = 0.5;
  EXPECT_NEAR(scal, x*x*y*y, 1e-12);

  Vec3 sg = f2Dg_scalar(X);
  EXPECT_NEAR(sg[0], 2*x*y*y, 1e-12);
  EXPECT_NEAR(sg[1], x*x*2*y, 1e-12);

  Vec3 vec = f2D_vec(X);
  EXPECT_NEAR(vec[0],  x*x*x*y*y/3.0, 1e-13);
  EXPECT_NEAR(vec[1], -x*x*y*y*y/3, 1e-13);

  Tensor vecg = f2Dg_vec(X);
  EXPECT_NEAR(vecg(1,1), x*x*y*y, 1e-13);
  EXPECT_NEAR(vecg(1,2), x*x*x*2.0*y/3.0, 1e-13);
  EXPECT_NEAR(vecg(2,1), -2.0*x*y*y*y/3.0, 1e-13);
  EXPECT_NEAR(vecg(2,2), -x*x*y*y, 1e-13);

  Tensor lap = f2Dl_vec(X);
  EXPECT_NEAR(lap(1,1), 2.0*x*y*y, 1e-12);
  EXPECT_NEAR(lap(1,2), x*x*x*2.0/3.0, 1e-12);
  EXPECT_NEAR(lap(2,1), -2.0*x*x*x/3.0, 1e-12);
  EXPECT_NEAR(lap(2,2), -x*x*2.0*y, 1e-12);
}


TEST(TestFieldFunctions, 2D2PLR)
{
  FieldFunction f2D_scalar("src/Utility/Test/refdata/Field2D-2PLR",
                           "Stokes-2", "v", 0);

  VecFieldFunction f2D_vec("src/Utility/Test/refdata/Field2D-2PLR",
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


TEST(TestFieldFunctions, 3D1PLR)
{
  // Fields contains: u_x = x^3*y^2*z^3 / 3
  //                  u_y = -x^2*y^3*z^3 / 3
  //                  u_z = -2*x^2*y^2*z^3 / 3
  //                    p = x^2*y^2*z^2
  FieldFunction f3D_scalar("src/Utility/Test/refdata/Field3D-1PLR",
                           "Stokes-2", "p", 0);

  ScalarGradFieldFunction f3Dg_scalar("src/Utility/Test/refdata/Field3D-1PLR",
                                      "Stokes-2", "p", 0);

  VecFieldFunction f3D_vec("src/Utility/Test/refdata/Field3D-1PLR",
                           "Stokes-1", "u", 0);

  ScalarLaplacianFieldFunction f3Dl_scalar("src/Utility/Test/refdata/Field3D-1PLR",
                                           "Stokes-2", "p", 0);

  TensorFieldFunction f3D_ten("src/Utility/Test/refdata/Field3D-1PLR",
                              "Stokes-1", "u_x,x|u_y,x|u_z,x|u_x,y|u_y,y|u_z,y|u_x,z|u_x,y|u_z,z", 0);

  VecGradFieldFunction f3Dg_vec("src/Utility/Test/refdata/Field3D-1PLR",
                                "Stokes-1", "u", 0);

  VecLaplacianFieldFunction f3Dl_vec("src/Utility/Test/refdata/Field3D-1PLR",
                                     "Stokes-1", "u", 0);

  STensorFieldFunction f3D_sten("src/Utility/Test/refdata/Field3D-1PLR",
                                "Stokes-1", "u_x,x|u_y,y|u_z,z|u_x,y|u_z,x|u_x,y", 0);

  double param[3] = {0.5, 0.5, 0.5};
  Vec4 X(param);
  double scal = f3D_scalar(X);
  double x = 0.5, y = 0.5, z = 0.5;
  EXPECT_NEAR(scal, x*x*y*y*z*z, 1e-13);

  Vec3 sg = f3Dg_scalar(X);
  EXPECT_NEAR(sg[0], 2.0*x*y*y*z*z, 1e-12);
  EXPECT_NEAR(sg[1], x*x*2.0*y*z*z, 1e-12);
  EXPECT_NEAR(sg[2], x*x*y*y*2.0*z, 1e-12);

  Vec3 vec = f3D_vec(X);
  EXPECT_NEAR(vec[0], x*x*x*y*y*z*z/3.0, 1e-14);
  EXPECT_NEAR(vec[1], x*x*y*y*y*z*z/3.0, 1e-14);
  EXPECT_NEAR(vec[2], -2.0*x*x*y*y*z*z*z/3.0, 1e-14);

  Vec3 slap = f3Dl_scalar(X);
  EXPECT_NEAR(slap[0], 2.0*y*y*z*z, 1e-11);
  EXPECT_NEAR(slap[1], x*x*2.0*z*z, 1e-11);
  EXPECT_NEAR(slap[2], x*x*y*y*2.0, 1e-11);

  Tensor vecg = f3Dg_vec(X);
  EXPECT_NEAR(vecg(1,1), x*x*y*y*z*z, 1e-13);
  EXPECT_NEAR(vecg(1,2), x*x*x*2.0*y*z*z / 3.0, 1e-13);
  EXPECT_NEAR(vecg(1,3), x*x*x*y*y*2.0*z / 3.0, 1e-13);
  EXPECT_NEAR(vecg(2,1), 2.0*x*y*y*y*z*z / 3.0, 1e-13);
  EXPECT_NEAR(vecg(2,2), x*x*y*y*z*z, 1e-13);
  EXPECT_NEAR(vecg(2,3), x*x*y*y*y*2.0*z / 3.0, 1e-13);
  EXPECT_NEAR(vecg(3,1), -2.0*2.0*x*y*y*z*z*z / 3.0, 1e-13);
  EXPECT_NEAR(vecg(3,2), -2.0*x*x*2.0*y*z*z*z / 3.0, 1e-13);
  EXPECT_NEAR(vecg(3,3), -2.0*x*x*y*y*z*z, 1e-13);

  Tensor ten = f3D_ten(X);
  EXPECT_NEAR(ten(1,1), x*x*y*y*z*z, 1e-13);
  EXPECT_NEAR(ten(1,2), x*x*x*2.0*y*z*z / 3.0, 1e-13);
  EXPECT_NEAR(ten(1,3), x*x*x*y*y*2.0*z / 3.0, 1e-13);
  EXPECT_NEAR(ten(2,1), 2*x*y*y*y*z*z / 3.0, 1e-13);
  EXPECT_NEAR(ten(2,2), x*x*y*y*z*z, 1e-13);
  EXPECT_NEAR(ten(2,3), x*x*y*y*y*2.0*z / 3.0, 1e-13);
  EXPECT_NEAR(ten(3,1), -2.0*2.0*x*y*y*z*z*z / 3.0, 1e-13);
  EXPECT_NEAR(ten(3,2), -2.0*x*x*2.0*y*z*z*z / 3.0, 1e-13);
  EXPECT_NEAR(ten(3,3), -2.0*x*x*y*y*z*z, 1e-13);

  Tensor lap = f3Dl_vec(X);
  EXPECT_NEAR(lap(1,1), 2.0*x*y*y*z*z, 1e-12);
  EXPECT_NEAR(lap(1,2), x*x*x*2.0*z*z / 3.0, 1e-12);
  EXPECT_NEAR(lap(1,3), x*x*x*y*y*2.0 / 3.0, 1e-12);
  EXPECT_NEAR(lap(2,1), 2.0*y*y*y*z*z / 3.0, 1e-12);
  EXPECT_NEAR(lap(2,2), x*x*2.0*y*z*z, 1e-12);
  EXPECT_NEAR(lap(2,3), x*x*y*y*y*2.0 / 3.0, 1e-12);
  EXPECT_NEAR(lap(3,1), -2.0*2.0*y*y*z*z*z / 3.0, 1e-11);
  EXPECT_NEAR(lap(3,2), -2.0*x*x*2.0*z*z*z / 3.0, 1e-11);
  EXPECT_NEAR(lap(3,3), -2.0*x*x*y*y*2.0*z, 1e-11);

  SymmTensor sten = f3D_sten(X);
  EXPECT_NEAR(sten(1,1), x*x*y*y*z*z, 1e-13);
  EXPECT_NEAR(sten(1,2), x*x*x*2.0*y*z*z / 3.0, 1e-13);
  EXPECT_NEAR(sten(1,3), x*x*x*y*y*2.0*z / 3.0, 1e-13);
  EXPECT_NEAR(sten(2,1), 2.0*x*y*y*y*z*z / 3.0, 1e-13);
  EXPECT_NEAR(sten(2,2), x*x*y*y*z*z, 1e-13);
  EXPECT_NEAR(sten(2,3), -2.0*x*x*2.0*y*z*z*z / 3.0, 1e-13);
  EXPECT_NEAR(sten(3,1), x*x*x*y*y*2.0*z / 3.0, 1e-13);
  EXPECT_NEAR(sten(3,2), -2.0*x*x*2.0*y*z*z*z / 3.0, 1e-13);
  EXPECT_NEAR(sten(3,3), -2.0*x*x*y*y*z*z, 1e-13);
}


TEST(TestFieldFunctions, 3D1PLRmx)
{
  // Fields contains: u_x = x^3*y^2*z^3 / 3
  //                  u_y = -x^2*y^3*z^3 / 3
  //                  u_z = -2*x^2*y^2*z^3 / 3
  //                    p = x^2*y^2*z^2
  FieldFunction f3D_scalar("src/Utility/Test/refdata/Field3D-1PLRmx",
                           "Stokes-4", "p", 0);

  ScalarGradFieldFunction f3Dg_scalar("src/Utility/Test/refdata/Field3D-1PLRmx",
                                      "Stokes-4", "p", 0);

  VecFieldFunction f3D_vec("src/Utility/Test/refdata/Field3D-1PLRmx",
                           "Stokes", "u_x|u_y|u_z", 0);

  ScalarLaplacianFieldFunction f3Dl_scalar("src/Utility/Test/refdata/Field3D-1PLRmx",
                                           "Stokes-4", "p", 0);

  VecGradFieldFunction f3Dg_vec("src/Utility/Test/refdata/Field3D-1PLRmx",
                                "Stokes", "u_x|u_y|u_z", 0);

  VecLaplacianFieldFunction f3Dl_vec("src/Utility/Test/refdata/Field3D-1PLRmx",
                                     "Stokes", "u_x|u_y|u_z", 0);

  double param[3] = {0.5, 0.5, 0.5};
  Vec4 X(param);
  double scal = f3D_scalar(X);
  double x = 0.5, y = 0.5, z = 0.5;
  EXPECT_NEAR(scal, x*x*y*y*z*z, 1e-14);

  Vec3 sg = f3Dg_scalar(X);
  EXPECT_NEAR(sg[0], 2.0*x*y*y*z*z, 1e-13);
  EXPECT_NEAR(sg[1], x*x*2.0*y*z*z, 1e-13);
  EXPECT_NEAR(sg[2], x*x*y*y*2.0*z, 1e-13);

  Vec3 vec = f3D_vec(X);
  EXPECT_NEAR(vec[0], x*x*x*y*y*z*z/3.0, 1e-14);
  EXPECT_NEAR(vec[1], x*x*y*y*y*z*z/3.0, 1e-14);
  EXPECT_NEAR(vec[2], -2.0*x*x*y*y*z*z*z/3.0, 1e-14);

  Vec3 slap = f3Dl_scalar(X);
  EXPECT_NEAR(slap[0], 2.0*y*y*z*z, 1e-13);
  EXPECT_NEAR(slap[1], x*x*2.0*z*z, 1e-13);
  EXPECT_NEAR(slap[2], x*x*y*y*2.0, 1e-13);

  Tensor vecg = f3Dg_vec(X);
  EXPECT_NEAR(vecg(1,1), x*x*y*y*z*z, 1e-14);
  EXPECT_NEAR(vecg(1,2), x*x*x*2.0*y*z*z / 3.0, 1e-14);
  EXPECT_NEAR(vecg(1,3), x*x*x*y*y*2.0*z / 3.0, 1e-14);
  EXPECT_NEAR(vecg(2,1), 2.0*x*y*y*y*z*z / 3.0, 1e-14);
  EXPECT_NEAR(vecg(2,2), x*x*y*y*z*z, 1e-14);
  EXPECT_NEAR(vecg(2,3), x*x*y*y*y*2.0*z / 3.0, 1e-14);
  EXPECT_NEAR(vecg(3,1), -2.0*2.0*x*y*y*z*z*z / 3.0, 1e-14);
  EXPECT_NEAR(vecg(3,2), -2.0*x*x*2.0*y*z*z*z / 3.0, 1e-14);
  EXPECT_NEAR(vecg(3,3), -2.0*x*x*y*y*z*z, 1e-14);

  Tensor lap = f3Dl_vec(X);
  EXPECT_NEAR(lap(1,1), 2.0*x*y*y*z*z, 1e-13);
  EXPECT_NEAR(lap(1,2), x*x*x*2.0*z*z / 3.0, 1e-13);
  EXPECT_NEAR(lap(1,3), x*x*x*y*y*2.0 / 3.0, 1e-13);
  EXPECT_NEAR(lap(2,1), 2.0*y*y*y*z*z / 3.0, 1e-13);
  EXPECT_NEAR(lap(2,2), x*x*2.0*y*z*z, 1e-13);
  EXPECT_NEAR(lap(2,3), x*x*y*y*y*2.0 / 3.0, 1e-13);
  EXPECT_NEAR(lap(3,1), -2.0*2.0*y*y*z*z*z / 3.0, 1e-13);
  EXPECT_NEAR(lap(3,2), -2.0*x*x*2.0*z*z*z / 3.0, 1e-13);
  EXPECT_NEAR(lap(3,3), -2.0*x*x*y*y*2.0*z, 1e-13);
}


TEST(TestFieldFunctions, 3D2PLR)
{
  FieldFunction f3D_scalar("src/Utility/Test/refdata/Field3D-2PLR",
                           "Stokes-2", "v", 0);

  VecFieldFunction f3D_vec("src/Utility/Test/refdata/Field3D-2PLR",
                           "Stokes-1", "u", 0);

  TensorFieldFunction f3D_ten("src/Utility/Test/refdata/Field3D-2PLR",
                              "Stokes-1", "u_x,x|u_y,x|u_z,x|u_x,y|u_y,y|u_z,y|u_x,z|u_x,y|u_z,z", 0);

  STensorFieldFunction f3D_sten("src/Utility/Test/refdata/Field3D-2PLR",
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
