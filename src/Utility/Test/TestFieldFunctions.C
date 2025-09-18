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

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;


TEST_CASE("TestFieldFunctions.2D1P")
{
  // Fields contains: u_x = x^3*y^2 / 3
  //                  u_y = -x^2*y^3 / 3
  //                    p = x^2*y^2
  FieldFunction f2D_scalar("src/Utility/Test/refdata/Field2D-1P",
                           "Stokes-2", "p", 0);

  ScalarGradFieldFunction f2Dg_scalar("src/Utility/Test/refdata/Field2D-1P",
                                      "Stokes-2", "p", 0);

  VecFieldFunction f2D_vec("src/Utility/Test/refdata/Field2D-1P",
                           "Stokes-1", "u", 0);

  ScalarLaplacianFieldFunction f2Dl_scalar("src/Utility/Test/refdata/Field2D-1P",
                                           "Stokes-2", "p", 0);

  TensorFieldFunction f2D_ten("src/Utility/Test/refdata/Field2D-1P",
                              "Stokes-1", "u_x,x|u_x,y|u_y,x|u_y,y", 0);

  VecGradFieldFunction f2Dg_vec("src/Utility/Test/refdata/Field2D-1P",
                               "Stokes-1", "u", 0);

  VecLaplacianFieldFunction f2Dl_vec("src/Utility/Test/refdata/Field2D-1P",
                                     "Stokes-1", "u", 0);

  STensorFieldFunction f2D_sten("src/Utility/Test/refdata/Field2D-1P",
                                "Stokes-1", "u_x,x|u_y,y|u_x,y", 0);

  double param[3] = {0.5, 0.5, 0.0};
  Vec4 X(param);
  double scal = f2D_scalar(X);
  double x = 0.5, y = 0.5;
  REQUIRE_THAT(scal, WithinRel(x*x*y*y, 1e-14));
  Vec3 sg = f2D_scalar.gradient(X);
  REQUIRE_THAT(sg[0], WithinRel(2*x*y*y, 1e-13));
  REQUIRE_THAT(sg[1], WithinRel(x*x*2*y, 1e-14));
  sg = f2Dg_scalar(X);
  REQUIRE_THAT(sg[0], WithinRel(2*x*y*y, 1e-13));
  REQUIRE_THAT(sg[1], WithinRel(x*x*2*y, 1e-14));

  Vec3 vec = f2D_vec(X);
  REQUIRE_THAT(vec[0], WithinRel( x*x*x*y*y/3.0, 1e-14));
  REQUIRE_THAT(vec[1], WithinRel(-x*x*y*y*y/3, 1e-14));

  Vec3 slap = f2Dl_scalar(X);
  REQUIRE_THAT(slap[0], WithinRel(2*y*y, 1e-13));
  REQUIRE_THAT(slap[1], WithinRel(2*x*x, 1e-14));

  Tensor vecg = f2D_vec.gradient(X);
  REQUIRE_THAT(vecg(1,1), WithinRel(x*x*y*y, 1e-14));
  REQUIRE_THAT(vecg(1,2), WithinRel(x*x*x*2.0*y/3.0, 1e-14));
  REQUIRE_THAT(vecg(2,1), WithinRel(-2.0*x*y*y*y/3.0, 1e-14));
  REQUIRE_THAT(vecg(2,2), WithinRel(-x*x*y*y, 1e-14));
  vecg = f2Dg_vec(X);
  REQUIRE_THAT(vecg(1,1), WithinRel(x*x*y*y, 1e-14));
  REQUIRE_THAT(vecg(1,2), WithinRel(x*x*x*2.0*y/3.0, 1e-14));
  REQUIRE_THAT(vecg(2,1), WithinRel(-2.0*x*y*y*y/3.0, 1e-14));
  REQUIRE_THAT(vecg(2,2), WithinRel(-x*x*y*y, 1e-14));

  Tensor ten = f2D_ten(X);
  REQUIRE_THAT(vecg(1,1), WithinRel(x*x*y*y, 1e-14));
  REQUIRE_THAT(vecg(1,2), WithinRel(x*x*x*2.0*y/3.0, 1e-14));
  REQUIRE_THAT(vecg(2,1), WithinRel(-2.0*x*y*y*y/3.0, 1e-14));
  REQUIRE_THAT(vecg(2,2), WithinRel(-x*x*y*y, 1e-14));

  Tensor lap = f2Dl_vec(X);
  REQUIRE_THAT(lap(1,1), WithinRel(2.0*x*y*y, 1e-13));
  REQUIRE_THAT(lap(1,2), WithinRel(x*x*x*2.0/3.0, 1e-13));
  REQUIRE_THAT(lap(2,1), WithinRel(-2.0*x*x*x/3.0, 1e-13));
  REQUIRE_THAT(lap(2,2), WithinRel(-x*x*2.0*y, 1e-14));

  SymmTensor sten = f2D_sten(X);
  REQUIRE_THAT(sten(1,1), WithinRel(x*x*y*y, 1e-14));
  REQUIRE_THAT(sten(1,2), WithinRel(x*x*x*2.0*y/3.0, 1e-14));
  REQUIRE_THAT(sten(2,1), WithinRel(x*x*x*2.0*y/3.0, 1e-14));
  REQUIRE_THAT(sten(2,2), WithinRel(-x*x*y*y, 1e-14));
}


TEST_CASE("TestFieldFunctions.2D1Pmx")
{
  // Fields contains: u_x = x^3*y^2 / 3
  //                  u_y = -x^2*y^3 / 3
  //                    p = x^2*y^2
  FieldFunction f2D_scalar("src/Utility/Test/refdata/Field2D-1Pmx",
                           "Stokes-3", "p", 0);

  ScalarGradFieldFunction f2Dg_scalar("src/Utility/Test/refdata/Field2D-1Pmx",
                                      "Stokes-3", "p", 0);

  VecFieldFunction f2D_vec("src/Utility/Test/refdata/Field2D-1Pmx",
                           "Stokes", "u_x|u_y", 0);

  VecGradFieldFunction f2Dg_vec("src/Utility/Test/refdata/Field2D-1Pmx",
                               "Stokes", "u_x|u_y", 0);

  VecLaplacianFieldFunction f2Dl_vec("src/Utility/Test/refdata/Field2D-1Pmx",
                                     "Stokes", "u_x|u_y", 0);

  double param[3] = {0.5, 0.5, 0.0};
  Vec4 X(param);
  double scal = f2D_scalar(X);
  double x = 0.5, y = 0.5;
  REQUIRE_THAT(scal, WithinRel(x*x*y*y, 1e-11));
  Vec3 sg = f2D_scalar.gradient(X);
  REQUIRE_THAT(sg[0], WithinRel(2*x*y*y, 1e-11));
  REQUIRE_THAT(sg[1], WithinRel(x*x*2*y, 1e-11));
  sg = f2Dg_scalar(X);
  REQUIRE_THAT(sg[0], WithinRel(2*x*y*y, 1e-11));
  REQUIRE_THAT(sg[1], WithinRel(x*x*2*y, 1e-11));

  Vec3 vec = f2D_vec(X);
  REQUIRE_THAT(vec[0], WithinRel( x*x*x*y*y/3.0, 1e-12));
  REQUIRE_THAT(vec[1], WithinRel(-x*x*y*y*y/3, 1e-13));

  Tensor vecg = f2D_vec.gradient(X);
  REQUIRE_THAT(vecg(1,1), WithinRel(x*x*y*y, 1e-11));
  REQUIRE_THAT(vecg(1,2), WithinRel(x*x*x*2.0*y/3.0, 1e-12));
  REQUIRE_THAT(vecg(2,1), WithinRel(-2.0*x*y*y*y/3.0, 1e-12));
  REQUIRE_THAT(vecg(2,2), WithinRel(-x*x*y*y, 1e-12));

  vecg = f2Dg_vec(X);
  REQUIRE_THAT(vecg(1,1), WithinRel(x*x*y*y, 1e-11));
  REQUIRE_THAT(vecg(1,2), WithinRel(x*x*x*2.0*y/3.0, 1e-12));
  REQUIRE_THAT(vecg(2,1), WithinRel(-2.0*x*y*y*y/3.0, 1e-12));
  REQUIRE_THAT(vecg(2,2), WithinRel(-x*x*y*y, 1e-12));

  Tensor lap = f2Dl_vec(X);
  REQUIRE_THAT(lap(1,1), WithinRel(2.0*x*y*y, 1e-12));
  REQUIRE_THAT(lap(1,2), WithinRel(x*x*x*2.0/3.0, 1e-11));
  REQUIRE_THAT(lap(2,1), WithinRel(-2.0*x*x*x/3.0, 1e-11));
  REQUIRE_THAT(lap(2,2), WithinRel(-x*x*2.0*y, 1e-11));
}


TEST_CASE("TestFieldFunctions.2D2P")
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
    REQUIRE_THAT(scal, WithinRel(0.5 + pIdx, 1e-14));
    f2D_vec.initPatch(pIdx);
    Vec3 vec = f2D_vec(X);
    REQUIRE_THAT(vec[0], WithinRel(0.5 + pIdx, 1e-14));
    REQUIRE_THAT(vec[1], WithinRel(-1.0, 1e-14));
  }
}


TEST_CASE("TestFieldFunctions.3D1P")
{
  // Fields contains: u_x = x^3*y^2*z^3 / 3
  //                  u_y = -x^2*y^3*z^3 / 3
  //                  u_z = -2*x^2*y^2*z^3 / 3
  //                    p = x^2*y^2*z^2
  FieldFunction f3D_scalar("src/Utility/Test/refdata/Field3D-1P",
                           "Stokes-2", "p", 0);

  ScalarGradFieldFunction f3Dg_scalar("src/Utility/Test/refdata/Field3D-1P",
                                      "Stokes-2", "p", 0);

  VecFieldFunction f3D_vec("src/Utility/Test/refdata/Field3D-1P",
                           "Stokes-1", "u", 0);

  ScalarLaplacianFieldFunction f3Dl_scalar("src/Utility/Test/refdata/Field3D-1P",
                                           "Stokes-2", "p", 0);

  TensorFieldFunction f3D_ten("src/Utility/Test/refdata/Field3D-1P",
                              "Stokes-1", "u_x,x|u_y,x|u_z,x|u_x,y|u_y,y|u_z,y|u_x,z|u_x,y|u_z,z", 0);

  VecGradFieldFunction f3Dg_vec("src/Utility/Test/refdata/Field3D-1P",
                                "Stokes-1", "u", 0);

  VecLaplacianFieldFunction f3Dl_vec("src/Utility/Test/refdata/Field3D-1P",
                                     "Stokes-1", "u", 0);

  STensorFieldFunction f3D_sten("src/Utility/Test/refdata/Field3D-1P",
                                "Stokes-1", "u_x,x|u_y,y|u_z,z|u_x,y|u_z,x|u_x,y", 0);

  double param[3] = {0.5, 0.5, 0.5};
  Vec4 X(param);
  double scal = f3D_scalar(X);
  double x = 0.5, y = 0.5, z = 0.5;
  REQUIRE_THAT(scal, WithinRel(x*x*y*y*z*z, 1e-13));
  Vec3 sg = f3D_scalar.gradient(X);
  REQUIRE_THAT(sg[0], WithinRel(2.0*x*y*y*z*z, 1e-12));
  REQUIRE_THAT(sg[1], WithinRel(x*x*2.0*y*z*z, 1e-12));
  REQUIRE_THAT(sg[2], WithinRel(x*x*y*y*2.0*z, 1e-12));

  sg = f3Dg_scalar(X);
  REQUIRE_THAT(sg[0], WithinRel(2.0*x*y*y*z*z, 1e-12));
  REQUIRE_THAT(sg[1], WithinRel(x*x*2.0*y*z*z, 1e-12));
  REQUIRE_THAT(sg[2], WithinRel(x*x*y*y*2.0*z, 1e-12));

  Vec3 vec = f3D_vec(X);
  REQUIRE_THAT(vec[0], WithinRel(x*x*x*y*y*z*z/3.0, 1e-13));
  REQUIRE_THAT(vec[1], WithinRel(x*x*y*y*y*z*z/3.0, 1e-13));
  REQUIRE_THAT(vec[2], WithinRel(-2.0*x*x*y*y*z*z*z/3.0, 1e-14));

  Tensor vecg = f3D_vec.gradient(X);
  REQUIRE_THAT(vecg(1,1), WithinRel(x*x*y*y*z*z, 1e-13));
  REQUIRE_THAT(vecg(1,2), WithinRel(x*x*x*2.0*y*z*z / 3.0, 1e-13));
  REQUIRE_THAT(vecg(1,3), WithinRel(x*x*x*y*y*2.0*z / 3.0, 1e-13));
  REQUIRE_THAT(vecg(2,1), WithinRel(2.0*x*y*y*y*z*z / 3.0, 1e-12));
  REQUIRE_THAT(vecg(2,2), WithinRel(x*x*y*y*z*z, 1e-13));
  REQUIRE_THAT(vecg(2,3), WithinRel(x*x*y*y*y*2.0*z / 3.0, 1e-13));
  REQUIRE_THAT(vecg(3,1), WithinRel(-2.0*2.0*x*y*y*z*z*z / 3.0, 1e-13));
  REQUIRE_THAT(vecg(3,2), WithinRel(-2.0*x*x*2.0*y*z*z*z / 3.0, 1e-13));
  REQUIRE_THAT(vecg(3,3), WithinRel(-2.0*x*x*y*y*z*z, 1e-14));

  Vec3 slap = f3Dl_scalar(X);
  REQUIRE_THAT(slap[0], WithinRel(2.0*y*y*z*z, 1e-12));
  REQUIRE_THAT(slap[1], WithinRel(x*x*2.0*z*z, 1e-12));
  REQUIRE_THAT(slap[2], WithinRel(x*x*y*y*2.0, 1e-12));

  vecg = f3Dg_vec(X);
  REQUIRE_THAT(vecg(1,1), WithinRel(x*x*y*y*z*z, 1e-13));
  REQUIRE_THAT(vecg(1,2), WithinRel(x*x*x*2.0*y*z*z / 3.0, 1e-13));
  REQUIRE_THAT(vecg(1,3), WithinRel(x*x*x*y*y*2.0*z / 3.0, 1e-13));
  REQUIRE_THAT(vecg(2,1), WithinRel(2.0*x*y*y*y*z*z / 3.0, 1e-12));
  REQUIRE_THAT(vecg(2,2), WithinRel(x*x*y*y*z*z, 1e-13));
  REQUIRE_THAT(vecg(2,3), WithinRel(x*x*y*y*y*2.0*z / 3.0, 1e-13));
  REQUIRE_THAT(vecg(3,1), WithinRel(-2.0*2.0*x*y*y*z*z*z / 3.0, 1e-13));
  REQUIRE_THAT(vecg(3,2), WithinRel(-2.0*x*x*2.0*y*z*z*z / 3.0, 1e-13));
  REQUIRE_THAT(vecg(3,3), WithinRel(-2.0*x*x*y*y*z*z, 1e-14));

  Tensor ten = f3D_ten(X);
  REQUIRE_THAT(ten(1,1), WithinRel(x*x*y*y*z*z, 1e-13));
  REQUIRE_THAT(ten(1,2), WithinRel(x*x*x*2.0*y*z*z / 3.0, 1e-13));
  REQUIRE_THAT(ten(1,3), WithinRel(x*x*x*y*y*2.0*z / 3.0, 1e-13));
  REQUIRE_THAT(ten(2,1), WithinRel(2*x*y*y*y*z*z / 3.0, 1e-12));
  REQUIRE_THAT(ten(2,2), WithinRel(x*x*y*y*z*z, 1e-13));
  REQUIRE_THAT(ten(2,3), WithinRel(x*x*y*y*y*2.0*z / 3.0, 1e-13));
  REQUIRE_THAT(ten(3,1), WithinRel(-2.0*2.0*x*y*y*z*z*z / 3.0, 1e-13));
  REQUIRE_THAT(ten(3,2), WithinRel(-2.0*x*x*2.0*y*z*z*z / 3.0, 1e-13));
  REQUIRE_THAT(ten(3,3), WithinRel(-2.0*x*x*y*y*z*z, 1e-14));

  Tensor lap = f3Dl_vec(X);
  REQUIRE_THAT(lap(1,1), WithinRel(2.0*x*y*y*z*z, 1e-12));
  REQUIRE_THAT(lap(1,2), WithinRel(x*x*x*2.0*z*z / 3.0, 1e-12));
  REQUIRE_THAT(lap(1,3), WithinRel(x*x*x*y*y*2.0 / 3.0, 1e-12));
  REQUIRE_THAT(lap(2,1), WithinRel(2.0*y*y*y*z*z / 3.0, 1e-12));
  REQUIRE_THAT(lap(2,2), WithinRel(x*x*2.0*y*z*z, 1e-12));
  REQUIRE_THAT(lap(2,3), WithinRel(x*x*y*y*y*2.0 / 3.0, 1e-12));
  REQUIRE_THAT(lap(3,1), WithinRel(-2.0*2.0*y*y*z*z*z / 3.0, 1e-12));
  REQUIRE_THAT(lap(3,2), WithinRel(-2.0*x*x*2.0*z*z*z / 3.0, 1e-12));
  REQUIRE_THAT(lap(3,3), WithinRel(-2.0*x*x*y*y*2.0*z, 1e-13));

  SymmTensor sten = f3D_sten(X);
  REQUIRE_THAT(sten(1,1), WithinRel(x*x*y*y*z*z, 1e-13));
  REQUIRE_THAT(sten(1,2), WithinRel(x*x*x*2.0*y*z*z / 3.0, 1e-13));
  REQUIRE_THAT(sten(1,3), WithinRel(x*x*x*y*y*2.0*z / 3.0, 1e-13));
  REQUIRE_THAT(sten(2,1), WithinRel(2.0*x*y*y*y*z*z / 3.0, 1e-13));
  REQUIRE_THAT(sten(2,2), WithinRel(x*x*y*y*z*z, 1e-13));
  REQUIRE_THAT(sten(2,3), WithinRel(-2.0*x*x*2.0*y*z*z*z / 3.0, 1e-13));
  REQUIRE_THAT(sten(3,1), WithinRel(x*x*x*y*y*2.0*z / 3.0, 1e-13));
  REQUIRE_THAT(sten(3,2), WithinRel(-2.0*x*x*2.0*y*z*z*z / 3.0, 1e-13));
  REQUIRE_THAT(sten(3,3), WithinRel(-2.0*x*x*y*y*z*z, 1e-14));
}


TEST_CASE("TestFieldFunctions.3D1Pmx")
{
  // Fields contains: u_x = x^3*y^2*z^3 / 3
  //                  u_y = -x^2*y^3*z^3 / 3
  //                  u_z = -2*x^2*y^2*z^3 / 3
  //                    p = x^2*y^2*z^2
  FieldFunction f3D_scalar("src/Utility/Test/refdata/Field3D-1Pmx",
                           "Stokes-4", "p", 0);

  ScalarGradFieldFunction f3Dg_scalar("src/Utility/Test/refdata/Field3D-1Pmx",
                                      "Stokes-4", "p", 0);

  VecFieldFunction f3D_vec("src/Utility/Test/refdata/Field3D-1Pmx",
                           "Stokes", "u_x|u_y|u_z", 0);

  ScalarLaplacianFieldFunction f3Dl_scalar("src/Utility/Test/refdata/Field3D-1Pmx",
                                           "Stokes-4", "p", 0);

  VecGradFieldFunction f3Dg_vec("src/Utility/Test/refdata/Field3D-1Pmx",
                                "Stokes", "u_x|u_y|u_z", 0);

  VecLaplacianFieldFunction f3Dl_vec("src/Utility/Test/refdata/Field3D-1Pmx",
                                     "Stokes", "u_x|u_y|u_z", 0);

  double param[3] = {0.5, 0.5, 0.5};
  Vec4 X(param);
  double scal = f3D_scalar(X);
  double x = 0.5, y = 0.5, z = 0.5;
  REQUIRE_THAT(scal, WithinRel(x*x*y*y*z*z, 1e-14));

  Vec3 sg = f3D_scalar.gradient(X);
  REQUIRE_THAT(sg[0], WithinRel(2.0*x*y*y*z*z, 1e-13));
  REQUIRE_THAT(sg[1], WithinRel(x*x*2.0*y*z*z, 1e-13));
  REQUIRE_THAT(sg[2], WithinRel(x*x*y*y*2.0*z, 1e-13));

  sg = f3Dg_scalar(X);
  REQUIRE_THAT(sg[0], WithinRel(2.0*x*y*y*z*z, 1e-13));
  REQUIRE_THAT(sg[1], WithinRel(x*x*2.0*y*z*z, 1e-13));
  REQUIRE_THAT(sg[2], WithinRel(x*x*y*y*2.0*z, 1e-13));

  Vec3 vec = f3D_vec(X);
  REQUIRE_THAT(vec[0], WithinRel(x*x*x*y*y*z*z/3.0, 1e-14));
  REQUIRE_THAT(vec[1], WithinRel(x*x*y*y*y*z*z/3.0, 1e-14));
  REQUIRE_THAT(vec[2], WithinRel(-2.0*x*x*y*y*z*z*z/3.0, 1e-14));

  Tensor vecg = f3D_vec.gradient(X);
  REQUIRE_THAT(vecg(1,1), WithinRel(x*x*y*y*z*z, 1e-14));
  REQUIRE_THAT(vecg(1,2), WithinRel(x*x*x*2.0*y*z*z / 3.0, 1e-14));
  REQUIRE_THAT(vecg(1,3), WithinRel(x*x*x*y*y*2.0*z / 3.0, 1e-13));
  REQUIRE_THAT(vecg(2,1), WithinRel(2.0*x*y*y*y*z*z / 3.0, 1e-13));
  REQUIRE_THAT(vecg(2,2), WithinRel(x*x*y*y*z*z, 1e-14));
  REQUIRE_THAT(vecg(2,3), WithinRel(x*x*y*y*y*2.0*z / 3.0, 1e-13));
  REQUIRE_THAT(vecg(3,1), WithinRel(-2.0*2.0*x*y*y*z*z*z / 3.0, 1e-14));
  REQUIRE_THAT(vecg(3,2), WithinRel(-2.0*x*x*2.0*y*z*z*z / 3.0, 1e-14));
  REQUIRE_THAT(vecg(3,3), WithinRel(-2.0*x*x*y*y*z*z, 1e-14));

  Vec3 slap = f3Dl_scalar(X);
  REQUIRE_THAT(slap[0], WithinRel(2.0*y*y*z*z, 1e-13));
  REQUIRE_THAT(slap[1], WithinRel(x*x*2.0*z*z, 1e-13));
  REQUIRE_THAT(slap[2], WithinRel(x*x*y*y*2.0, 1e-13));

  vecg = f3Dg_vec(X);
  REQUIRE_THAT(vecg(1,1), WithinRel(x*x*y*y*z*z, 1e-14));
  REQUIRE_THAT(vecg(1,2), WithinRel(x*x*x*2.0*y*z*z / 3.0, 1e-14));
  REQUIRE_THAT(vecg(1,3), WithinRel(x*x*x*y*y*2.0*z / 3.0, 1e-13));
  REQUIRE_THAT(vecg(2,1), WithinRel(2.0*x*y*y*y*z*z / 3.0, 1e-13));
  REQUIRE_THAT(vecg(2,2), WithinRel(x*x*y*y*z*z, 1e-14));
  REQUIRE_THAT(vecg(2,3), WithinRel(x*x*y*y*y*2.0*z / 3.0, 1e-13));
  REQUIRE_THAT(vecg(3,1), WithinRel(-2.0*2.0*x*y*y*z*z*z / 3.0, 1e-14));
  REQUIRE_THAT(vecg(3,2), WithinRel(-2.0*x*x*2.0*y*z*z*z / 3.0, 1e-14));
  REQUIRE_THAT(vecg(3,3), WithinRel(-2.0*x*x*y*y*z*z, 1e-14));

  Tensor lap = f3Dl_vec(X);
  REQUIRE_THAT(lap(1,1), WithinRel(2.0*x*y*y*z*z, 1e-13));
  REQUIRE_THAT(lap(1,2), WithinRel(x*x*x*2.0*z*z / 3.0, 1e-13));
  REQUIRE_THAT(lap(1,3), WithinRel(x*x*x*y*y*2.0 / 3.0, 1e-13));
  REQUIRE_THAT(lap(2,1), WithinRel(2.0*y*y*y*z*z / 3.0, 1e-13));
  REQUIRE_THAT(lap(2,2), WithinRel(x*x*2.0*y*z*z, 1e-13));
  REQUIRE_THAT(lap(2,3), WithinRel(x*x*y*y*y*2.0 / 3.0, 1e-13));
  REQUIRE_THAT(lap(3,1), WithinRel(-2.0*2.0*y*y*z*z*z / 3.0, 1e-13));
  REQUIRE_THAT(lap(3,2), WithinRel(-2.0*x*x*2.0*z*z*z / 3.0, 1e-13));
  REQUIRE_THAT(lap(3,3), WithinRel(-2.0*x*x*y*y*2.0*z, 1e-13));
}


TEST_CASE("TestFieldFunctions.3D2P")
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
    REQUIRE_THAT(scal, WithinRel(1.0-2*x, 1e-14));
    f3D_vec.initPatch(pIdx);
    Vec3 vec = f3D_vec(X);
    REQUIRE_THAT(vec[0], WithinRel(y*(1.0-y), 1e-14)); // y*(1-y)
    REQUIRE_THAT(vec[1], WithinRel(z*(1.0-z), 1e-14)); // z*(1-z)
    REQUIRE_THAT(vec[2], WithinRel(x*(1.0-x), 1e-14)); // x*(1-x)
    f3D_ten.initPatch(pIdx);
    Tensor ten = f3D_ten(X);
    REQUIRE_THAT(ten(1,1), WithinAbs(0.0, 1e-14));
    REQUIRE_THAT(ten(1,2), WithinAbs(1.0-2*y, 1e-12));
    REQUIRE_THAT(ten(1,3), WithinAbs(0.0, 1e-14));
    REQUIRE_THAT(ten(2,1), WithinAbs(0.0, 1e-14));
    REQUIRE_THAT(ten(2,2), WithinAbs(0.0, 1e-14));
    REQUIRE_THAT(ten(2,3), WithinAbs(1.0-2*z, 1e-14));
    REQUIRE_THAT(ten(3,1), WithinRel(1.0-2*x, 1e-14));
    REQUIRE_THAT(ten(3,2), WithinAbs(0.0, 1e-14));
    REQUIRE_THAT(ten(3,3), WithinAbs(0.0, 1e-14));

    f3D_sten.initPatch(pIdx);
    SymmTensor sten = f3D_sten(X);
    REQUIRE_THAT(sten(1,1), WithinAbs(0.0, 1e-14));
    REQUIRE_THAT(sten(1,2), WithinAbs(1.0-2*y, 1e-14));
    REQUIRE_THAT(sten(1,3), WithinRel(1.0-2*x, 1e-14));
    REQUIRE_THAT(sten(2,1), WithinAbs(1.0-2*y, 1e-14));
    REQUIRE_THAT(sten(2,2), WithinAbs(0.0, 1e-14));
    REQUIRE_THAT(sten(2,3), WithinAbs(1.0-2*z, 1e-14));
    REQUIRE_THAT(sten(3,1), WithinRel(1.0-2*x, 1e-14));
    REQUIRE_THAT(sten(3,2), WithinAbs(1.0-2*z, 1e-14));
    REQUIRE_THAT(sten(3,3), WithinAbs(0.0, 1e-14));
  }
}
