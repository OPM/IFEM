//==============================================================================
//!
//! \file TestFuncAutoDiff.C
//!
//! \date Sep 20 2023
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for expression functions with auto-differentation.
//!
//==============================================================================

#include "ExprFunctions.h"
#include "TensorFunction.h"
#include "Functions.h"

#include <autodiff/reverse/var.hpp>

#include <cstdlib>
#include <cmath>

#include "Catch2Support.h"


using RealFuncAd    = EvalFuncSpatial<autodiff::var>;
using VecFuncAd     = EvalMultiFunction<VecFunc,Vec3,autodiff::var>;
using TensorFuncAd  = EvalMultiFunction<TensorFunc,Tensor,autodiff::var>;
using STensorFuncAd = EvalMultiFunction<STensorFunc,SymmTensor,autodiff::var>;


TEST_CASE("TestScalarFunc.AutoDiff")
{
  EvalFuncScalar<autodiff::var> f("sin(1.5*t)*t","t");

  CHECK(!f.isConstant());

  double t = 0.0;
  for (int i = 0; i < 20; i++)
  {
    t += 0.314*(double)random()/(double)RAND_MAX;
    CHECK_THAT(f.deriv(t), WithinRel(1.5*cos(1.5*t)*t+sin(1.5*t)));
  }
}


TEST_CASE("TestRealFunc.Gradient")
{
  const char* g   = "sin(x)*sin(y)*sin(z)";
  const char* g_x = "cos(x)*sin(y)*sin(z)";
  const char* g_y = "sin(x)*cos(y)*sin(z)";
  const char* g_z = "sin(x)*sin(y)*cos(z)";

  EvalFunction f1(g);
  f1.addDerivative(g_x, "", 1);
  f1.addDerivative(g_y, "", 2);
  f1.addDerivative(g_z, "", 3);

  RealFuncAd f2(g);

  CHECK(f1.isConstant());
  CHECK(f2.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const Vec3 r(cos(x)*sin(y)*sin(z),
                     sin(x)*cos(y)*sin(z),
                     sin(x)*sin(y)*cos(z));

        for (const RealFunc* fp : {static_cast<const RealFunc*>(&f1),
                                   static_cast<const RealFunc*>(&f2)}) {
          const Vec3 grad = fp->gradient(X);
          for (size_t i = 1; i <= 3; ++i) {
            CHECK_THAT(fp->deriv(X, i), WithinRel(r[i-1]));
            CHECK_THAT(grad[i-1], WithinRel(r[i-1]));
          }
        }
      }
}


TEST_CASE("TestRealFunc.Hessian")
{
  const char* g    = "sin(x)*sin(y)*sin(z)";
  const char* g_xx = "-sin(x)*sin(y)*sin(z)";
  const char* g_xy = "cos(x)*cos(y)*sin(z)";
  const char* g_xz = "cos(x)*sin(y)*cos(z)";
  const char* g_yy = "-sin(x)*sin(y)*sin(z)";
  const char* g_zz = "-sin(x)*sin(y)*sin(z)";
  const char* g_yz = "sin(x)*cos(y)*cos(z)";

  EvalFunction f1(g);
  f1.addDerivative(g_xx,"",1,1);
  f1.addDerivative(g_xy,"",1,2);
  f1.addDerivative(g_xz,"",1,3);
  f1.addDerivative(g_yy,"",2,2);
  f1.addDerivative(g_yz,"",2,3);
  f1.addDerivative(g_zz,"",3,3);

  RealFuncAd f2(g);

  CHECK(f1.isConstant());
  CHECK(f2.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const SymmTensor r({-sin(x)*sin(y)*sin(z),
                            -sin(x)*sin(y)*sin(z),
                            -sin(x)*sin(y)*sin(z),
                             cos(x)*cos(y)*sin(z),
                             sin(x)*cos(y)*cos(z),
                             cos(x)*sin(y)*cos(z)});

        for (const RealFunc* fp : {static_cast<const RealFunc*>(&f1),
                                   static_cast<const RealFunc*>(&f2)}) {
          const SymmTensor hess = fp->hessian(X);
          for (size_t i = 1; i <= 3; ++i)
            for (size_t j = 1; j <= 3; ++j) {
              CHECK_THAT(fp->dderiv(X,i,j), WithinRel(r(i,j)));
              CHECK_THAT(hess(i,j), WithinRel(r(i,j)));
            }
        }
      }
}


TEST_CASE("TestVecFunc.Evaluate")
{
  const char* func = "sin(x) | cos (y) | exp(z)";

  VecFuncExpr f1(func);
  VecFuncAd f2(func);

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const Vec3 r(sin(x), cos(y), exp(z));

        for (const VecFunc* fp : {static_cast<const VecFunc*>(&f1),
                                  static_cast<const VecFunc*>(&f2)}) {
          const Vec3 fx = (*fp)(X);
          CHECK_THAT(fx.x, WithinRel(r.x));
          CHECK_THAT(fx.y, WithinRel(r.y));
          CHECK_THAT(fx.z, WithinRel(r.z));
        }
      }
}


TEST_CASE("TestVecFunction.Gradient2D")
{
  const char* g   = "sin(x)*sin(y) | x*x*y*y";
  const char* g_x = "cos(x)*sin(y) | 2*x*y*y";
  const char* g_y = "sin(x)*cos(y) | 2*x*x*y";

  VecFuncExpr f1(g);
  f1.addDerivative(g_x,"",1);
  f1.addDerivative(g_y,"",2);

  VecFuncAd f2(g);

  CHECK(f1.isConstant());
  CHECK(f2.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7}) {
      const Vec3 X(x,y);
      const Tensor r({cos(x)*sin(y), 2*x*y*y,
                      sin(x)*cos(y), 2*x*x*y});

      for (const VecFunc* fp : {static_cast<const VecFunc*>(&f1),
                                static_cast<const VecFunc*>(&f2)}) {
        const Tensor grad = fp->gradient(X);
        for (size_t d = 1; d <= 2; ++d) {
          const Vec3 dx = fp->deriv(X,d);
          for (size_t i = 1; i <= 2; ++i) {
            CHECK_THAT(dx[i-1], WithinRel(r(i,d)));
            CHECK_THAT(grad(i,d), WithinRel(r(i,d)));
          }
        }
      }
    }
}


TEST_CASE("TestVecFunction.Gradient3D")
{
  const char* g   = "sin(x)*sin(y)*sin(z) | x*x*y*y*z*z*z | exp(-x)*exp(2*y)*exp(z*z)";
  const char* g_x = "cos(x)*sin(y)*sin(z) | 2*x*y*y*z*z*z | -exp(-x)*exp(2*y)*exp(z*z)";
  const char* g_y = "sin(x)*cos(y)*sin(z) | 2*x*x*y*z*z*z | 2*exp(-x)*exp(2*y)*exp(z*z)";
  const char* g_z = "sin(x)*sin(y)*cos(z) | 3*x*x*y*y*z*z | 2*z*exp(-x)*exp(2*y)*exp(z*z)";

  VecFuncExpr f1(g);
  f1.addDerivative(g_x,"",1);
  f1.addDerivative(g_y,"",2);
  f1.addDerivative(g_z,"",3);

  VecFuncAd f2(g);

  CHECK(f1.isConstant());
  CHECK(f2.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const Tensor r({cos(x)*sin(y)*sin(z), 2*x*y*y*z*z*z, -exp(-x)*exp(2*y)*exp(z*z),
                        sin(x)*cos(y)*sin(z), 2*x*x*y*z*z*z, 2*exp(-x)*exp(2*y)*exp(z*z),
                        sin(x)*sin(y)*cos(z), 3*x*x*y*y*z*z, 2*z*exp(-x)*exp(2*y)*exp(z*z)});

        for (const VecFunc* fp : {static_cast<const VecFunc*>(&f1),
                                  static_cast<const VecFunc*>(&f2)}) {
          const Tensor grad = fp->gradient(X);
          for (size_t d = 1; d <= 3; ++d) {
            const Vec3 dx = fp->deriv(X,d);
            for (size_t i = 1; i <= 3; ++i) {
              CHECK_THAT(dx[i-1], WithinRel(r(i,d)));
              CHECK_THAT(grad(i,d), WithinRel(r(i,d)));
            }
          }
        }
      }
}


TEST_CASE("TestVecFunction.Hessian2D")
{
  const char* g    = "sin(x)*sin(y)  | x*x*y*y";
  const char* g_xx = "-sin(x)*sin(y) | 2*y*y";
  const char* g_xy = "cos(x)*cos(y)  | 2*x*2*y";
  const char* g_yy = "-sin(x)*sin(y) | 2*x*x";

  VecFuncExpr f1(g);
  f1.addDerivative(g_xx,"",1,1);
  f1.addDerivative(g_yy,"",2,2);
  f1.addDerivative(g_xy,"",1,2);

  VecFuncAd f2(g);

  CHECK(f1.isConstant());
  CHECK(f2.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7}) {
      const Vec3 X(x,y);
      utl::matrix3d<Real> r(2,2,2);
      r.fill(std::array{-sin(x)*sin(y), 2*y*y,
                         cos(x)*cos(y), 2*x*2*y,
                         cos(x)*cos(y), 2*x*2*y,
                        -sin(x)*sin(y), 2*x*x}.data());

      for (const VecFunc* fp : {static_cast<const VecFunc*>(&f1),
                                static_cast<const VecFunc*>(&f2)}) {
        const utl::matrix3d<Real> hess = fp->hessian(X);
        for (size_t d1 = 1; d1 <= 2; ++d1)
          for (size_t d2 = 1; d2 <= 2; ++d2) {
            const Vec3 d2x = fp->dderiv(X,d1,d2);
            for (size_t i = 1; i <= 2; ++i) {
              CHECK_THAT(hess(i,d1,d2), WithinRel(r(i,d1,d2)));
              CHECK_THAT(d2x[i-1], WithinRel(r(i,d1,d2)));
            }
          }
      }
    }
}


TEST_CASE("TestVecFunction.Hessian3D")
{
  const char* g    = "sin(x)*sin(y)*sin(z)  | x*x*y*y*z*z | exp(x)*exp(2*y)*exp(3*z)";
  const char* g_xx = "-sin(x)*sin(y)*sin(z) | 2*y*y*z*z   | exp(x)*exp(2*y)*exp(3*z)";
  const char* g_xy = "cos(x)*cos(y)*sin(z)  | 2*x*2*y*z*z | exp(x)*2*exp(2*y)*exp(3*z)";
  const char* g_xz = "cos(x)*sin(y)*cos(z)  | 2*x*y*y*2*z | exp(x)*exp(2*y)*3*exp(3*z)";
  const char* g_yy = "-sin(x)*sin(y)*sin(z) | x*x*2*z*z   | exp(x)*4*exp(2*y)*exp(3*z)";
  const char* g_yz = "sin(x)*cos(y)*cos(z)  | x*x*2*y*2*z | exp(x)*2*exp(2*y)*3*exp(3*z)";
  const char* g_zz = "-sin(x)*sin(y)*sin(z) | x*x*y*y*2   | exp(x)*exp(2*y)*9*exp(3*z)";

  VecFuncExpr f1(g);
  f1.addDerivative(g_xx,"",1,1);
  f1.addDerivative(g_yy,"",2,2);
  f1.addDerivative(g_zz,"",3,3);
  f1.addDerivative(g_xy,"",1,2);
  f1.addDerivative(g_xz,"",1,3);
  f1.addDerivative(g_yz,"",2,3);

  VecFuncAd f2(g);

  CHECK(f1.isConstant());
  CHECK(f2.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        utl::matrix3d<Real> r(3,3,3);
        r.fill(std::array{-sin(x)*sin(y)*sin(z), 2*y*y*z*z,   exp(x)*exp(2*y)*exp(3*z),
                           cos(x)*cos(y)*sin(z), 2*x*2*y*z*z, exp(x)*2*exp(2*y)*exp(3*z),
                           cos(x)*sin(y)*cos(z), 2*x*y*y*2*z, exp(x)*exp(2*y)*3*exp(3*z),

                           cos(x)*cos(y)*sin(z), 2*x*2*y*z*z, exp(x)*2*exp(2*y)*exp(3*z),
                          -sin(x)*sin(y)*sin(z), x*x*2*z*z,   exp(x)*4*exp(2*y)*exp(3*z),
                           sin(x)*cos(y)*cos(z), x*x*2*y*2*z, exp(x)*2*exp(2*y)*3*exp(3*z),

                           cos(x)*sin(y)*cos(z), 2*x*y*y*2*z, exp(x)*exp(2*y)*3*exp(3*z),
                           sin(x)*cos(y)*cos(z), x*x*2*y*2*z, exp(x)*2*exp(2*y)*3*exp(3*z),
                          -sin(x)*sin(y)*sin(z), x*x*y*y*2,   exp(x)*exp(2*y)*9*exp(3*z)}.data());

        for (const VecFunc* fp : {static_cast<const VecFunc*>(&f1),
                                  static_cast<const VecFunc*>(&f2)}) {
          const utl::matrix3d<Real> hess = fp->hessian(X);
          for (size_t d1 = 1; d1 <= 3; ++d1)
            for (size_t d2 = 1; d2 <= 3; ++d2) {
              const Vec3 d2x = fp->dderiv(X,d1,d2);
              for (size_t i = 1; i <= 3; ++i) {
                CHECK_THAT(d2x[i-1], WithinRel(r(i,d1,d2)));
                CHECK_THAT(hess(i,d1,d2), WithinRel(r(i,d1,d2)));
              }
            }
        }
      }
}


TEST_CASE("TestTensorFunc.Evaluate")
{
  const char* func = "sin(x) | cos (y) | exp(z) | sin(x)*cos(y)";

  TensorFuncExpr f1(func);
  TensorFuncAd f2(func);

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const Tensor r({sin(x), cos(y), exp(z), sin(x)*cos(y)});

        for (const TensorFunc* fp : {static_cast<const TensorFunc*>(&f1),
                                     static_cast<const TensorFunc*>(&f2)}) {
          const Tensor fx = (*fp)(X);
          for (size_t i = 1; i <= 2; ++i)
            for (size_t j = 1; j <= 2; ++j)
              CHECK_THAT(fx(i,j), WithinRel(r(i,j)));
        }
      }
}


TEST_CASE("TestTensorFunction.Gradient2D")
{
  const char* g   = "sin(x)*sin(y) | x*x*y*y | exp(x)*exp(2*y)   | exp(-2*x)*exp(y)";
  const char* g_x = "cos(x)*sin(y) | 2*x*y*y | exp(x)*exp(2*y)   | -2*exp(-2*x)*exp(y)";
  const char* g_y = "sin(x)*cos(y) | 2*x*x*y | 2*exp(x)*exp(2*y) | exp(-2*x)*exp(y)";

  TensorFuncExpr f1(g);
  f1.addDerivative(g_x,"",1);
  f1.addDerivative(g_y,"",2);

  TensorFuncAd f2(g);

  CHECK(f1.isConstant());
  CHECK(f2.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7}) {
      const Vec3 X(x,y);
      utl::matrix3d<Real> r(2,2,2);
      r.fill(std::array{cos(x)*sin(y), 2*x*y*y, exp(x)*exp(2*y), -2*exp(-2*x)*exp(y),
                        sin(x)*cos(y), 2*x*x*y, 2*exp(x)*exp(2*y), exp(-2*x)*exp(y)}.data());

      for (const TensorFunc* fp : {static_cast<const TensorFunc*>(&f1),
                                   static_cast<const TensorFunc*>(&f2)}) {
        const utl::matrix3d<Real> grad = fp->gradient(X);
        for (size_t d = 1; d <= 2; ++d) {
          const Tensor dx = fp->deriv(X,d);
          for (size_t i = 1; i <= 2; ++i)
            for (size_t j = 1; j <= 2; ++j) {
              CHECK_THAT(dx(i,j), WithinRel(r(i,j,d)));
              CHECK_THAT(grad(i,j,d), WithinRel(r(i,j,d)));
            }
        }
      }
    }
}


TEST_CASE("TestTensorFunction.Gradient3D")
{
  const char* g   = "sin(x)*sin(y)*sin(z)  | x*x*y*y*z | exp(x)*exp(2*y)*z*z |"
                    "exp(-2*x)*exp(y)*z    | x*y*z     | x*y*z*z |"
                    "x                     | y         | z";
  const char* g_x = "cos(x)*sin(y)*sin(z)  | 2*x*y*y*z | exp(x)*exp(2*y)*z*z |"
                    "-2*exp(-2*x)*exp(y)*z | y*z       | y*z*z |"
                    "1.0                   | 0.0       | 0.0";
  const char* g_y = "sin(x)*cos(y)*sin(z)  | 2*x*x*y*z | 2*exp(x)*exp(2*y)*z*z |"
                    "exp(-2*x)*exp(y)*z    | x*z       | x*z*z |"
                    "0.0                   | 1.0       | 0.0";
  const char* g_z = "sin(x)*sin(y)*cos(z)  | x*x*y*y   | exp(x)*exp(2*y)*2*z |"
                    "exp(-2*x)*exp(y)      | x*y       | x*y*2*z |"
                    "0.0                   | 0.0       | 1.0";

  TensorFuncExpr f1(g);
  f1.addDerivative(g_x,"",1);
  f1.addDerivative(g_y,"",2);
  f1.addDerivative(g_z,"",3);

  TensorFuncAd f2(g);

  CHECK(f1.isConstant());
  CHECK(f2.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        utl::matrix3d<Real> r(3,3,3);
        r.fill(std::array{cos(x)*sin(y)*sin(z),  2*x*y*y*z, exp(x)*exp(2*y)*z*z,
                          -2*exp(-2*x)*exp(y)*z, y*z,       y*z*z,
                          1.0,                   0.0,       0.0,

                          sin(x)*cos(y)*sin(z),  2*x*x*y*z, 2*exp(x)*exp(2*y)*z*z,
                          exp(-2*x)*exp(y)*z,    x*z,       x*z*z,
                          0.0,                   1.0,       0.0,

                          sin(x)*sin(y)*cos(z),  x*x*y*y,   exp(x)*exp(2*y)*2*z,
                          exp(-2*x)*exp(y),      x*y,       x*y*2*z,
                          0.0,                   0.0,       1.0}.data());

        for (const TensorFunc* fp : {static_cast<const TensorFunc*>(&f1),
                                     static_cast<const TensorFunc*>(&f2)}) {
          const utl::matrix3d<Real> grad = fp->gradient(X);
          for (size_t d = 1; d <= 3; ++d) {
            const Tensor dx = fp->deriv(X,d);
            for (size_t i = 1; i <= 3; ++i)
              for (size_t j = 1; j <= 3; ++j) {
                CHECK_THAT(dx(i,j), WithinRel(r(i,j,d)));
                CHECK_THAT(grad(i,j,d), WithinRel(r(i,j,d)));
              }
          }
        }
      }
}


TEST_CASE("TestTensorFunction.Hessian2D")
{
  const char* g    = "sin(x)*sin(y)  | x*x*y*y | exp(x)*exp(2*y)   | exp(-2*x)*exp(y)";
  const char* g_xx = "-sin(x)*sin(y) | 2*y*y   | exp(x)*exp(2*y)   | 4*exp(-2*x)*exp(y)";
  const char* g_xy = "cos(x)*cos(y)  | 2*x*2*y | exp(x)*2*exp(2*y) | -2*exp(-2*x)*exp(y)";
  const char* g_yy = "-sin(x)*sin(y) | 2*x*x   | exp(x)*4*exp(2*y) | exp(-2*x)*exp(y)";

  TensorFuncExpr f1(g);
  f1.addDerivative(g_xx,"",1,1);
  f1.addDerivative(g_yy,"",2,2);
  f1.addDerivative(g_xy,"",1,2);

  TensorFuncAd f2(g);

  CHECK(f1.isConstant());
  CHECK(f2.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7}) {
      const Vec3 X(x,y);
      utl::matrix4d<Real> r(2,2,2,2);
      r.fill(std::array{-sin(x)*sin(y), 2.0*y*y,
                        exp(x)*exp(2*y), 4.0*exp(-2*x)*exp(y),

                        cos(x)*cos(y), 2*x*2*y,
                        exp(x)*2*exp(2*y), -2*exp(-2*x)*exp(y),

                        cos(x)*cos(y), 2*x*2*y,
                        exp(x)*2*exp(2*y), -2*exp(-2*x)*exp(y),

                        -sin(x)*sin(y), 2*x*x,
                        exp(x)*4*exp(2*y), exp(-2*x)*exp(y)}.data());

      for (const TensorFunc* fp : {static_cast<const TensorFunc*>(&f1),
                                   static_cast<const TensorFunc*>(&f2)}) {
        const utl::matrix4d<Real> hess = fp->hessian(X);
        for (size_t d1 = 1; d1 <= 2; ++d1)
          for (size_t d2 = 1; d2 <= 2; ++d2) {
            const Tensor dx = fp->dderiv(X,d1,d2);
            for (size_t i = 1; i <= 2; ++i)
              for (size_t j = 1; j <= 2; ++j) {
                CHECK_THAT(dx(i,j), WithinRel(r(i,j,d1,d2)));
                CHECK_THAT(hess(i,j,d1,d2), WithinRel(r(i,j,d1,d2)));
              }
          }
      }
    }
}


TEST_CASE("TestTensorFunction.Hessian3D")
{
  const char* g    = "sin(x)*sin(y)*sin(z)  | x*x*y*y*z | exp(x)*exp(2*y)*z*z |"
                     "exp(-2*x)*exp(y)*z    | x*y*z     | x*y*z*z |"
                     "x*x                   | y*y       | z*z";
  const char* g_xx = "-sin(x)*sin(y)*sin(z) | 2*y*y*z   | exp(x)*exp(2*y)*z*z |"
                     "4*exp(-2*x)*exp(y)*z  | 0.0       | 0.0 |"
                     "2.0                   | 0.0       | 0.0";
  const char* g_xy = "cos(x)*cos(y)*sin(z)  | 2*x*2*y*z | exp(x)*2*exp(2*y)*z*z |"
                     "-2*exp(-2*x)*exp(y)*z | z         | z*z |"
                     "0.0                   | 0.0       | 0.0";
  const char* g_xz = "cos(x)*sin(y)*cos(z)  | 2*x*y*y   | exp(x)*exp(2*y)*2*z |"
                     "-2*exp(-2*x)*exp(y)   | y         | y*2*z |"
                     "0.0                   | 0.0       | 0.0";
  const char* g_yy = "-sin(x)*sin(y)*sin(z) | x*x*2*z   | exp(x)*4*exp(2*y)*z*z |"
                     "exp(-2*x)*exp(y)*z    | 0.0       | 0.0 |"
                     "0.0                   | 2.0       | 0.0";
  const char* g_yz = "sin(x)*cos(y)*cos(z)  | x*x*2*y   | exp(x)*2*exp(2*y)*2*z |"
                     "exp(-2*x)*exp(y)      | x         | x*2*z |"
                     "0.0                   | 0.0       | 0.0";
  const char* g_zz = "-sin(x)*sin(y)*sin(z) | 0.0       | exp(x)*exp(2*y)*2 |"
                     "0.0                   | 0.0       | 2.0*x*y |"
                     "0.0                   | 0.0       | 2.0";

  TensorFuncExpr f1(g);
  f1.addDerivative(g_xx,"",1,1);
  f1.addDerivative(g_yy,"",2,2);
  f1.addDerivative(g_zz,"",3,3);
  f1.addDerivative(g_xy,"",1,2);
  f1.addDerivative(g_xz,"",1,3);
  f1.addDerivative(g_yz,"",2,3);

  TensorFuncAd f2(g);

  CHECK(f1.isConstant());
  CHECK(f2.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        utl::matrix4d<Real> r(3,3,3,3);
        r.fill(std::array{-sin(x)*sin(y)*sin(z), 2*y*y*z,   exp(x)*exp(2*y)*z*z,
                          4*exp(-2*x)*exp(y)*z,  0.0,       0.0,
                          2.0,                   0.0,       0.0,

                          cos(x)*cos(y)*sin(z),  2*x*2*y*z, exp(x)*2*exp(2*y)*z*z,
                          -2*exp(-2*x)*exp(y)*z, z,         z*z,
                          0.0,                   0.0,       0.0,

                          cos(x)*sin(y)*cos(z),  2*x*y*y,   exp(x)*exp(2*y)*2*z,
                          -2*exp(-2*x)*exp(y),   y,         y*2*z,
                          0.0,                   0.0,       0.0,

                          cos(x)*cos(y)*sin(z),  2*x*2*y*z, exp(x)*2*exp(2*y)*z*z,
                          -2*exp(-2*x)*exp(y)*z, z,         z*z,
                          0.0,                   0.0,       0.0,

                          -sin(x)*sin(y)*sin(z), x*x*2*z,   exp(x)*4*exp(2*y)*z*z,
                          exp(-2*x)*exp(y)*z,    0.0,       0.0,
                          0.0,                   2.0,       0.0,

                          sin(x)*cos(y)*cos(z),  x*x*2*y,   exp(x)*2*exp(2*y)*2*z,
                          exp(-2*x)*exp(y),      x,         x*2*z,
                          0.0,                   0.0,       0.0,

                          cos(x)*sin(y)*cos(z),  2*x*y*y,   exp(x)*exp(2*y)*2*z,
                          -2*exp(-2*x)*exp(y),   y,         y*2*z,
                          0.0,                   0.0,       0.0,

                          sin(x)*cos(y)*cos(z),  x*x*2*y,   exp(x)*2*exp(2*y)*2*z,
                          exp(-2*x)*exp(y),      x,         x*2*z,
                          0.0,                   0.0,       0.0,

                          -sin(x)*sin(y)*sin(z), 0.0,      exp(x)*exp(2*y)*2,
                          0.0,                   0.0,      2.0*x*y,
                          0.0,                   0.0,      2.0}.data());

        for (const TensorFunc* fp : {static_cast<const TensorFunc*>(&f1),
                                     static_cast<const TensorFunc*>(&f2)}) {
          const utl::matrix4d<Real> hess = fp->hessian(X);
          for (size_t d1 = 1; d1 <= 3; ++d1)
            for (size_t d2 = 1; d2 <= 3; ++d2) {
              const Tensor dx = fp->dderiv(X,d1,d2);
              for (size_t i = 1; i <= 3; ++i)
                for (size_t j = 1; j <= 3; ++j) {
                  CHECK_THAT(dx(i,j), WithinRel(r(i,j,d1,d2)));
                  CHECK_THAT(hess(i,j,d1,d2), WithinRel(r(i,j,d1,d2)));
                }
            }
        }
      }
}


TEST_CASE("TestSTensorFunc.Evaluate2D")
{
  const char* func = "sin(x) | cos (y) | exp(z)";

  STensorFuncExpr f1(func);
  STensorFuncAd f2(func);

  CHECK(f1.isConstant());
  CHECK(f2.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const Tensor r({sin(x), exp(z), exp(z), cos(y)});

        for (const STensorFunc* fp : {static_cast<const STensorFunc*>(&f1),
                                      static_cast<const STensorFunc*>(&f2)}) {
          const SymmTensor fx = (*fp)(X);
          for (size_t i = 1; i <= 2; ++i)
            for (size_t j = 1; j <= 2; ++j)
              CHECK_THAT(fx(i,j), WithinRel(r(i,j)));
        }
      }
}


TEST_CASE("TestSTensorFunc.Evaluate3D")
{
  const char* func = "sin(x) | cos (y) | exp(z) |"
                     "sin(x)*sin(y) | sin(x)*cos(y) | exp(x)*exp(y)";

  STensorFuncExpr f1(func);
  STensorFuncAd f2(func);

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const Tensor r({sin(x), sin(x)*sin(y), exp(x)*exp(y),
                        sin(x)*sin(y), cos(y), sin(x)*cos(y),
                        exp(x)*exp(y), sin(x)*cos(y), exp(z)});

        for (const STensorFunc* fp : {static_cast<const STensorFunc*>(&f1),
                                      static_cast<const STensorFunc*>(&f2)}) {
          const SymmTensor fx = (*fp)(X);
          for (size_t i = 1; i <= 3; ++i)
            for (size_t j = 1; j <= 3; ++j)
              CHECK_THAT(fx(i,j), WithinRel(r(i,j)));
        }
      }
}


TEST_CASE("TestSTensorFunction.Gradient2D")
{
  const char* g   = "sin(x)*sin(y) | exp(x)*exp(2*y) | x*x*y*y";
  const char* g_x = "cos(x)*sin(y) | exp(x)*exp(2*y) | 2*x*y*y";
  const char* g_y = "sin(x)*cos(y) | 2*exp(x)*exp(2*y) | 2*x*x*y";

  STensorFuncExpr f1(g);
  f1.addDerivative(g_x,"",1);
  f1.addDerivative(g_y,"",2);

  STensorFuncAd f2(g);

  CHECK(f1.isConstant());
  CHECK(f2.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7}) {
      const Vec3 X(x,y);
      utl::matrix3d<Real> r(2,2,2);
      r.fill(std::array{cos(x)*sin(y), 2*x*y*y, 2*x*y*y, exp(x)*exp(2*y),
                        sin(x)*cos(y), 2*x*x*y, 2*x*x*y, 2*exp(x)*exp(2*y)}.data());

      for (const STensorFunc* fp : {static_cast<const STensorFunc*>(&f1),
                                    static_cast<const STensorFunc*>(&f2)}) {
        const utl::matrix3d<Real> grad = fp->gradient(X);
        for (size_t d = 1; d <= 2; ++d) {
          const SymmTensor dx = fp->deriv(X,d);
          for (size_t i = 1; i <= 2; ++i)
            for (size_t j = 1; j <= 2; ++j) {
              CHECK_THAT(dx(i,j), WithinRel(r(i,j,d)));
              CHECK_THAT(grad(i,j,d), WithinRel(r(i,j,d)));
            }
        }
      }
    }
}


TEST_CASE("TestSTensorFunction.Gradient3D")
{
  const char* g   = "sin(x)*sin(y)*sin(z) | exp(x)*exp(2*y)*exp(z) | x*x*y*y*z*z |"
                    "x*y*z | x*x*y*z | z";
  const char* g_x = "cos(x)*sin(y)*sin(z) | exp(x)*exp(2*y)*exp(z) | 2*x*y*y*z*z |"
                    "y*z | 2*x*y*z | 0";
  const char* g_y = "sin(x)*cos(y)*sin(z) | exp(x)*2*exp(2*y)*exp(z) | x*x*2*y*z*z |"
                    "x*z | x*x*z | 0";
  const char* g_z = "sin(x)*sin(y)*cos(z) | exp(x)*exp(2*y)*exp(z) | x*x*y*y*2*z |"
                    "x*y | x*x*y | 1";

  STensorFuncExpr f1(g);
  f1.addDerivative(g_x,"",1);
  f1.addDerivative(g_y,"",2);
  f1.addDerivative(g_z,"",3);

  STensorFuncAd f2(g);

  CHECK(f1.isConstant());
  CHECK(f2.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        utl::matrix3d<Real> r(3,3,3);
        r.fill(std::array{cos(x)*sin(y)*sin(z), y*z,                      0.0,
                          y*z,                  exp(x)*exp(2*y)*exp(z),   2*x*y*z,
                          0.0,                  2*x*y*z,                  2*x*y*y*z*z,

                          sin(x)*cos(y)*sin(z), x*z,                      0.0,
                          x*z,                  exp(x)*2*exp(2*y)*exp(z), x*x*z,
                          0.0,                  x*x*z,                    x*x*2*y*z*z,

                          sin(x)*sin(y)*cos(z), x*y,                      1.0,
                          x*y,                  exp(x)*exp(2*y)*exp(z),   x*x*y,
                          1.0,                  x*x*y,                    x*x*y*y*2*z}.data());

        for (const STensorFunc* fp : {static_cast<const STensorFunc*>(&f1),
                                      static_cast<const STensorFunc*>(&f2)}) {
          const utl::matrix3d<Real> grad = fp->gradient(X);
          for (size_t d = 1; d <= 3; ++d) {
            const SymmTensor dx = fp->deriv(X,d);
            for (size_t i = 1; i <= 3; ++i)
              for (size_t j = 1; j <= 3; ++j) {
                CHECK_THAT(dx(i,j), WithinRel(r(i,j,d)));
                CHECK_THAT(grad(i,j,d), WithinRel(r(i,j,d)));
              }
          }
        }
      }
}


TEST_CASE("TestSTensorFunction.Hessian2D")
{
  const char* g    = " sin(x)*sin(y) | exp(x)*exp(2*y)   | x*x*y*y";
  const char* g_xx = "-sin(x)*sin(y) | exp(x)*exp(2*y)   | 2*y*y";
  const char* g_xy = " cos(x)*cos(y) | exp(x)*2*exp(2*y) | 2*x*2*y";
  const char* g_yy = "-sin(x)*sin(y) | exp(x)*4*exp(2*y) | x*x*2";

  STensorFuncExpr f1(g);
  f1.addDerivative(g_xx,"",1,1);
  f1.addDerivative(g_xy,"",1,2);
  f1.addDerivative(g_yy,"",2,2);

  STensorFuncAd f2(g);

  CHECK(f1.isConstant());
  CHECK(f2.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7}) {
      const Vec3 X(x,y);
      utl::matrix4d<Real> r(2,2,2,2);
      r.fill(std::array{-sin(x)*sin(y), 2*y*y,
                         2*y*y,         exp(x)*exp(2*y),

                         cos(x)*cos(y), 2*x*2*y,
                         2*x*2*y,       exp(x)*2*exp(2*y),

                         cos(x)*cos(y), 2*x*2*y,
                         2*x*2*y,       exp(x)*2*exp(2*y),

                        -sin(x)*sin(y), x*x*2,
                         x*x*2,         exp(x)*4*exp(2*y)}.data());

      for (const STensorFunc* fp : {static_cast<const STensorFunc*>(&f1),
                                    static_cast<const STensorFunc*>(&f2)}) {
        const utl::matrix4d<Real> hess = fp->hessian(X);
        for (size_t d1 = 1; d1 <= 2; ++d1)
          for (size_t d2 = 1; d2 <= 2; ++d2) {
            const SymmTensor dx = fp->dderiv(X,d1,d2);
            for (size_t i = 1; i <= 2; ++i)
              for (size_t j = 1; j <= 2; ++j) {
                CHECK_THAT(dx(i,j), WithinRel(r(i,j,d1,d2)));
                CHECK_THAT(hess(i,j,d1,d2), WithinRel(r(i,j,d1,d2)));
              }
          }
      }
    }
}


TEST_CASE("TestSTensorFunction.Hessian3D")
{
  const char* g    = "sin(x)*sin(y)*sin(z) | exp(x)*exp(2*y)*exp(z) | x*x*y*y*z*z |"
                     "x*y*z | x*x*y*z | z";
  const char* g_xx = "-sin(x)*sin(y)*sin(z) | exp(x)*exp(2*y)*exp(z) | 2*y*y*z*z |"
                     "0.0 | 2*y*z | 0.0";
  const char* g_xy = "cos(x)*cos(y)*sin(z) | exp(x)*2*exp(2*y)*exp(z) | 2*x*2*y*z*z |"
                     "z | 2*x*z | 0.0";
  const char* g_xz = "cos(x)*sin(y)*cos(z) | exp(x)*exp(2*y)*exp(z) | 2*x*y*y*2*z |"
                     "y | 2*x*y | 0.0";
  const char* g_yy = "-sin(x)*sin(y)*sin(z) | exp(x)*4*exp(2*y)*exp(z) | x*x*2*z*z |"
                     "0.0 | 0.0 | 0.0";
  const char* g_yz = "sin(x)*cos(y)*cos(z) | exp(x)*2*exp(2*y)*exp(z) | x*x*2*y*2*z |"
                     "x | x*x | 0.0";
  const char* g_zz = "-sin(x)*sin(y)*sin(z) | exp(x)*exp(2*y)*exp(z) | x*x*y*y*2 |"
                     "0.0 | 0.0 | 0.0";

  STensorFuncExpr f1(g);
  f1.addDerivative(g_xx,"",1,1);
  f1.addDerivative(g_xy,"",1,2);
  f1.addDerivative(g_xz,"",1,3);
  f1.addDerivative(g_yy,"",2,2);
  f1.addDerivative(g_yz,"",2,3);
  f1.addDerivative(g_zz,"",3,3);

  STensorFuncAd f2(g);

  CHECK(f1.isConstant());
  CHECK(f2.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        utl::matrix4d<Real> r(3,3,3,3);
        r.fill(std::array{-sin(x)*sin(y)*sin(z), 0.0, 0.0,
                           0.0,                  exp(x)*exp(2*y)*exp(z),   2*y*z,
                           0.0,                  2*y*z,                    2*y*y*z*z,

                           cos(x)*cos(y)*sin(z), z,                        0.0,
                           z,                    exp(x)*2*exp(2*y)*exp(z), 2*x*z,
                           0.0,                  2*x*z,                    2*x*2*y*z*z,

                           cos(x)*sin(y)*cos(z), y,                        0.0,
                           y,                    exp(x)*exp(2*y)*exp(z),   2*x*y,
                           0.0,                  2*x*y,                    2*x*y*y*2*z,

                           cos(x)*cos(y)*sin(z), z,                        0.0,
                           z,                    exp(x)*2*exp(2*y)*exp(z), 2*x*z,
                           0.0,                  2*x*z,                    2*x*2*y*z*z,

                          -sin(x)*sin(y)*sin(z), 0.0,                      0.0,
                           0.0,                  exp(x)*4*exp(2*y)*exp(z), 0.0,
                           0.0,                  0.0,                      x*x*2*z*z,

                           sin(x)*cos(y)*cos(z), x,                        0.0,
                           x,                    exp(x)*2*exp(2*y)*exp(z), x*x,
                           0.0,                  x*x,                      x*x*2*y*2*z,

                           cos(x)*sin(y)*cos(z), y,                        0.0,
                           y,                    exp(x)*exp(2*y)*exp(z),   2*x*y,
                           0.0,                  2*x*y,                    2*x*y*y*2*z,

                           sin(x)*cos(y)*cos(z), x,                        0.0,
                           x,                    exp(x)*2*exp(2*y)*exp(z), x*x,
                           0.0,                  x*x,                      x*x*2*y*2*z,

                          -sin(x)*sin(y)*sin(z), 0.0,                      0.0,
                           0.0,                  exp(x)*exp(2*y)*exp(z),   0.0,
                           0.0,                  0.0,                      x*x*y*y*2}.data());

        for (const STensorFunc* fp : {static_cast<const STensorFunc*>(&f1),
                                      static_cast<const STensorFunc*>(&f2)}) {
          const utl::matrix4d<Real> hess = fp->hessian(X);
          for (size_t d1 = 1; d1 <= 3; ++d1)
            for (size_t d2 = 1; d2 <= 3; ++d2) {
              const SymmTensor dx = fp->dderiv(X,d1,d2);
              for (size_t i = 1; i <= 3; ++i)
                for (size_t j = 1; j <= 3; ++j) {
                  CHECK_THAT(dx(i,j), WithinRel(r(i,j,d1,d2)));
                  CHECK_THAT(hess(i,j,d1,d2), WithinRel(r(i,j,d1,d2)));
                }
            }
        }
      }
}
