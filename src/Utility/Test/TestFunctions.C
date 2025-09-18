//==============================================================================
//!
//! \file TestFunctions.C
//!
//! \date Apr 28 2018
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Tests for parsing of functions.
//!
//==============================================================================

#include "ExprFunctions.h"
#include "TensorFunction.h"
#include "Functions.h"

#include <autodiff/reverse/var.hpp>

#include <cstdlib>
#include <cmath>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;


using EvalFuncAd = EvalFuncScalar<autodiff::var>;
using EvalFunctionAd = EvalFuncSpatial<autodiff::var>;
using VecFuncExprAd = EvalMultiFunction<VecFunc,Vec3,autodiff::var>;
using TensorFuncExprAd = EvalMultiFunction<TensorFunc,Tensor,autodiff::var>;
using STensorFuncExprAd = EvalMultiFunction<STensorFunc,SymmTensor,autodiff::var>;


TEST_CASE("TestScalarFunc.ParseDerivative")
{
  const char* func1 = "sin(1.5*t)*t";
  const char* func2 = "sin(1.5*t)*t:1.5*cos(1.5*t)*t+sin(1.5*t)";

  ScalarFunc* f1 = utl::parseTimeFunc(func1,"expression");
  ScalarFunc* f2 = utl::parseTimeFunc(func2,"expression");

  EvalFuncAd f3(func1,"t");

  REQUIRE(f1 != nullptr);
  REQUIRE(f2 != nullptr);

  REQUIRE(!f1->isConstant());
  REQUIRE(!f2->isConstant());

  double t = 0.0;
  for (int i = 0; i < 20; i++)
  {
    t += 0.314*(double)random()/(double)RAND_MAX;
    std::cout <<"f("<< t <<") = "<< (*f1)(t)
              <<"  f'("<< t <<") = "<< f1->deriv(t) << std::endl;
    REQUIRE_THAT((*f1)(t), WithinRel(sin(1.5*t)*t));
    REQUIRE_THAT((*f2)(t), WithinRel(sin(1.5*t)*t));
    REQUIRE_THAT(f1->deriv(t), WithinRel(1.5*cos(1.5*t)*t+sin(1.5*t), 1e-6));
    REQUIRE_THAT(f2->deriv(t), WithinRel(1.5*cos(1.5*t)*t+sin(1.5*t)));
    REQUIRE_THAT(f3.deriv(t), WithinRel(1.5*cos(1.5*t)*t+sin(1.5*t)));
  }
}


TEST_CASE("TestScalarFunc.ParseFunction")
{
  std::cout <<"Parsing scalar function: ";
  ScalarFunc* f1 = utl::parseTimeFunc("1.2 100.0","Dirac");
  std::cout <<"Parsing scalar function: ";
  ScalarFunc* f2 = utl::parseTimeFunc("5.0 100.0","Ramp");

  REQUIRE_THAT((*f1)(1.1), WithinAbs(0.0, 1e-12));
  REQUIRE_THAT((*f1)(1.2), WithinRel(100.0));
  REQUIRE_THAT((*f1)(1.3), WithinAbs(0.0, 1e-12));
  REQUIRE_THAT((*f2)(2.5), WithinRel(50.0));
  REQUIRE_THAT((*f2)(5.0), WithinRel(100.0));
  REQUIRE_THAT((*f2)(7.0), WithinRel(100.0));
}


TEST_CASE("TestRealFunc.ParseFunction")
{
  std::cout <<"Parsing real function";
  RealFunc* f1 = utl::parseRealFunc("100.0 1.2","Dirac");
  std::cout <<"\nParsing real function";
  RealFunc* f2 = utl::parseRealFunc("100.0 5.0","Ramp",false);
  std::cout << std::endl;

  REQUIRE_THAT((*f1)(Vec4(1.1)), WithinAbs(0.0, 1e-12));
  REQUIRE_THAT((*f1)(Vec4(1.2)), WithinRel(100.0));
  REQUIRE_THAT((*f1)(Vec4(1.3)), WithinAbs(0.0, 1e-12));
  REQUIRE_THAT((*f2)(Vec4(2.5)), WithinRel(50.0));
  REQUIRE_THAT((*f2)(Vec4(5.0)), WithinRel(100.0));
  REQUIRE_THAT((*f2)(Vec4(7.0)), WithinRel(100.0));
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

  EvalFunctionAd f2(g);

  REQUIRE(f1.isConstant());
  REQUIRE(f2.isConstant());

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
            REQUIRE_THAT(fp->deriv(X, i), WithinRel(r[i-1]));
            REQUIRE_THAT(grad[i-1], WithinRel(r[i-1]));
          }
        }
      }
}


TEST_CASE("TestRealFunc.GradientFD")
{
  const char* f1 = "sin(x)*sin(y)*sin(z)";

  const double eps = 1e-6;

  EvalFunction f(f1, eps);

  REQUIRE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const Vec3 Xp(x+0.5*eps, y+0.5*eps, z+0.5*eps);
        const Vec3 Xm(x-0.5*eps, y-0.5*eps, z-0.5*eps);
        const Vec3 r((sin(Xp.x) - sin(Xm.x))*sin(y)*sin(z) / eps,
                     sin(x)*(sin(Xp.y) - sin(Xm.y))*sin(z) / eps,
                     sin(x)*sin(y)*(sin(Xp.z) - sin(Xm.z)) / eps);

        const Vec3 grad = f.gradient(X);
        for (size_t i = 1; i <= 3; ++i) {
          REQUIRE_THAT(f.deriv(X, i), WithinRel(r[i-1], 1e-8));
          REQUIRE_THAT(grad[i-1], WithinRel(r[i-1], 1e-8));
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

  EvalFunctionAd f2(g);

  REQUIRE(f1.isConstant());
  REQUIRE(f2.isConstant());

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
              REQUIRE_THAT(fp->dderiv(X,i,j), WithinRel(r(i,j)));
              REQUIRE_THAT(hess(i,j), WithinRel(r(i,j)));
            }
          }
      }
}


TEST_CASE("TestVecFunc.Evaluate")
{
  const char* func = "sin(x) | cos (y) | exp(z)";

  VecFuncExpr f1(func);
  VecFuncExprAd f2(func);

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const Vec3 r(sin(x), cos(y), exp(z));
        for (const VecFunc* fp : {static_cast<const VecFunc*>(&f1),
                                  static_cast<const VecFunc*>(&f2)}) {
          const Vec3 fx = (*fp)(X);
          REQUIRE_THAT(fx.x, WithinRel(r.x));
          REQUIRE_THAT(fx.y, WithinRel(r.y));
          REQUIRE_THAT(fx.z, WithinRel(r.z));
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

  VecFuncExprAd f2(g);

  REQUIRE(f1.isConstant());
  REQUIRE(f2.isConstant());

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
            REQUIRE_THAT(dx[i-1], WithinRel(r(i,d)));
            REQUIRE_THAT(grad(i,d), WithinRel(r(i,d)));
          }
        }
      }
    }
}


TEST_CASE("TestVecFunction.Gradient2DFD")
{
  const char* g = "sin(x)*sin(y) | x*x*y*y";

  const double eps = 1e-6;

  VecFuncExpr f(g,"",eps);

  REQUIRE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7}) {
      const Vec3 X(x,y);
      const Vec3 Xp(x + 0.5*eps, y + 0.5*eps);
      const Vec3 Xm(x - 0.5*eps, y - 0.5*eps);
      const Tensor r({(sin(Xp.x) - sin(Xm.x))*sin(y) / eps,
                      (Xp.x*Xp.x - Xm.x*Xm.x)*y*y / eps,
                      sin(x)*(sin(Xp.y) - sin(Xm.y)) / eps,
                      x*x*(Xp.y*Xp.y - Xm.y*Xm.y) / eps});

      const Tensor grad = f.gradient(X);
      for (size_t d = 1; d <= 2; ++d) {
        const Vec3 dx = f.deriv(X,d);
        for (size_t i = 1; i <= 2; ++i) {
          REQUIRE_THAT(dx[i-1], WithinRel(r(i,d), 1e-8));
          REQUIRE_THAT(grad(i,d), WithinRel(r(i,d), 1e-8));
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

  VecFuncExprAd f2(g);

  REQUIRE(f1.isConstant());
  REQUIRE(f2.isConstant());

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
              REQUIRE_THAT(dx[i-1], WithinRel(r(i,d)));
              REQUIRE_THAT(grad(i,d), WithinRel(r(i,d)));
            }
          }
        }
      }
}


TEST_CASE("TestVecFunction.Gradient3DFD")
{
  const char* g   = "sin(x)*sin(y)*sin(z) | x*x*y*y*z*z*z | exp(-x)*exp(2*y)*exp(z*z)";

  const double eps = 1e-6;
  VecFuncExpr f(g,"",eps);

  REQUIRE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const Vec3 Xp(x + 0.5*eps, y + 0.5*eps, z + 0.5*eps);
        const Vec3 Xm(x - 0.5*eps, y - 0.5*eps, z - 0.5*eps);
        const Tensor r({(sin(Xp.x) - sin(Xm.x))*sin(y)*sin(z) / eps,
                        (Xp.x*Xp.x - Xm.x*Xm.x)*X.y*X.y*X.z*X.z*X.z / eps,
                        (exp(-Xp.x) - exp(-Xm.x))*exp(2*y)*exp(z*z) / eps,
                        sin(x)*(sin(Xp.y) - sin(Xm.y))*sin(z) / eps,
                        x*x*(Xp.y*Xp.y - Xm.y*Xm.y)*z*z*z / eps,
                        exp(-x)*(exp(2*Xp.y) - exp(2*Xm.y))*exp(z*z) / eps,
                        sin(x)*sin(y)*(sin(Xp.z) - sin(Xm.z)) / eps,
                        x*x*y*y*(Xp.z*Xp.z*Xp.z - Xm.z*Xm.z*Xm.z) / eps,
                        exp(-x)*exp(2*y)*(exp(Xp.z*Xp.z) - exp(Xm.z*Xm.z)) / eps});

        const Tensor grad = f.gradient(X);
        for (size_t d = 1; d <= 3; ++d) {
          const Vec3 dx = f.deriv(X,d);
          for (size_t i = 1; i <= 3; ++i) {
            REQUIRE_THAT(dx[i-1], WithinRel(r(i,d), 1e-8));
            REQUIRE_THAT(grad(i,d), WithinRel(r(i,d), 1e-8));
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

  VecFuncExprAd f2(g);

  REQUIRE(f1.isConstant());
  REQUIRE(f2.isConstant());

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
              REQUIRE_THAT(hess(i,d1,d2), WithinRel(r(i,d1,d2)));
              REQUIRE_THAT(d2x[i-1], WithinRel(r(i,d1,d2)));
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

  VecFuncExprAd f2(g);

  REQUIRE(f1.isConstant());
  REQUIRE(f2.isConstant());

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
                REQUIRE_THAT(d2x[i-1], WithinRel(r(i,d1,d2)));
                REQUIRE_THAT(hess(i,d1,d2), WithinRel(r(i,d1,d2)));
              }
            }
        }
      }
 }


TEST_CASE("TestVecFuncExpr.NumDimensions")
{
  const char* func1 = "x";
  const char* func2 = "x | y";
  const char* func3 = "x | y | z";

  VecFuncExpr f1(func1);
  REQUIRE(f1.getNoSpaceDim() == 1);
  REQUIRE(f1.dim() ==  1);
  VecFuncExpr f2(func2);
  REQUIRE(f2.getNoSpaceDim() == 2);
  REQUIRE(f2.dim() == 2);
  VecFuncExpr f3(func3);
  REQUIRE(f3.getNoSpaceDim() == 3);
  REQUIRE(f3.dim() == 3);
}


TEST_CASE("TestVecFuncExpr.TimeDerivative")
{
  const char* g   = "sin(x)*sin(y)*sin(z)*sin(t) | x*x*y*y*z*t*t | exp(-2*t)";
  const char* g_t = "sin(x)*sin(y)*sin(z)*cos(t) | x*x*y*y*z*2*t | -2*exp(-2*t)";

  VecFuncExpr f(g);
  f.addDerivative(g_t,"",4);

  REQUIRE(!f.isConstant());

  for (double t : {0.1, 0.2, 0.3})
    for (double x : {0.1, 0.2, 0.3})
      for (double y : {0.5, 0.6, 0.7})
        for (double z : {0.8, 0.9, 1.0}) {
          const Vec4 X(x,y,z,t);
          const Vec3 r(sin(x)*sin(y)*sin(z)*cos(t),
                       x*x*y*y*z*2*t,
                       -2*exp(-2*t));
          const Vec3 dt = f.timeDerivative(X);
          REQUIRE_THAT(dt[0], WithinRel(r[0]));
          REQUIRE_THAT(dt[1], WithinRel(r[1]));
          REQUIRE_THAT(dt[2], WithinRel(r[2]));
        }
}


TEST_CASE("TestVecFuncExpr.TimeDerivativeFD")
{
  const char* g   = "sin(x)*sin(y)*sin(z)*sin(t) | x*x*y*y*z*t*t | exp(-2*t)";

  const double eps = 1e-6;
  VecFuncExpr f(g,"",1e-8,eps);

  REQUIRE(!f.isConstant());

  for (double t : {0.1, 0.2, 0.3}) {
    const double tp = t + 0.5*eps;
    const double tm = t - 0.5*eps;
    for (double x : {0.1, 0.2, 0.3})
      for (double y : {0.5, 0.6, 0.7})
        for (double z : {0.8, 0.9, 1.0}) {
          const Vec4 X(x,y,z,t);
          Vec3 r(sin(x)*sin(y)*sin(z)*(sin(tp) - sin(tm)),
                 x*x*y*y*z*(tp*tp - tm*tm),
                 exp(-2*tp) - exp(-2*tm));
          r *= 1.0 / eps;

          const Vec3 dt = f.timeDerivative(X);
          REQUIRE_THAT(dt[0], WithinRel(r[0], 1e-8));
          REQUIRE_THAT(dt[1], WithinRel(r[1], 1e-8));
          REQUIRE_THAT(dt[2], WithinRel(r[2], 1e-8));
        }
  }
}


TEST_CASE("TestTensorFunc.Evaluate")
{
  const char* func = "sin(x) | cos (y) | exp(z) | sin(x)*cos(y)";

  TensorFuncExpr f1(func);
  TensorFuncExprAd f2(func);

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
              REQUIRE_THAT(fx(i,j), WithinRel(r(i,j)));
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

  TensorFuncExprAd f2(g);

  REQUIRE(f1.isConstant());

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
              REQUIRE_THAT(dx(i,j), WithinRel(r(i,j,d)));
              REQUIRE_THAT(grad(i,j,d), WithinRel(r(i,j,d)));
            }
        }
      }
    }
}


TEST_CASE("TestTensorFunction.Gradient2DFD")
{
  const char* g   = "sin(x)*sin(y) | x*x*y*y | exp(x)*exp(2*y) | exp(-2*x)*exp(y)";

  const double eps = 1e-6;
  TensorFuncExpr f(g,"",eps);

  REQUIRE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7}) {
      const Vec3 X(x,y);
      const Vec3 Xp(x + 0.5*eps, y + 0.5*eps);
      const Vec3 Xm(x - 0.5*eps, y - 0.5*eps);
      utl::matrix3d<Real> r(2,2,2);
      r.fill(std::array{(sin(Xp.x) - sin(Xm.x))*sin(y),
                        (Xp.x*Xp.x - Xm.x*Xm.x)*y*y,
                        (exp(Xp.x) - exp(Xm.x))*exp(2*y),
                        (exp(-2*Xp.x) - exp(-2*Xm.x))*exp(y),
                        sin(x)*(sin(Xp.y) - sin(Xm.y)),
                        x*x*(Xp.y*Xp.y - Xm.y*Xm.y),
                        exp(x)*(exp(2*Xp.y) - exp(2*Xm.y)),
                        exp(-2*x)*(exp(Xp.y) - exp(Xm.y))}.data());
      r.multiply(1.0 / eps);

      const utl::matrix3d<Real> grad = f.gradient(X);
      for (size_t d = 1; d <= 2; ++d) {
        const Tensor dx = f.deriv(X,d);
        for (size_t i = 1; i <= 2; ++i)
          for (size_t j = 1; j <= 2; ++j) {
            REQUIRE_THAT(dx(i,j), WithinRel(r(i,j,d), 1e-8));
            REQUIRE_THAT(grad(i,j,d), WithinRel(r(i,j,d), 1e-8));
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

  TensorFuncExprAd f2(g);

  REQUIRE(f1.isConstant());
  REQUIRE(f2.isConstant());

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
                REQUIRE_THAT(dx(i,j), WithinRel(r(i,j,d)));
                REQUIRE_THAT(grad(i,j,d), WithinRel(r(i,j,d)));
              }
          }
        }
      }
}


TEST_CASE("TestTensorFunction.Gradient3DFD")
{
  const char* g   = "sin(x)*sin(y)*sin(z)  | x*x*y*y*z | exp(x)*exp(2*y)*z*z |"
                    "exp(-2*x)*exp(y)*z    | x*y*z     | x*y*z*z |"
                    "x                     | y         | z";

  const double eps = 1e-6;
  TensorFuncExpr f(g,"",eps);

  REQUIRE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const Vec3 Xp(x + 0.5*eps, y + 0.5*eps, z + 0.5*eps);
        const Vec3 Xm(x - 0.5*eps, y - 0.5*eps, z - 0.5*eps);
        utl::matrix3d<Real> r(3,3,3);
        r.fill(std::array{(sin(Xp.x) - sin(Xm.x))*sin(y)*sin(z),
                          (Xp.x*Xp.x - Xm.x*Xm.x)*y*y*z,
                          (exp(Xp.x) - exp(Xm.x))*exp(2*y)*z*z,
                          (exp(-2*Xp.x) - exp(-2*Xm.x))*exp(y)*z,
                          (Xp.x - Xm.x)*y*z,
                          (Xp.x - Xm.x)*y*z*z,
                          Xp.x -Xm.x,
                          0.0,
                          0.0,
                          sin(x)*(sin(Xp.y) - sin(Xm.y))*sin(z),
                          x*x*(Xp.y*Xp.y - Xm.y*Xm.y)*z,
                          exp(x)*(exp(2*Xp.y) - exp(2*Xm.y))*z*z,
                          exp(-2*x)*(exp(Xp.y) - exp(Xm.y))*z,
                          x*(Xp.y - Xm.y)*z,
                          x*(Xp.y - Xm.y)*z*z,
                          0.0,
                          Xp.y - Xm.y,
                          0.0,
                          sin(x)*sin(y)*(sin(Xp.z) - sin(Xm.z)),
                          x*x*y*y*(Xp.z - Xm.z),
                          exp(x)*exp(2*y)*(Xp.z*Xp.z - Xm.z*Xm.z),
                          exp(-2*x)*exp(y)*(Xp.z - Xm.z),
                          x*y*(Xp.z - Xm.z),
                          x*y*(Xp.z*Xp.z - Xm.z*Xm.z),
                          0.0,
                          0.0,
                          Xp.z - Xm.z}.data());
        r.multiply(1.0 / eps);

        const utl::matrix3d<Real> grad = f.gradient(X);
        for (size_t d = 1; d <= 3; ++d) {
          const Tensor dx = f.deriv(X,d);
          for (size_t i = 1; i <= 3; ++i)
            for (size_t j = 1; j <= 3; ++j) {
              REQUIRE_THAT(dx(i,j), WithinRel(r(i,j,d), 1e-8));
              REQUIRE_THAT(grad(i,j,d), WithinRel(r(i,j,d), 1e-8));
            }
        }
    }
}


TEST_CASE("TestTensorFunction.Hessian2D")
{
  const char* g    = "sin(x)*sin(y)   | x*x*y*y | exp(x)*exp(2*y)   | exp(-2*x)*exp(y)";
  const char* g_xx = "-sin(x)*sin(y) | 2*y*y   | exp(x)*exp(2*y)   | 4*exp(-2*x)*exp(y)";
  const char* g_xy = "cos(x)*cos(y)  | 2*x*2*y | exp(x)*2*exp(2*y) | -2*exp(-2*x)*exp(y)";
  const char* g_yy = "-sin(x)*sin(y) | 2*x*x   | exp(x)*4*exp(2*y) | exp(-2*x)*exp(y)";

  TensorFuncExpr f1(g);
  f1.addDerivative(g_xx,"",1,1);
  f1.addDerivative(g_yy,"",2,2);
  f1.addDerivative(g_xy,"",1,2);

  TensorFuncExprAd f2(g);

  REQUIRE(f1.isConstant());
  REQUIRE(f2.isConstant());

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
                REQUIRE_THAT(dx(i,j), WithinRel(r(i,j,d1,d2)));
                REQUIRE_THAT(hess(i,j,d1,d2), WithinRel(r(i,j,d1,d2)));
              }
          }
      }
    }
}


TEST_CASE("TestTensorFunction.Hessian3D")
{
  const char* g    = "sin(x)*sin(y)*sin(z)  | x*x*y*y*z | exp(x)*exp(2*y)*z*z |"
                     "exp(-2*x)*exp(y)*z     | x*y*z     | x*y*z*z |"
                     "x*x                    | y*y       | z*z";
  const char* g_xx = "-sin(x)*sin(y)*sin(z) | 2*y*y*z   | exp(x)*exp(2*y)*z*z |"
                     "4*exp(-2*x)*exp(y)*z   | 0.0       | 0.0 |"
                     "2.0                    | 0.0       | 0.0";
  const char* g_xy = "cos(x)*cos(y)*sin(z)  | 2*x*2*y*z | exp(x)*2*exp(2*y)*z*z |"
                     "-2*exp(-2*x)*exp(y)*z  | z         | z*z |"
                     "0.0                    | 0.0       | 0.0";
  const char* g_xz = "cos(x)*sin(y)*cos(z)  | 2*x*y*y   | exp(x)*exp(2*y)*2*z |"
                     "-2*exp(-2*x)*exp(y)    | y         | y*2*z |"
                     "0.0                    | 0.0       | 0.0";
  const char* g_yy = "-sin(x)*sin(y)*sin(z) | x*x*2*z   | exp(x)*4*exp(2*y)*z*z |"
                     "exp(-2*x)*exp(y)*z     | 0.0       | 0.0 |"
                     "0.0                    | 2.0       | 0.0";
  const char* g_yz = "sin(x)*cos(y)*cos(z)  | x*x*2*y   | exp(x)*2*exp(2*y)*2*z |"
                     "exp(-2*x)*exp(y)       | x         | x*2*z |"
                     "0.0                    | 0.0       | 0.0";
  const char* g_zz = "-sin(x)*sin(y)*sin(z) | 0.0       | exp(x)*exp(2*y)*2 |"
                     "0.0                    | 0.0       | 2.0*x*y |"
                     "0.0                    | 0.0       | 2.0";

  TensorFuncExpr f1(g);
  f1.addDerivative(g_xx,"",1,1);
  f1.addDerivative(g_yy,"",2,2);
  f1.addDerivative(g_zz,"",3,3);
  f1.addDerivative(g_xy,"",1,2);
  f1.addDerivative(g_xz,"",1,3);
  f1.addDerivative(g_yz,"",2,3);

  TensorFuncExprAd f2(g);

  REQUIRE(f1.isConstant());
  REQUIRE(f2.isConstant());

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
                  REQUIRE_THAT(dx(i,j), WithinRel(r(i,j,d1,d2)));
                  REQUIRE_THAT(hess(i,j,d1,d2), WithinRel(r(i,j,d1,d2)));
                }
            }
        }
      }
}


TEST_CASE("TestTensorFunction.TimeDerivative")
{
  const char* g   = "sin(x)*sin(y)*sin(z)*sin(t)  | x*x*y*y*z*t*t | exp(x)*exp(2*y)*z*z*t |"
                    "exp(-2*x)*exp(y)*z*sin(t)    | x*y*z*t*t     | x*y*z*z*t |"
                    "x*sin(t)                     | y*t*t         | z*t";
  const char* g_t = "sin(x)*sin(y)*sin(z)*cos(t)  | x*x*y*y*z*2*t | exp(x)*exp(2*y)*z*z |"
                    "exp(-2*x)*exp(y)*z*cos(t)    | x*y*z*2*t     | x*y*z*z |"
                    "x*cos(t)                     | y*2*t         | z";

  TensorFuncExpr f(g);
  f.addDerivative(g_t,"",4);

  REQUIRE(!f.isConstant());

  for (double t : {0.1, 0.2, 0.3})
    for (double x : {0.1, 0.2, 0.3})
      for (double y : {0.5, 0.6, 0.7})
        for (double z : {0.8, 0.9, 1.0}) {
          const Vec4 X(x,y,z,t);
          const Tensor r({sin(x)*sin(y)*sin(z)*cos(t), x*x*y*y*z*2*t, exp(x)*exp(2*y)*z*z,
                          exp(-2*x)*exp(y)*z*cos(t),   x*y*z*2*t,     x*y*z*z,
                          x*cos(t),                    y*2*t,         z});

          const Tensor dt = f.timeDerivative(X);
          for (size_t i = 1; i <= 3; ++i)
            for (size_t j = 1; j <= 3; ++j)
              REQUIRE_THAT(dt(i,j), WithinRel(r(i,j)));
      }
}


TEST_CASE("TestTensorFuncExpr.NumDimensions")
{
  const char* func1 = "x";
  const char* func2 = "x | y | z | x";
  const char* func3 = "x | y | z | x | y | z | x | y | z";

  TensorFuncExpr f1(func1);
  REQUIRE(f1.getNoSpaceDim() == 1);
  REQUIRE(f1.dim() == 1);
  TensorFuncExpr f2(func2);
  REQUIRE(f2.getNoSpaceDim() == 2);
  REQUIRE(f2.dim() == 4);
  TensorFuncExpr f3(func3);
  REQUIRE(f3.getNoSpaceDim() == 3);
  REQUIRE(f3.dim() == 9);
}


TEST_CASE("TestSTensorFunc.Evaluate2D")
{
  const char* func = "sin(x) | cos (y) | exp(z)";

  STensorFuncExpr f1(func);
  STensorFuncExprAd f2(func);

  REQUIRE(f1.isConstant());
  REQUIRE(f2.isConstant());

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
              REQUIRE_THAT(fx(i,j), WithinRel(r(i,j)));
        }
      }
}


TEST_CASE("TestSTensorFunc.Evaluate2Dzz")
{
  const char* func = "sin(x) | cos (y) | sin(x)*cos(y) | exp(z)";

  STensorFuncExpr f(func);

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const SymmTensor fx = f(X);
        const Tensor r({sin(x), exp(z), exp(z), cos(y)});
        for (size_t i = 1; i <= 2; ++i)
          for (size_t j = 1; j <= 2; ++j)
            REQUIRE_THAT(fx(i,j), WithinRel(r(i,j)));
        REQUIRE_THAT(fx(3,3), WithinRel(sin(x)*cos(y)));
      }
}


TEST_CASE("TestSTensorFunc.Evaluate3D")
{
  const char* func = "sin(x) | cos (y) | exp(z) |"
                     "sin(x)*sin(y) | sin(x)*cos(y) | exp(x)*exp(y)";

  STensorFuncExpr f1(func);
  STensorFuncExprAd f2(func);

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
              REQUIRE_THAT(fx(i,j), WithinRel(r(i,j)));
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

  STensorFuncExprAd f2(g);

  REQUIRE(f1.isConstant());
  REQUIRE(f2.isConstant());

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
              REQUIRE_THAT(dx(i,j), WithinRel(r(i,j,d)));
              REQUIRE_THAT(grad(i,j,d), WithinRel(r(i,j,d)));
            }
        }
      }
    }
}


TEST_CASE("TestSTensorFunction.Gradient2Dzz")
{
  const char* g   = "sin(x)*sin(y) | exp(x)*exp(2*y) | sin(x)*sin(y) | x*x*y*y";
  const char* g_x = "cos(x)*sin(y) | exp(x)*exp(2*y) | cos(x)*sin(y) | 2*x*y*y";
  const char* g_y = "sin(x)*cos(y) | 2*exp(x)*exp(2*y) | sin(x)*cos(y) | 2*x*x*y";

  STensorFuncExpr f(g);
  f.addDerivative(g_x,"",1);
  f.addDerivative(g_y,"",2);

  REQUIRE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7}) {
      const Vec3 X(x,y);
      utl::matrix3d<Real> r(2,2,2);
      r.fill(std::array{cos(x)*sin(y), 2*x*y*y, 2*x*y*y, exp(x)*exp(2*y),
                        sin(x)*cos(y), 2*x*x*y, 2*x*x*y, 2*exp(x)*exp(2*y)}.data());

      for (size_t d = 1; d <= 2; ++d) {
        const SymmTensor dx = f.deriv(X,d);
        for (size_t i = 1; i <= 2; ++i)
          for (size_t j = 1; j <= 2; ++j)
            REQUIRE_THAT(dx(i,j), WithinRel(r(i,j,d), 1e-12));
        REQUIRE_THAT(dx(3,3), WithinRel((d == 1 ? cos(x) : sin(x)) * (d == 2 ? cos(y) : sin(y)), 1e-12));
      }
    }
}


TEST_CASE("TestSTensorFunction.Gradient2DFD")
{
  const char* g   = "sin(x)*sin(y) | exp(x)*exp(2*y) | x*x*y*y";

  const double eps = 1e-6;
  STensorFuncExpr f(g,"",eps);

  REQUIRE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7}) {
      const Vec3 X(x,y);
      const Vec3 Xp(x + 0.5*eps, y + 0.5*eps);
      const Vec3 Xm(x - 0.5*eps, y - 0.5*eps);
      utl::matrix3d<Real> r(2,2,2);
      r.fill(std::array{(sin(Xp.x) - sin(Xm.x))*sin(y),
                        (Xp.x*Xp.x - Xm.x*Xm.x)*y*y,
                        (Xp.x*Xp.x - Xm.x*Xm.x)*y*y,
                        (exp(Xp.x) - exp(Xm.x))*exp(2*y),
                        sin(x)*(sin(Xp.y) - sin(Xm.y)),
                        x*x*(Xp.y*Xp.y - Xm.y*Xm.y),
                        x*x*(Xp.y*Xp.y - Xm.y*Xm.y),
                        exp(x)*(exp(2*Xp.y) - exp(2*Xm.y))}.data());
      r.multiply(1.0 / eps);

      const utl::matrix3d<Real> grad = f.gradient(X);
      for (size_t d = 1; d <= 2; ++d) {
        const SymmTensor dx = f.deriv(X,d);
        for (size_t i = 1; i <= 2; ++i)
          for (size_t j = 1; j <= 2; ++j) {
            REQUIRE_THAT(dx(i,j), WithinRel(r(i,j,d), 1e-8));
            REQUIRE_THAT(grad(i,j,d), WithinRel(r(i,j,d), 1e-8));
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

  STensorFuncExprAd f2(g);

  REQUIRE(f1.isConstant());
  REQUIRE(f2.isConstant());

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
                REQUIRE_THAT(dx(i,j), WithinRel(r(i,j,d)));
                REQUIRE_THAT(grad(i,j,d), WithinRel(r(i,j,d)));
              }
          }
        }
      }
}


TEST_CASE("TestSTensorFunction.Gradient3DFD")
{
  const char* g   = "sin(x)*sin(y)*sin(z) | exp(x)*exp(2*y)*exp(z) | x*x*y*y*z*z |"
                    "x*y*z | x*x*y*z | z";

  const double eps = 1e-6;
  STensorFuncExpr f(g,"",eps);

  REQUIRE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const Vec3 Xp(x + 0.5*eps, y + 0.5*eps, z + 0.5*eps);
        const Vec3 Xm(x - 0.5*eps, y - 0.5*eps, z - 0.5*eps);
        utl::matrix3d<Real> r(3,3,3);
        r.fill(std::array{(sin(Xp.x) - sin(Xm.x))*sin(y)*sin(z),
                          (Xp.x - Xm.x)*y*z,
                          0.0,
                          (Xp.x - Xm.x)*y*z,
                          (exp(Xp.x) - exp(Xm.x))*exp(2*y)*exp(z),
                          (Xp.x*Xp.x - Xm.x*Xm.x)*y*z,
                          0.0,
                          (Xp.x*Xp.x - Xm.x*Xm.x)*y*z,
                          (Xp.x*Xp.x - Xm.x*Xm.x)*y*y*z*z,

                          sin(x)*(sin(Xp.y) - sin(Xm.y))*sin(z),
                          x*(Xp.y - Xm.y)*z,
                          0.0,
                          x*(Xp.y - Xm.y)*z,
                          exp(x)*(exp(2*Xp.y) - exp(2*Xm.y))*exp(z),
                          x*x*(Xp.y - Xm.y)*z,
                          0.0,
                          x*x*(Xp.y - Xm.y)*z,
                          x*x*(Xp.y*Xp.y - Xm.y*Xm.y)*z*z,

                          sin(x)*sin(y)*(sin(Xp.z) - sin(Xm.z)),
                          x*y*(Xp.z - Xm.z),
                          Xp.z - Xm.z,
                          x*y*(Xp.z - Xm.z),
                          exp(x)*exp(2*y)*(exp(Xp.z) -  exp(Xm.z)),
                          x*x*y*(Xp.z - Xm.z),
                          Xp.z - Xm.z,
                          x*x*y*(Xp.z - Xm.z),
                          x*x*y*y*(Xp.z*Xp.z - Xm.z*Xm.z)}.data());
        r.multiply(1.0 / eps);

        const utl::matrix3d<Real> grad = f.gradient(X);
        for (size_t d = 1; d <= 3; ++d) {
          const SymmTensor dx = f.deriv(X,d);
          for (size_t i = 1; i <= 3; ++i)
            for (size_t j = 1; j <= 3; ++j) {
              REQUIRE_THAT(dx(i,j), WithinRel(r(i,j,d), 1e-8));
              REQUIRE_THAT(grad(i,j,d), WithinRel(r(i,j,d), 1e-8));
            }
        }
      }
}


TEST_CASE("TestSTensorFuncExpr.NumDimensions")
{
  const char* func1 = "x";
  const char* func2 = "x | y | z";
  const char* func3 = "x | y | z | x | y | z";

  STensorFuncExpr f1(func1);
  REQUIRE(f1.getNoSpaceDim() == 1);
  REQUIRE(f1.dim() == 1);
  STensorFuncExpr f2(func2);
  REQUIRE(f2.getNoSpaceDim() == 2);
  REQUIRE(f2.dim() == 3);
  STensorFuncExpr f3(func3);
  REQUIRE(f3.getNoSpaceDim() == 3);
  REQUIRE(f3.dim() == 6);
}


TEST_CASE("TestSTensorFunction.Hessian2D")
{
  const char* g    = "sin(x)*sin(y) | exp(x)*exp(2*y) | x*x*y*y";
  const char* g_xx = "-sin(x)*sin(y) | exp(x)*exp(2*y) | 2*y*y";
  const char* g_xy = "cos(x)*cos(y) | exp(x)*2*exp(2*y) | 2*x*2*y";
  const char* g_yy = "-sin(x)*sin(y) | exp(x)*4*exp(2*y) | x*x*2";

  STensorFuncExpr f1(g);
  f1.addDerivative(g_xx,"",1,1);
  f1.addDerivative(g_xy,"",1,2);
  f1.addDerivative(g_yy,"",2,2);

  STensorFuncExprAd f2(g);

  REQUIRE(f1.isConstant());
  REQUIRE(f2.isConstant());

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
                REQUIRE_THAT(dx(i,j), WithinRel(r(i,j,d1,d2)));
                REQUIRE_THAT(hess(i,j,d1,d2), WithinRel(r(i,j,d1,d2)));
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

  STensorFuncExprAd f2(g);

  REQUIRE(f1.isConstant());
  REQUIRE(f2.isConstant());

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
                  REQUIRE_THAT(dx(i,j), WithinRel(r(i,j,d1,d2)));
                  REQUIRE_THAT(hess(i,j,d1,d2), WithinRel(r(i,j,d1,d2)));
                }
            }
        }
      }
 }


TEST_CASE("TestSTensorFunction.TimeDerivative")
{
  const char* g   = "sin(x)*sin(y)*sin(t) | exp(x)*exp(2*y)*exp(-4*t) | x*x*y*y*t*t";
  const char* g_t = "sin(x)*sin(y)*cos(t) | exp(x)*exp(2*y)*-4*exp(-4*t) | x*x*y*y*2*t";

  STensorFuncExpr f(g);
  f.addDerivative(g_t,"",4);

  REQUIRE(!f.isConstant());

  for (double t : {0.1, 0.2, 0.3})
    for (double x : {0.1, 0.2, 0.3})
      for (double y : {0.5, 0.6, 0.7}) {
        const Vec4 X(x,y,0,t);
        const SymmTensor r({sin(x)*sin(y)*cos(t), exp(x)*exp(2*y)*-4*exp(-4*t), x*x*y*y*2*t});

        const SymmTensor dt = f.timeDerivative(X);
        for (size_t i = 1; i <= 2; ++i)
          for (size_t j = 1; j <= 2; ++j)
            REQUIRE_THAT(dt(i,j), WithinRel(r(i,j)));
      }
}


TEST_CASE("TestSTensorFunction.TimeDerivativeFD")
{
  const char* g   = "sin(x)*sin(y)*sin(t) | exp(x)*exp(2*y)*exp(-4*t) | x*x*y*y*t*t";

  const double eps = 1e-6;
  STensorFuncExpr f(g,"",1e-8,eps);

  REQUIRE(!f.isConstant());

  for (double t : {0.1, 0.2, 0.3}) {
    const double tp = t + 0.5*eps;
    const double tm = t - 0.5*eps;
    for (double x : {0.1, 0.2, 0.3})
      for (double y : {0.5, 0.6, 0.7}) {
        const Vec4 X(x,y,0,t);
        SymmTensor r({sin(x)*sin(y)*(sin(tp) - sin(tm)),
                      exp(x)*exp(2*y)*(exp(-4*tp) - exp(-4*tm)),
                      x*x*y*y*(tp*tp - tm*tm)});
        r *= 1.0 / eps;

        const SymmTensor dt = f.timeDerivative(X);
        for (size_t i = 1; i <= 2; ++i)
          for (size_t j = 1; j <= 2; ++j)
            REQUIRE_THAT(dt(i,j), WithinRel(r(i,j), 1e-8));
      }
  }
}


TEST_CASE("TestEvalFunction.ExtraParam")
{
  EvalFunction f("x*foo");
  f.setParam("foo", 2.0);
  Vec3 X(1.0,0.0,0.0);
  REQUIRE_THAT(f(X), WithinRel(2.0));
  X.x = 0.5;
  f.setParam("foo", 4.0);
  REQUIRE_THAT(f(X), WithinRel(2.0));
}


TEST_CASE("TestEvalFunction.isConstant")
{
  REQUIRE(EvalFunction("2.0*x*y").isConstant());
  REQUIRE(!EvalFunction("2.0*x*t").isConstant());
  REQUIRE(EvalFunction("1.8*tan(x)*x").isConstant());
  REQUIRE(!EvalFunction("2.0*x*tan(t*3)+y").isConstant());
}


TEST_CASE("TestEvalFunction.Derivatives")
{
  const char* g    = "sin(x)*sin(y)*sin(z)*sin(t)";
  const char* g_x  = "cos(x)*sin(y)*sin(z)*sin(t)";
  const char* g_y  = "sin(x)*cos(y)*sin(z)*sin(t)";
  const char* g_z  = "sin(x)*sin(y)*cos(z)*sin(t)";
  const char* g_t  = "sin(x)*sin(y)*sin(z)*cos(t)";
  const char* g_xx = "-sin(x)*sin(y)*sin(z)*sin(t)";
  const char* g_xy = "cos(x)*cos(y)*sin(z)*sin(t)";
  const char* g_xz = "cos(x)*sin(y)*cos(z)*sin(t)";
  const char* g_yy = "-sin(x)*sin(y)*sin(z)*sin(t)";
  const char* g_yz = "sin(x)*cos(y)*cos(z)*sin(t)";
  const char* g_zz = "-sin(x)*sin(y)*sin(z)*sin(t)";

  EvalFunction f(g);
  f.addDerivative(g_x,"",1);
  f.addDerivative(g_y,"",2);
  f.addDerivative(g_z,"",3);
  f.addDerivative(g_t,"",4);
  f.addDerivative(g_xx,"",1,1);
  f.addDerivative(g_xy,"",1,2);
  f.addDerivative(g_xz,"",1,3);
  f.addDerivative(g_yy,"",2,2);
  f.addDerivative(g_yz,"",2,3);
  f.addDerivative(g_zz,"",3,3);

  for (const double t : {0.1, 0.2, 0.3})
    for (const double x : {0.4, 0.5, 0.6})
      for (const double y : {0.7, 0.8, 0.9})
        for (const double z : {1.0, 1.1, 1.2}) {
          const Vec4 X(x,y,z,t);
          REQUIRE_THAT(f(X),            WithinRel(sin(x)*sin(y)*sin(z)*sin(t)));
          REQUIRE_THAT(f.deriv(X,1),    WithinRel(cos(x)*sin(y)*sin(z)*sin(t)));
          REQUIRE_THAT(f.deriv(X,2),    WithinRel(sin(x)*cos(y)*sin(z)*sin(t)));
          REQUIRE_THAT(f.deriv(X,3),    WithinRel(sin(x)*sin(y)*cos(z)*sin(t)));
          REQUIRE_THAT(f.deriv(X,4),    WithinRel( sin(x)*sin(y)*sin(z)*cos(t)));
          REQUIRE_THAT(f.dderiv(X,1,1), WithinRel(-sin(x)*sin(y)*sin(z)*sin(t)));
          REQUIRE_THAT(f.dderiv(X,1,2), WithinRel(cos(x)*cos(y)*sin(z)*sin(t)));
          REQUIRE_THAT(f.dderiv(X,1,3), WithinRel(cos(x)*sin(y)*cos(z)*sin(t)));
          REQUIRE_THAT(f.dderiv(X,2,1), WithinRel(cos(x)*cos(y)*sin(z)*sin(t)));
          REQUIRE_THAT(f.dderiv(X,2,2), WithinRel(-sin(x)*sin(y)*sin(z)*sin(t)));
          REQUIRE_THAT(f.dderiv(X,2,3), WithinRel(sin(x)*cos(y)*cos(z)*sin(t)));
          REQUIRE_THAT(f.dderiv(X,3,1), WithinRel(cos(x)*sin(y)*cos(z)*sin(t)));
          REQUIRE_THAT(f.dderiv(X,3,2), WithinRel( sin(x)*cos(y)*cos(z)*sin(t)));
          REQUIRE_THAT(f.dderiv(X,3,3), WithinRel(-sin(x)*sin(y)*sin(z)*sin(t)));
        }
}
