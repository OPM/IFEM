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

#include <cstdlib>
#include <cmath>

#include "Catch2Support.h"


TEST_CASE("TestScalarFunc.ParseDerivative")
{
  const char* func1 = "sin(1.5*t)*t";
  const char* func2 = "sin(1.5*t)*t:1.5*cos(1.5*t)*t+sin(1.5*t)";

  ScalarFunc* f1 = utl::parseTimeFunc(func1,"expression");
  ScalarFunc* f2 = utl::parseTimeFunc(func2,"expression");

  REQUIRE(f1 != nullptr);
  REQUIRE(f2 != nullptr);

  CHECK(!f1->isConstant());
  CHECK(!f2->isConstant());

  double t = 0.0;
  for (int i = 0; i < 20; i++)
  {
    t += 0.314*(double)random()/(double)RAND_MAX;
    std::cout <<"f("<< t <<") = "<< (*f1)(t)
              <<"  f'("<< t <<") = "<< f1->deriv(t) << std::endl;
    CHECK_THAT((*f1)(t), WithinRel(sin(1.5*t)*t));
    CHECK_THAT((*f2)(t), WithinRel(sin(1.5*t)*t));
    CHECK_THAT(f1->deriv(t), WithinRel(1.5*cos(1.5*t)*t+sin(1.5*t), 1e-6));
    CHECK_THAT(f2->deriv(t), WithinRel(1.5*cos(1.5*t)*t+sin(1.5*t)));
  }
}


TEST_CASE("TestScalarFunc.Parse")
{
  std::cout <<"Parsing scalar function: ";
  ScalarFunc* f1 = utl::parseTimeFunc("1.2 100.0","Dirac");
  std::cout <<"Parsing scalar function: ";
  ScalarFunc* f2 = utl::parseTimeFunc("5.0 100.0","Ramp");
  std::cout <<"Parsing scalar function: ";
  ScalarFunc* f3 = utl::parseTimeFunc("100.0 0.3 1.0 2.0 3.5 5.0 8.0","Pulse");
  std::cout <<"Parsing scalar function: ";
  ScalarFunc* f4 = utl::parseTimeFunc("100.0 1.0 0.5 2.0 0.8 5.0 0.1","Pulses");

  REQUIRE(f1 != nullptr);
  REQUIRE(f2 != nullptr);
  REQUIRE(f3 != nullptr);
  REQUIRE(f4 != nullptr);

  CHECK_THAT((*f1)(1.1), WithinAbs(0.0, 1e-12));
  CHECK_THAT((*f1)(1.2), WithinRel(100.0));
  CHECK_THAT((*f1)(1.3), WithinAbs(0.0, 1e-12));
  CHECK_THAT((*f2)(2.5), WithinRel(50.0));
  CHECK_THAT((*f2)(5.0), WithinRel(100.0));
  CHECK_THAT((*f2)(7.0), WithinRel(100.0));
  CHECK_THAT(f2->integrate(10.0), WithinRel(750.0));

  CHECK_THAT((*f3)(0.8), WithinRel(0.0));
  CHECK_THAT((*f3)(1.2), WithinRel(100.0));
  CHECK_THAT((*f3)(1.6), WithinRel(0.0));
  CHECK_THAT((*f3)(2.1), WithinRel(100.0));
  CHECK_THAT((*f3)(3.0), WithinRel(0.0));
  CHECK_THAT((*f3)(3.7), WithinRel(100.0));
  CHECK_THAT((*f3)(5.1), WithinRel(100.0));
  CHECK_THAT((*f3)(5.4), WithinRel(0.0));
  CHECK_THAT(f3->integrate(10.0), WithinRel(150.0));

  CHECK_THAT((*f4)(0.8), WithinRel(0.0));
  CHECK_THAT((*f4)(1.2), WithinRel(100.0));
  CHECK_THAT((*f4)(1.6), WithinRel(0.0));
  CHECK_THAT((*f4)(2.6), WithinRel(100.0));
  CHECK_THAT((*f4)(3.0), WithinRel(0.0));
  CHECK_THAT((*f4)(5.1), WithinRel(100.0));
  CHECK_THAT((*f4)(5.2), WithinRel(0.0));
  CHECK_THAT(f4->integrate(10.0), WithinRel(140.0));
}


TEST_CASE("TestSineFunc.Integrate")
{
  SineFunc f(10.0,0.2,0.5);

  auto integrate = [&f](double x)
  {
    double fInt = f.integrate(x);
    std::cout <<"f.integrate(0,"<< x <<") = "<< fInt << std::endl;
    return fInt;
  };

  CHECK_THAT(integrate(1.0), WithinAbs(f.ScalarFunc::integrate(0.0,1.0), 1e-6));
  CHECK_THAT(integrate(3.5), WithinAbs(f.ScalarFunc::integrate(0.0,3.5), 1e-6));
  CHECK_THAT(integrate(4.2), WithinAbs(f.ScalarFunc::integrate(0.0,4.2), 1e-6));
}


TEST_CASE("TestLinearFunc.Integrate")
{
  LinearFunc f({ 1.0, 2.0, 3.0, 4.0 },
               { 0.0, 1.0, 0.5, 1.2 });

  auto integrate = [&f](double x0, double x1)
  {
    constexpr double dx = 0.001;
    double fInt = f.integrate(x0,x1,dx);
    std::cout <<"f.integrate("<< x0 <<","<< x1 <<") = "<< fInt << std::endl;
    return fInt;
  };

  CHECK_THAT(integrate(0.0,5.0), WithinAbs(f.ScalarFunc::integrate(0.0,5.0), 1e-12));
  CHECK_THAT(integrate(1.2,1.5), WithinAbs(f.ScalarFunc::integrate(1.2,1.5), 1e-12));
  CHECK_THAT(integrate(1.5,2.6), WithinAbs(f.ScalarFunc::integrate(1.5,2.6), 1e-12));
  CHECK_THAT(integrate(1.8,3.3), WithinAbs(f.ScalarFunc::integrate(1.8,3.3), 1e-12));
  CHECK_THAT(integrate(2.2,4.5), WithinAbs(f.ScalarFunc::integrate(2.2,4.5), 1e-12));
}


TEST_CASE("TestRealFunc.Parse")
{
  std::cout <<"Parsing real function";
  RealFunc* f1 = utl::parseRealFunc("100.0 1.2","Dirac");
  std::cout <<"\nParsing real function";
  RealFunc* f2 = utl::parseRealFunc("100.0 5.0","Ramp",false);
  std::cout << std::endl;

  REQUIRE(f1 != nullptr);
  REQUIRE(f2 != nullptr);

  CHECK_THAT((*f1)(Vec4(1.1)), WithinAbs(0.0, 1e-12));
  CHECK_THAT((*f1)(Vec4(1.2)), WithinRel(100.0));
  CHECK_THAT((*f1)(Vec4(1.3)), WithinAbs(0.0, 1e-12));
  CHECK_THAT((*f2)(Vec4(2.5)), WithinRel(50.0));
  CHECK_THAT((*f2)(Vec4(5.0)), WithinRel(100.0));
  CHECK_THAT((*f2)(Vec4(7.0)), WithinRel(100.0));
}


TEST_CASE("TestTractionFunc.Parse")
{
  std::cout <<"Parsing traction function";
  TractionFunc* f1 = utl::parseTracFunc("1.0+x*x-y","expression");
  std::cout << std::endl;
  REQUIRE(f1 != nullptr);

  class VecX : public Vec3
  {
  public:
    explicit VecX(double x) : Vec3(x,0.0,0.0) {}
  };

  const VecX n(1.0);
  CHECK((*f1)(VecX(0.0),n).equal(VecX(1.00)));
  CHECK((*f1)(VecX(0.5),n).equal(VecX(1.25)));
  CHECK((*f1)(Vec3(0.7,0.3),n).equal(VecX(1.19)));
}


TEST_CASE("TestRealFunc.GradientFD")
{
  const char* g = "sin(x)*sin(y)*sin(z)";

  constexpr double eps = 1.0e-6;

  EvalFunction f(g,eps);

  CHECK(f.isConstant());

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
          CHECK_THAT(f.deriv(X, i), WithinRel(r[i-1], 1e-8));
          CHECK_THAT(grad[i-1], WithinRel(r[i-1], 1e-8));
        }
      }
}


TEST_CASE("TestVecFunction.Gradient2DFD")
{
  const char* g = "sin(x)*sin(y) | x*x*y*y";

  constexpr double eps = 1.0e-6;

  VecFuncExpr f(g,"",eps);

  CHECK(f.isConstant());

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
          CHECK_THAT(dx[i-1], WithinRel(r(i,d), 1e-8));
          CHECK_THAT(grad(i,d), WithinRel(r(i,d), 1e-8));
        }
      }
    }
}


TEST_CASE("TestVecFunction.Gradient3DFD")
{
  const char* g = "sin(x)*sin(y)*sin(z) | x*x*y*y*z*z*z | exp(-x)*exp(2*y)*exp(z*z)";

  constexpr double eps = 1.0e-6;

  VecFuncExpr f(g,"",eps);

  CHECK(f.isConstant());

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
            CHECK_THAT(dx[i-1], WithinRel(r(i,d), 1e-8));
            CHECK_THAT(grad(i,d), WithinRel(r(i,d), 1e-8));
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
  CHECK(f1.getNoSpaceDim() == 1);
  CHECK(f1.dim() ==  1);
  VecFuncExpr f2(func2);
  CHECK(f2.getNoSpaceDim() == 2);
  CHECK(f2.dim() == 2);
  VecFuncExpr f3(func3);
  CHECK(f3.getNoSpaceDim() == 3);
  CHECK(f3.dim() == 3);
}


TEST_CASE("TestVecFuncExpr.TimeDerivative")
{
  const char* g = "sin(x)*sin(y)*sin(z)*sin(t) | x*x*y*y*z*t*t | exp(-2*t)";
  const char* d = "sin(x)*sin(y)*sin(z)*cos(t) | x*x*y*y*z*2*t | -2*exp(-2*t)";

  VecFuncExpr f(g);
  f.addDerivative(d,"",4);

  CHECK(!f.isConstant());

  for (double t : {0.1, 0.2, 0.3})
    for (double x : {0.1, 0.2, 0.3})
      for (double y : {0.5, 0.6, 0.7})
        for (double z : {0.8, 0.9, 1.0}) {
          const Vec4 X(x,y,z,t);
          const Vec3 r(sin(x)*sin(y)*sin(z)*cos(t),
                       x*x*y*y*z*2*t,
                       -2*exp(-2*t));
          const Vec3 dt = f.timeDerivative(X);
          CHECK_THAT(dt[0], WithinRel(r[0]));
          CHECK_THAT(dt[1], WithinRel(r[1]));
          CHECK_THAT(dt[2], WithinRel(r[2]));
        }
}


TEST_CASE("TestVecFuncExpr.TimeDerivativeFD")
{
  const char* g = "sin(x)*sin(y)*sin(z)*sin(t) | x*x*y*y*z*t*t | exp(-2*t)";

  constexpr double eps = 1.0e-6;

  VecFuncExpr f(g,"",1.0e-8,eps);

  CHECK(!f.isConstant());

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
          CHECK_THAT(dt[0], WithinRel(r[0], 1e-8));
          CHECK_THAT(dt[1], WithinRel(r[1], 1e-8));
          CHECK_THAT(dt[2], WithinRel(r[2], 1e-8));
        }
  }
}


TEST_CASE("TestTensorFunction.Gradient2DFD")
{
  const char* g = "sin(x)*sin(y) | x*x*y*y | exp(x)*exp(2*y) | exp(-2*x)*exp(y)";

  constexpr double eps = 1.0e-6;

  TensorFuncExpr f(g,"",eps);

  CHECK(f.isConstant());

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
            CHECK_THAT(dx(i,j), WithinRel(r(i,j,d), 1e-8));
            CHECK_THAT(grad(i,j,d), WithinRel(r(i,j,d), 1e-8));
          }
      }
    }
}


TEST_CASE("TestTensorFunction.Gradient3DFD")
{
  const char* g = "sin(x)*sin(y)*sin(z) | x*x*y*y*z | exp(x)*exp(2*y)*z*z |"
                  "exp(-2*x)*exp(y)*z   | x*y*z     | x*y*z*z |"
                  "x                    | y         | z";

  constexpr double eps = 1.0e-6;

  TensorFuncExpr f(g,"",eps);

  CHECK(f.isConstant());

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
              CHECK_THAT(dx(i,j), WithinRel(r(i,j,d), 1e-8));
              CHECK_THAT(grad(i,j,d), WithinRel(r(i,j,d), 1e-8));
            }
        }
    }
}


TEST_CASE("TestTensorFunction.TimeDerivative")
{
  const char* g = "sin(x)*sin(y)*sin(z)*sin(t) | x*x*y*y*z*t*t | exp(x)*exp(2*y)*z*z*t |"
                  "exp(-2*x)*exp(y)*z*sin(t)   | x*y*z*t*t     | x*y*z*z*t |"
                  "x*sin(t)                    | y*t*t         | z*t";
  const char* d = "sin(x)*sin(y)*sin(z)*cos(t) | x*x*y*y*z*2*t | exp(x)*exp(2*y)*z*z |"
                  "exp(-2*x)*exp(y)*z*cos(t)   | x*y*z*2*t     | x*y*z*z |"
                  "x*cos(t)                    | y*2*t         | z";

  TensorFuncExpr f(g);
  f.addDerivative(d,"",4);

  CHECK(!f.isConstant());

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
              CHECK_THAT(dt(i,j), WithinRel(r(i,j)));
      }
}


TEST_CASE("TestTensorFuncExpr.NumDimensions")
{
  const char* func1 = "x";
  const char* func2 = "x | y | z | x";
  const char* func3 = "x | y | z | x | y | z | x | y | z";

  TensorFuncExpr f1(func1);
  CHECK(f1.getNoSpaceDim() == 1);
  CHECK(f1.dim() == 1);
  TensorFuncExpr f2(func2);
  CHECK(f2.getNoSpaceDim() == 2);
  CHECK(f2.dim() == 4);
  TensorFuncExpr f3(func3);
  CHECK(f3.getNoSpaceDim() == 3);
  CHECK(f3.dim() == 9);
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
            CHECK_THAT(fx(i,j), WithinRel(r(i,j)));
        CHECK_THAT(fx(3,3), WithinRel(sin(x)*cos(y)));
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

  CHECK(f.isConstant());

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
            CHECK_THAT(dx(i,j), WithinRel(r(i,j,d), 1e-12));
        CHECK_THAT(dx(3,3), WithinRel((d == 1 ? cos(x) : sin(x)) * (d == 2 ? cos(y) : sin(y)), 1e-12));
      }
    }
}


TEST_CASE("TestSTensorFunction.Gradient2DFD")
{
  const char* g = "sin(x)*sin(y) | exp(x)*exp(2*y) | x*x*y*y";

  constexpr double eps = 1.0e-6;

  STensorFuncExpr f(g,"",eps);

  CHECK(f.isConstant());

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
            CHECK_THAT(dx(i,j), WithinRel(r(i,j,d), 1e-8));
            CHECK_THAT(grad(i,j,d), WithinRel(r(i,j,d), 1e-8));
          }
      }
    }
}


TEST_CASE("TestSTensorFunction.Gradient3DFD")
{
  const char* g = "sin(x)*sin(y)*sin(z) | exp(x)*exp(2*y)*exp(z) | x*x*y*y*z*z |"
                  "x*y*z | x*x*y*z | z";

  constexpr double eps = 1e-6;

  STensorFuncExpr f(g,"",eps);

  CHECK(f.isConstant());

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
              CHECK_THAT(dx(i,j), WithinRel(r(i,j,d), 1e-8));
              CHECK_THAT(grad(i,j,d), WithinRel(r(i,j,d), 1e-8));
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
  CHECK(f1.getNoSpaceDim() == 1);
  CHECK(f1.dim() == 1);
  STensorFuncExpr f2(func2);
  CHECK(f2.getNoSpaceDim() == 2);
  CHECK(f2.dim() == 3);
  STensorFuncExpr f3(func3);
  CHECK(f3.getNoSpaceDim() == 3);
  CHECK(f3.dim() == 6);
}


TEST_CASE("TestSTensorFunction.TimeDerivative")
{
  const char* g = "sin(x)*sin(y)*sin(t) | exp(x)*exp(2*y)*exp(-4*t) | x*x*y*y*t*t";
  const char* d = "sin(x)*sin(y)*cos(t) | exp(x)*exp(2*y)*-4*exp(-4*t) | x*x*y*y*2*t";

  STensorFuncExpr f(g);
  f.addDerivative(d,"",4);

  CHECK(!f.isConstant());

  for (double t : {0.1, 0.2, 0.3})
    for (double x : {0.1, 0.2, 0.3})
      for (double y : {0.5, 0.6, 0.7}) {
        const Vec4 X(x,y,0,t);
        const SymmTensor r({sin(x)*sin(y)*cos(t), exp(x)*exp(2*y)*-4*exp(-4*t), x*x*y*y*2*t});

        const SymmTensor dt = f.timeDerivative(X);
        for (size_t i = 1; i <= 2; ++i)
          for (size_t j = 1; j <= 2; ++j)
            CHECK_THAT(dt(i,j), WithinRel(r(i,j)));
      }
}


TEST_CASE("TestSTensorFunction.TimeDerivativeFD")
{
  const char* g = "sin(x)*sin(y)*sin(t) | exp(x)*exp(2*y)*exp(-4*t) | x*x*y*y*t*t";

  constexpr double eps = 1.0e-6;

  STensorFuncExpr f(g,"",1.0e-8,eps);

  CHECK(!f.isConstant());

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
            CHECK_THAT(dt(i,j), WithinRel(r(i,j), 1e-8));
      }
  }
}


TEST_CASE("TestEvalFunction.ExtraParam")
{
  EvalFunction f("x*foo");
  f.setParam("foo", 2.0);
  Vec3 X(1.0,0.0,0.0);
  CHECK_THAT(f(X), WithinRel(2.0));
  X.x = 0.5;
  f.setParam("foo", 4.0);
  CHECK_THAT(f(X), WithinRel(2.0));
}


TEST_CASE("TestEvalFunction.isConstant")
{
  CHECK(EvalFunction("2.0*x*y").isConstant());
  CHECK(!EvalFunction("2.0*x*t").isConstant());
  CHECK(EvalFunction("1.8*tan(x)*x").isConstant());
  CHECK(!EvalFunction("2.0*x*tan(t*3)+y").isConstant());
}


TEST_CASE("TestEvalFunction.isTime")
{
  CHECK(!utl::isTimeExpression("2.0*x*y"));
  CHECK(utl::isTimeExpression("2.0*x*y*sin(t)"));
  CHECK(utl::isTimeExpression("2*t"));
  CHECK(!utl::isTimeExpression("time0*x"));
  CHECK(!utl::isTimeExpression("hat*x+y"));
  CHECK(utl::isTimeExpression("20+5*t"));
  CHECK(!utl::isTimeExpression("3+t_0"));
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
          CHECK_THAT(f(X),            WithinRel(sin(x)*sin(y)*sin(z)*sin(t)));
          CHECK_THAT(f.deriv(X,1),    WithinRel(cos(x)*sin(y)*sin(z)*sin(t)));
          CHECK_THAT(f.deriv(X,2),    WithinRel(sin(x)*cos(y)*sin(z)*sin(t)));
          CHECK_THAT(f.deriv(X,3),    WithinRel(sin(x)*sin(y)*cos(z)*sin(t)));
          CHECK_THAT(f.deriv(X,4),    WithinRel( sin(x)*sin(y)*sin(z)*cos(t)));
          CHECK_THAT(f.dderiv(X,1,1), WithinRel(-sin(x)*sin(y)*sin(z)*sin(t)));
          CHECK_THAT(f.dderiv(X,1,2), WithinRel(cos(x)*cos(y)*sin(z)*sin(t)));
          CHECK_THAT(f.dderiv(X,1,3), WithinRel(cos(x)*sin(y)*cos(z)*sin(t)));
          CHECK_THAT(f.dderiv(X,2,1), WithinRel(cos(x)*cos(y)*sin(z)*sin(t)));
          CHECK_THAT(f.dderiv(X,2,2), WithinRel(-sin(x)*sin(y)*sin(z)*sin(t)));
          CHECK_THAT(f.dderiv(X,2,3), WithinRel(sin(x)*cos(y)*cos(z)*sin(t)));
          CHECK_THAT(f.dderiv(X,3,1), WithinRel(cos(x)*sin(y)*cos(z)*sin(t)));
          CHECK_THAT(f.dderiv(X,3,2), WithinRel( sin(x)*cos(y)*cos(z)*sin(t)));
          CHECK_THAT(f.dderiv(X,3,3), WithinRel(-sin(x)*sin(y)*sin(z)*sin(t)));
        }
}
