//==============================================================================
//!
//! \file TestChebyshev.C
//!
//! \date Dec 6 2023
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for Chebyshev functions.
//!
//==============================================================================

#include "Chebyshev.h"

#include "Catch2Support.h"


namespace {

#include "refdata/cheby_strings.h"

}


TEST_CASE("TestChebyshevFunc.Value1D")
{
  ChebyshevFunc cheb(cheby_1, false);

  auto f = [](Real x) { return x*x*x + x*x + x - 1.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0}) {
    const Vec4 Xc({x,0.0, 0.0, 0.0}, &x);
    REQUIRE_THAT(cheb(Xc), WithinRel(f(x), 1e-12));
  }
}


TEST_CASE("TestChebyshevFunc.Gradient1D")
{
  ChebyshevFunc cheb(cheby_1, false);

  auto df = [](Real x) { return 3.0*x*x + 2.0*x + 1.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0}) {
    const Vec4 Xc({x,0.0, 0.0, 0.0}, &x);
    const Vec3 dx = cheb.gradient(Xc);
    REQUIRE_THAT(dx[0], WithinRel(df(x), 1e-12));
    REQUIRE_THAT(dx[1], WithinRel(0.0));
    REQUIRE_THAT(dx[2], WithinRel(0.0));
  }
}


TEST_CASE("TestChebyshevFunc.Hessian1D")
{
  ChebyshevFunc cheb(cheby_1, false);

  auto d2f = [](Real x) { return 6.0*x + 2.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0}) {
    const Vec4 Xc({x,0.0, 0.0, 0.0}, &x);
    const SymmTensor d2x = cheb.hessian(Xc);
    REQUIRE_THAT(d2x(1,1), WithinRel(d2f(x), 1e-12));
    REQUIRE_THAT(d2x(1,2), WithinRel(0.0));
    REQUIRE_THAT(d2x(1,3), WithinRel(0.0));
    REQUIRE_THAT(d2x(2,2), WithinRel(0.0));
    REQUIRE_THAT(d2x(2,3), WithinRel(0.0));
    REQUIRE_THAT(d2x(3,3), WithinRel(0.0));
  }
}


TEST_CASE("TestChebyshevFunc.Value2D")
{
  ChebyshevFunc cheb(cheby_2, false);

  auto f = [](Real x) { return x*x*x + x*x + x - 1.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88}) {
      const double p[] {x,y};
      const Vec4 Xc({x, y, 0.0, 0.0}, p);
      REQUIRE_THAT(cheb(Xc), WithinRel(f(x) * f(y), 1e-12));
    }
}


TEST_CASE("TestChebyshevFunc.Gradient2D")
{
  ChebyshevFunc cheb(cheby_2, false);

  auto  f = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto df = [](Real x) { return 3.0*x*x + 2.0*x + 1.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88}) {
      const double p[] {x,y};
      const Vec4 Xc({x, y, 0.0, 0.0}, p);
      const Vec3 dx = cheb.gradient(Xc);
      REQUIRE_THAT(dx[0], WithinRel(df(x) *  f(y), 1e-12));
      REQUIRE_THAT(dx[1], WithinRel(f(x) * df(y), 1e-12));
      REQUIRE_THAT(dx[2], WithinRel(0.0));
    }
}


TEST_CASE("TestChebyshevFunc.Hessian2D")
{
  ChebyshevFunc cheb(cheby_2, false);

  auto   f = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto  df = [](Real x) { return 3.0*x*x + 2.0*x + 1.0; };
  auto d2f = [](Real x) { return 6.0*x + 2.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88}) {
      const double p[] {x,y};
      const Vec4 Xc({x, y, 0.0, 0.0}, p);
      const SymmTensor d2x = cheb.hessian(Xc);
      REQUIRE_THAT(d2x(1,1), WithinRel(d2f(x)*f(y), 1e-12));
      REQUIRE_THAT(d2x(1,2), WithinRel(df(x)*df(y), 1e-12));
      REQUIRE_THAT(d2x(1,3), WithinRel(0.0));
      REQUIRE_THAT(d2x(2,2), WithinRel(f(x)*d2f(y), 1e-12));
      REQUIRE_THAT(d2x(2,3), WithinRel(0.0));
      REQUIRE_THAT(d2x(3,3), WithinRel(0.0));
    }
}


TEST_CASE("TestChebyshevFunc.Value3D")
{
  ChebyshevFunc cheb(cheby_3, false);

  auto f = [](Real x) { return x*x*x + x*x + x - 1.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88})
      for (double z : {0.2, 0.4, 0.6, 0.77}) {
        const double p[] {x,y,z};
        const Vec4 Xc({x, y, z, 0.0}, p);
        REQUIRE_THAT(cheb(Xc), WithinRel(f(x) * f(y) * f(z), 1e-12));
      }
}


TEST_CASE("TestChebyshevFunc3D.Gradient3D")
{
  ChebyshevFunc cheb(cheby_3, false);

  auto  f = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto df = [](Real x) { return 3.0*x*x + 2.0*x + 1.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88})
      for (double z : {0.2, 0.4, 0.6, 0.77}) {
        const double p[] {x,y,z};
        const Vec4 Xc({x, y, z, 0.0}, p);
        const Vec3 dx = cheb.gradient(Xc);
        REQUIRE_THAT(dx[0], WithinRel(df(x) *  f(y) *  f(z), 1e-12));
        REQUIRE_THAT(dx[1], WithinRel(f(x) * df(y) *  f(z), 1e-12));
        REQUIRE_THAT(dx[2], WithinRel(f(x) *  f(y) * df(z), 1e-12));
      }
}


TEST_CASE("TestChebyshevFunc.Hessian3D")
{
  ChebyshevFunc cheb(cheby_3, false);

  auto   f = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto  df = [](Real x) { return 3.0*x*x + 2.0*x + 1.0; };
  auto d2f = [](Real x) { return 6.0*x + 2.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88})
      for (double z : {0.2, 0.4, 0.6, 0.77}) {
        const double p[] {x,y,z};
        const Vec4 Xc({x, y, z, 0.0}, p);
        const SymmTensor d2x = cheb.hessian(Xc);
        REQUIRE_THAT(d2x(1,1), WithinRel(d2f(x) * f(y) * f(z), 1e-12));
        REQUIRE_THAT(d2x(1,2), WithinRel(df(x) * df(y) * f(z), 1e-12));
        REQUIRE_THAT(d2x(1,3), WithinRel(df(x) * f(y) * df(z), 1e-12));
        REQUIRE_THAT(d2x(2,2), WithinRel(f(x) * d2f(y) * f(z), 1e-12));
        REQUIRE_THAT(d2x(2,3), WithinRel(f(x) * df(y) * df(z), 1e-12));
        REQUIRE_THAT(d2x(3,3), WithinRel(f(x) * f(y) * d2f(z), 1e-12));
      }
}


TEST_CASE("TestChebyshevVecFunc.Value2D")
{
  ChebyshevVecFunc cheb({cheby_12, cheby_21}, false);

  auto f1 = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto f2 = [](Real x) { return x*x*x - 2.0*x*x + 3.0*x - 1.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88}) {
      const double p[] {x,y};
      const Vec4 Xc({x, y, 0.0, 0.0}, p);
      const Vec3 res = cheb(Xc);
      REQUIRE_THAT(res[0], WithinRel(f1(x) * f2(y), 1e-12));
      REQUIRE_THAT(res[1], WithinRel(f1(y) * f2(x), 1e-12));
      REQUIRE_THAT(res[2], WithinRel(0.0));
    }
}


TEST_CASE("TestChebyshevVecFunc.Gradient2D")
{
  ChebyshevVecFunc cheb({cheby_12, cheby_21}, false);

  auto  f1 = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto df1 = [](Real x) { return 3.0*x*x + 2.0*x + 1.0; };
  auto  f2 = [](Real x) { return x*x*x - 2.0*x*x + 3.0*x - 1.0; };
  auto df2 = [](Real x) { return 3.0*x*x - 4.0*x + 3.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88}) {
      const double p[] {x,y};
      const Vec4 Xc({x, y, 0.0, 0.0}, p);
      const Tensor res = cheb.gradient(Xc);
      REQUIRE_THAT(res(1,1), WithinRel(df1(x) *  f2(y), 1e-12));
      REQUIRE_THAT(res(1,2), WithinRel(f1(x) * df2(y), 1e-12));
      REQUIRE_THAT(res(2,1), WithinRel(df2(x) *  f1(y), 1e-12));
      REQUIRE_THAT(res(2,2), WithinRel(f2(x) * df1(y), 1e-12));
    }
}


TEST_CASE("TestChebyshevVecFunc.Hessian2D")
{
  ChebyshevVecFunc cheb({cheby_12, cheby_21}, false);

  auto   f1 = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto  df1 = [](Real x) { return 3.0*x*x + 2.0*x + 1.0; };
  auto d2f1 = [](Real x) { return 6*x + 2.0; };
  auto   f2 = [](Real x) { return x*x*x - 2.0*x*x + 3.0*x - 1.0; };
  auto  df2 = [](Real x) { return 3.0*x*x - 4.0*x + 3.0; };
  auto d2f2 = [](Real x) { return 6.0*x - 4.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88}) {
      const double p[] {x,y};
      const Vec4 Xc({x, y, 0.0, 0.0}, p);
      const utl::matrix3d<Real> res = cheb.hessian(Xc);
      REQUIRE_THAT(res(1,1,1), WithinRel(d2f1(x) *   f2(y), 1e-12));
      REQUIRE_THAT(res(1,1,2), WithinRel( df1(x) *  df2(y), 1e-12));
      REQUIRE_THAT(res(1,2,1), WithinRel( df1(x) *  df2(y), 1e-12));
      REQUIRE_THAT(res(1,2,2), WithinRel(  f1(x) * d2f2(y), 1e-12));
      REQUIRE_THAT(res(2,1,1), WithinRel(d2f2(x) *   f1(y), 1e-12));
      REQUIRE_THAT(res(2,1,2), WithinRel( df2(x) *  df1(y), 1e-12));
      REQUIRE_THAT(res(2,2,1), WithinRel( df2(x) *  df1(y), 1e-12));
      REQUIRE_THAT(res(2,2,2), WithinRel(  f2(x) * d2f1(y), 1e-12));
    }
}


TEST_CASE("TestChebyshevVecFunc.Value3D")
{
  ChebyshevVecFunc cheb({cheby_123, cheby_132, cheby_213}, false);

  auto f1 = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto f2 = [](Real x) { return x*x*x - 2.0*x*x + 3.0*x - 1.0; };
  auto f3 = [](Real x) { return 4.0*x*x*x - 7.0*x*x - 3.0*x + 5.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88})
      for (double z : {0.2, 0.4, 0.6, 0.77}) {
        const double p[] {x,y,z};
        const Vec4 Xc({x, y, z, 0.0}, p);
        const Vec3 res = cheb(Xc);
        REQUIRE_THAT(res[0], WithinRel(f1(x) * f2(y) * f3(z), 1e-12));
        REQUIRE_THAT(res[1], WithinRel(f1(x) * f2(z) * f3(y), 1e-12));
        REQUIRE_THAT(res[2], WithinRel(f2(x) * f1(y) * f3(z), 1e-12));
      }
}


TEST_CASE("TestChebyshevVecFunc.Gradient3D")
{
  ChebyshevVecFunc cheb({cheby_123, cheby_132, cheby_213}, false);

  auto  f1 = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto df1 = [](Real x) { return 3.0*x*x + 2.0*x + 1.0; };
  auto  f2 = [](Real x) { return x*x*x - 2.0*x*x + 3.0*x - 1.0; };
  auto df2 = [](Real x) { return 3.0*x*x - 4.0*x + 3.0; };
  auto  f3 = [](Real x) { return 4.0*x*x*x - 7.0*x*x - 3.0*x + 5.0; };
  auto df3 = [](Real x) { return 12.0*x*x - 14.0*x - 3.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88})
      for (double z : {0.2, 0.4, 0.6, 0.77}) {
        const double p[] {x,y,z};
        const Vec4 Xc({x, y, z, 0.0}, p);
        const Tensor res = cheb.gradient(Xc);
        REQUIRE_THAT(res(1,1), WithinRel(df1(x) *  f2(y) *  f3(z), 1e-12));
        REQUIRE_THAT(res(1,2), WithinRel( f1(x) * df2(y) *  f3(z), 1e-12));
        REQUIRE_THAT(res(1,3), WithinRel( f1(x) *  f2(y) * df3(z), 1e-12));
        REQUIRE_THAT(res(2,1), WithinRel(df1(x) *  f3(y) *  f2(z), 1e-12));
        REQUIRE_THAT(res(2,2), WithinRel( f1(x) * df3(y) *  f2(z), 1e-12));
        REQUIRE_THAT(res(2,3), WithinRel( f1(x) *  f3(y) * df2(z), 1e-12));
        REQUIRE_THAT(res(3,1), WithinRel(df2(x) *  f1(y) *  f3(z), 1e-12));
        REQUIRE_THAT(res(3,2), WithinRel( f2(x) * df1(y) *  f3(z), 1e-12));
        REQUIRE_THAT(res(3,3), WithinRel( f2(x) *  f1(y) * df3(z), 1e-12));
      }
}


TEST_CASE("TestChebyshevVecFunc.Hessian3D")
{
  ChebyshevVecFunc cheb({cheby_123, cheby_132, cheby_213}, false);

  auto   f1 = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto  df1 = [](Real x) { return 3.0*x*x + 2.0*x + 1.0; };
  auto d2f1 = [](Real x) { return 6.0*x + 2.0; };
  auto   f2 = [](Real x) { return x*x*x - 2.0*x*x + 3.0*x - 1.0; };
  auto  df2 = [](Real x) { return 3.0*x*x - 4.0*x + 3.0; };
  auto d2f2 = [](Real x) { return 6.0*x - 4.0; };
  auto   f3 = [](Real x) { return 4.0*x*x*x - 7.0*x*x - 3.0*x + 5.0; };
  auto  df3 = [](Real x) { return 12*x*x - 14.0*x - 3.0; };
  auto d2f3 = [](Real x) { return 24.0*x - 14.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88})
      for (double z : {0.2, 0.4, 0.6, 0.77}) {
        const double p[] {x,y,z};
        const Vec4 Xc({x, y, z, 0.0}, p);
        const utl::matrix3d<Real> res = cheb.hessian(Xc);
        REQUIRE_THAT(res(1,1,1), WithinRel(d2f1(x) *   f2(y) *   f3(z), 3e-12));
        REQUIRE_THAT(res(1,1,2), WithinRel( df1(x) *  df2(y) *   f3(z), 3e-12));
        REQUIRE_THAT(res(1,1,3), WithinRel( df1(x) *   f2(y) *  df3(z), 3e-12));
        REQUIRE_THAT(res(1,2,1), WithinRel( df1(x) *  df2(y) *   f3(z), 3e-12));
        REQUIRE_THAT(res(1,2,2), WithinRel(  f1(x) * d2f2(y) *   f3(z), 3e-12));
        REQUIRE_THAT(res(1,2,3), WithinRel(  f1(x) *  df2(y) *  df3(z), 3e-12));
        REQUIRE_THAT(res(1,3,1), WithinRel( df1(x) *   f2(y) *  df3(z), 3e-12));
        REQUIRE_THAT(res(1,3,2), WithinRel(  f1(x) *  df2(y) *  df3(z), 3e-12));
        REQUIRE_THAT(res(1,3,3), WithinRel(  f1(x) *   f2(y) * d2f3(z), 3e-12));

        REQUIRE_THAT(res(2,1,1), WithinRel(d2f1(x) *   f3(y) *   f2(z), 3e-12));
        REQUIRE_THAT(res(2,1,2), WithinRel( df1(x) *  df3(y) *   f2(z), 3e-12));
        REQUIRE_THAT(res(2,2,1), WithinRel( df1(x) *  df3(y) *   f2(z), 3e-12));
        REQUIRE_THAT(res(2,1,3), WithinRel( df1(x) *   f3(y) *  df2(z), 3e-12));
        REQUIRE_THAT(res(2,2,1), WithinRel( df1(x) *  df3(y) *   f2(z), 3e-12));
        REQUIRE_THAT(res(2,2,2), WithinRel(  f1(x) * d2f3(y) *   f2(z), 3e-12));
        REQUIRE_THAT(res(2,2,3), WithinRel(  f1(x) *  df3(y) *  df2(z), 3e-12));
        REQUIRE_THAT(res(2,3,1), WithinRel( df1(x) *   f3(y) *  df2(z), 3e-12));
        REQUIRE_THAT(res(2,3,2), WithinRel(  f1(x) *  df3(y) *  df2(z), 3e-12));
        REQUIRE_THAT(res(2,3,3), WithinRel(  f1(x) *   f3(y) * d2f2(z), 3e-12));

        REQUIRE_THAT(res(3,1,1), WithinRel(d2f2(x) *   f1(y) *   f3(z), 3e-12));
        REQUIRE_THAT(res(3,1,2), WithinRel( df2(x) *  df1(y) *   f3(z), 3e-12));
        REQUIRE_THAT(res(3,1,3), WithinRel( df2(x) *   f1(y) *  df3(z), 3e-12));
        REQUIRE_THAT(res(3,2,1), WithinRel( df2(x) *  df1(y) *   f3(z), 3e-12));
        REQUIRE_THAT(res(3,2,2), WithinRel(  f2(x) * d2f1(y) *   f3(z), 3e-12));
        REQUIRE_THAT(res(3,2,3), WithinRel(  f2(x) *  df1(y) *  df3(z), 3e-12));
        REQUIRE_THAT(res(3,3,1), WithinRel( df2(x) *   f1(y) *  df3(z), 3e-12));
        REQUIRE_THAT(res(3,3,2), WithinRel(  f2(x) *  df1(y) *  df3(z), 3e-12));
        REQUIRE_THAT(res(3,3,3), WithinRel(  f2(x) *   f1(y) * d2f3(z), 3e-12));
      }
}


TEST_CASE("TestChebyshevTensorFunc.Value2D")
{
  ChebyshevTensorFunc cheb({cheby_12, cheby_13, cheby_23, cheby_21}, false);

  auto f1 = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto f2 = [](Real x) { return x*x*x - 2.0*x*x + 3.0*x - 1.0; };
  auto f3 = [](Real x) { return 4.0*x*x*x - 7.0*x*x - 3.0*x + 5.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88}) {
      const double p[] {x,y};
      const Vec4 Xc({x, y, 0.0, 0.0}, p);
      const Tensor res = cheb(Xc);
      REQUIRE_THAT(res(1,1), WithinRel(f1(x) * f2(y), 1e-12));
      REQUIRE_THAT(res(2,1), WithinRel(f1(x) * f3(y), 1e-12));
      REQUIRE_THAT(res(1,2), WithinRel(f2(x) * f3(y), 1e-12));
      REQUIRE_THAT(res(2,2), WithinRel(f2(x) * f1(y), 1e-12));
    }
}


TEST_CASE("TestChebyshevTensorFunc.Gradient2D")
{
  ChebyshevTensorFunc cheb({cheby_12, cheby_13, cheby_23, cheby_21}, false);

  auto  f1 = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto df1 = [](Real x) { return 3.0*x*x + 2.0*x + 1.0; };
  auto  f2 = [](Real x) { return x*x*x - 2.0*x*x + 3.0*x - 1.0; };
  auto df2 = [](Real x) { return 3.0*x*x - 4.0*x + 3.0; };
  auto  f3 = [](Real x) { return 4.0*x*x*x - 7.0*x*x - 3.0*x + 5.0; };
  auto df3 = [](Real x) { return 12.0*x*x - 14.0*x - 3.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88}) {
      const double p[] {x,y};
      const Vec4 Xc({x, y, 0.0, 0.0}, p);
      const utl::matrix3d<Real> res = cheb.gradient(Xc);
      REQUIRE_THAT(res(1,1,1), WithinRel(df1(x) *  f2(y), 3e-12));
      REQUIRE_THAT(res(1,1,2), WithinRel( f1(x) * df2(y), 3e-12));
      REQUIRE_THAT(res(2,1,1), WithinRel(df1(x) *  f3(y), 3e-12));
      REQUIRE_THAT(res(2,1,2), WithinRel( f1(x) * df3(y), 3e-12));
      REQUIRE_THAT(res(1,2,1), WithinRel(df2(x) *  f3(y), 3e-12));
      REQUIRE_THAT(res(1,2,2), WithinRel( f2(x) * df3(y), 3e-12));
      REQUIRE_THAT(res(2,2,1), WithinRel(df2(x) *  f1(y), 3e-12));
      REQUIRE_THAT(res(2,2,2), WithinRel( f2(x) * df1(y), 3e-12));
    }
}


TEST_CASE("TestChebyshevTensorFunc.Hessian2D")
{
  ChebyshevTensorFunc cheb({cheby_12, cheby_13, cheby_23, cheby_21}, false);

  auto   f1 = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto  df1 = [](Real x) { return 3.0*x*x + 2.0*x + 1.0; };
  auto d2f1 = [](Real x) { return 6.0*x + 2.0; };
  auto   f2 = [](Real x) { return x*x*x - 2.0*x*x + 3.0*x - 1.0; };
  auto  df2 = [](Real x) { return 3.0*x*x - 4.0*x + 3.0; };
  auto d2f2 = [](Real x) { return 6.0*x - 4.0; };
  auto   f3 = [](Real x) { return 4.0*x*x*x - 7.0*x*x - 3.0*x + 5.0; };
  auto  df3 = [](Real x) { return 12*x*x - 14.0*x - 3.0; };
  auto d2f3 = [](Real x) { return 24.0*x - 14.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88}) {
      const double p[] {x,y};
      const Vec4 Xc({x, y, 0.0, 0.0}, p);
      const utl::matrix4d<Real> res = cheb.hessian(Xc);
      REQUIRE_THAT(res(1,1,1,1), WithinRel(d2f1(x) *   f2(y), 5e-12));
      REQUIRE_THAT(res(1,1,1,2), WithinRel( df1(x) *  df2(y), 5e-12));
      REQUIRE_THAT(res(1,1,2,1), WithinRel( df1(x) *  df2(y), 5e-12));
      REQUIRE_THAT(res(1,1,2,2), WithinRel(  f1(x) * d2f2(y), 5e-12));

      REQUIRE_THAT(res(2,1,1,1), WithinRel(d2f1(x) *   f3(y), 5e-12));
      REQUIRE_THAT(res(2,1,1,2), WithinRel( df1(x) *  df3(y), 5e-12));
      REQUIRE_THAT(res(2,1,2,1), WithinRel( df1(x) *  df3(y), 5e-12));
      REQUIRE_THAT(res(2,1,2,2), WithinRel(  f1(x) * d2f3(y), 5e-12));

      REQUIRE_THAT(res(1,2,1,1), WithinRel(d2f2(x) *   f3(y), 5e-12));
      REQUIRE_THAT(res(1,2,1,2), WithinRel( df2(x) *  df3(y), 5e-12));
      REQUIRE_THAT(res(1,2,2,1), WithinRel( df2(x) *  df3(y), 5e-12));
      REQUIRE_THAT(res(1,2,2,2), WithinRel(  f2(x) * d2f3(y), 5e-12));

      REQUIRE_THAT(res(2,2,1,1), WithinRel(d2f2(x) *   f1(y), 5e-12));
      REQUIRE_THAT(res(2,2,1,2), WithinRel( df2(x) *  df1(y), 5e-12));
      REQUIRE_THAT(res(2,2,2,1), WithinRel( df2(x) *  df1(y), 5e-12));
      REQUIRE_THAT(res(2,2,2,2), WithinRel(  f2(x) * d2f1(y), 5e-12));
    }
}


TEST_CASE("TestChebyshevTensorFunc.Value3D")
{
  ChebyshevTensorFunc cheb({cheby_123, cheby_213, cheby_132,
                            cheby_132, cheby_231, cheby_321,
                            cheby_213, cheby_312, cheby_312}, false);

  auto f1 = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto f2 = [](Real x) { return x*x*x - 2.0*x*x + 3.0*x - 1.0; };
  auto f3 = [](Real x) { return 4.0*x*x*x - 7.0*x*x - 3.0*x + 5.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88})
      for (double z : {0.2, 0.4, 0.6, 0.77}) {
        const double p[] {x,y,z};
        const Vec4 Xc({x, y, z, 0.0}, p);
        const Tensor res = cheb(Xc);
        REQUIRE_THAT(res(1,1), WithinRel(f1(x) * f2(y) * f3(z), 1e-12));
        REQUIRE_THAT(res(2,1), WithinRel(f2(x) * f1(y) * f3(z), 1e-12));
        REQUIRE_THAT(res(3,1), WithinRel(f1(x) * f3(y) * f2(z), 1e-12));

        REQUIRE_THAT(res(1,2), WithinRel(f1(x) * f3(y) * f2(z), 1e-12));
        REQUIRE_THAT(res(2,2), WithinRel(f2(x) * f3(y) * f1(z), 1e-12));
        REQUIRE_THAT(res(3,2), WithinRel(f3(x) * f2(y) * f1(z), 1e-12));

        REQUIRE_THAT(res(1,3), WithinRel(f2(x) * f1(y) * f3(z), 1e-12));
        REQUIRE_THAT(res(2,3), WithinRel(f3(x) * f1(y) * f2(z), 1e-12));
        REQUIRE_THAT(res(3,3), WithinRel(f3(x) * f1(y) * f2(z), 1e-12));
      }
}


TEST_CASE("TestChebyshevTensorFunc.Gradient3D")
{
  ChebyshevTensorFunc cheb({cheby_123, cheby_213, cheby_132,
                            cheby_132, cheby_231, cheby_321,
                            cheby_213, cheby_312, cheby_312}, false);

  auto  f1 = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto df1 = [](Real x) { return 3.0*x*x + 2.0*x + 1.0; };
  auto  f2 = [](Real x) { return x*x*x - 2.0*x*x + 3.0*x - 1.0; };
  auto df2 = [](Real x) { return 3.0*x*x - 4.0*x + 3.0; };
  auto  f3 = [](Real x) { return 4.0*x*x*x - 7.0*x*x - 3.0*x + 5.0; };
  auto df3 = [](Real x) { return 12.0*x*x - 14.0*x - 3.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88})
      for (double z : {0.2, 0.4, 0.6, 0.77}) {
        const double p[] {x,y,z};
        const Vec4 Xc({x, y, z, 0.0}, p);
        const utl::matrix3d<Real> res = cheb.gradient(Xc);
        REQUIRE_THAT(res(1,1,1), WithinRel(df1(x) *  f2(y) *  f3(z), 1e-12));
        REQUIRE_THAT(res(1,1,2), WithinRel( f1(x) * df2(y) *  f3(z), 1e-12));
        REQUIRE_THAT(res(1,1,3), WithinRel( f1(x) *  f2(y) * df3(z), 1e-12));

        REQUIRE_THAT(res(2,1,1), WithinRel(df2(x) *  f1(y) *  f3(z), 1e-12));
        REQUIRE_THAT(res(2,1,2), WithinRel( f2(x) * df1(y) *  f3(z), 1e-12));
        REQUIRE_THAT(res(2,1,3), WithinRel( f2(x) *  f1(y) * df3(z), 1e-12));

        REQUIRE_THAT(res(3,1,1), WithinRel(df1(x) *  f3(y) *  f2(z), 1e-12));
        REQUIRE_THAT(res(3,1,2), WithinRel( f1(x) * df3(y) *  f2(z), 1e-12));
        REQUIRE_THAT(res(3,1,3), WithinRel( f1(x) *  f3(y) * df2(z), 1e-12));

        REQUIRE_THAT(res(1,2,1), WithinRel(df1(x) *  f3(y) *  f2(z), 1e-12));
        REQUIRE_THAT(res(1,2,2), WithinRel( f1(x) * df3(y) *  f2(z), 1e-12));
        REQUIRE_THAT(res(1,2,3), WithinRel( f1(x) *  f3(y) * df2(z), 1e-12));

        REQUIRE_THAT(res(2,2,1), WithinRel(df2(x) *  f3(y) *  f1(z), 1e-12));
        REQUIRE_THAT(res(2,2,2), WithinRel( f2(x) * df3(y) *  f1(z), 1e-12));
        REQUIRE_THAT(res(2,2,3), WithinRel( f2(x) *  f3(y) * df1(z), 1e-12));

        REQUIRE_THAT(res(3,2,1), WithinRel(df3(x) *  f2(y) *  f1(z), 1e-12));
        REQUIRE_THAT(res(3,2,2), WithinRel( f3(x) * df2(y) *  f1(z), 1e-12));
        REQUIRE_THAT(res(3,2,3), WithinRel( f3(x) *  f2(y) * df1(z), 1e-12));

        REQUIRE_THAT(res(1,3,1), WithinRel(df2(x) *  f1(y) *  f3(z), 1e-12));
        REQUIRE_THAT(res(1,3,2), WithinRel( f2(x) * df1(y) *  f3(z), 1e-12));
        REQUIRE_THAT(res(1,3,3), WithinRel( f2(x) *  f1(y) * df3(z), 1e-12));

        REQUIRE_THAT(res(2,3,1), WithinRel(df3(x) *  f1(y) *  f2(z), 1e-12));
        REQUIRE_THAT(res(2,3,2), WithinRel( f3(x) * df1(y) *  f2(z), 1e-12));
        REQUIRE_THAT(res(2,3,3), WithinRel( f3(x) *  f1(y) * df2(z), 1e-12));

        REQUIRE_THAT(res(3,3,1), WithinRel(df3(x) *  f1(y) *  f2(z), 1e-12));
        REQUIRE_THAT(res(3,3,2), WithinRel( f3(x) * df1(y) *  f2(z), 1e-12));
        REQUIRE_THAT(res(3,3,3), WithinRel( f3(x) *  f1(y) * df2(z), 1e-12));
      }
}


TEST_CASE("TestChebyshevTensorFunc.Hessian3D")
{
  ChebyshevTensorFunc cheb({cheby_123, cheby_213, cheby_132,
                            cheby_132, cheby_231, cheby_321,
                            cheby_213, cheby_312, cheby_312}, false);

  auto   f1 = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto  df1 = [](Real x) { return 3.0*x*x + 2.0*x + 1.0; };
  auto d2f1 = [](Real x) { return 6.0*x + 2.0; };
  auto   f2 = [](Real x) { return x*x*x - 2.0*x*x + 3.0*x - 1.0; };
  auto  df2 = [](Real x) { return 3.0*x*x - 4.0*x + 3.0; };
  auto d2f2 = [](Real x) { return 6.0*x - 4.0; };
  auto   f3 = [](Real x) { return 4.0*x*x*x - 7.0*x*x - 3.0*x + 5.0; };
  auto  df3 = [](Real x) { return 12*x*x - 14.0*x - 3.0; };
  auto d2f3 = [](Real x) { return 24.0*x - 14.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88})
      for (double z : {0.2, 0.4, 0.6, 0.77}) {
        const double p[] {x,y,z};
        const Vec4 Xc({x, y, z, 0.0}, p);
        const utl::matrix4d<Real> res = cheb.hessian(Xc);
        REQUIRE_THAT(res(1,1,1,1), WithinRel(d2f1(x) *   f2(y) *   f3(z), 5e-12));
        REQUIRE_THAT(res(1,1,1,2), WithinRel( df1(x) *  df2(y) *   f3(z), 5e-12));
        REQUIRE_THAT(res(1,1,1,3), WithinRel( df1(x) *   f2(y) *  df3(z), 5e-12));
        REQUIRE_THAT(res(1,1,2,1), WithinRel( df1(x) *  df2(y) *   f3(z), 5e-12));
        REQUIRE_THAT(res(1,1,2,2), WithinRel(  f1(x) * d2f2(y) *   f3(z), 5e-12));
        REQUIRE_THAT(res(1,1,2,3), WithinRel(  f1(x) *  df2(y) *  df3(z), 5e-12));
        REQUIRE_THAT(res(1,1,3,1), WithinRel( df1(x) *   f2(y) *  df3(z), 5e-12));
        REQUIRE_THAT(res(1,1,3,2), WithinRel(  f1(x) *  df2(y) *  df3(z), 5e-12));
        REQUIRE_THAT(res(1,1,3,3), WithinRel(  f1(x) *   f2(y) * d2f3(z), 5e-12));

        REQUIRE_THAT(res(2,1,1,1), WithinRel(d2f2(x) *   f1(y) *   f3(z), 5e-12));
        REQUIRE_THAT(res(2,1,1,2), WithinRel( df2(x) *  df1(y) *   f3(z), 5e-12));
        REQUIRE_THAT(res(2,1,1,3), WithinRel( df2(x) *   f1(y) *  df3(z), 5e-12));
        REQUIRE_THAT(res(2,1,2,1), WithinRel( df2(x) *  df1(y) *   f3(z), 5e-12));
        REQUIRE_THAT(res(2,1,2,2), WithinRel(  f2(x) * d2f1(y) *   f3(z), 5e-12));
        REQUIRE_THAT(res(2,1,2,3), WithinRel(  f2(x) *  df1(y) *  df3(z), 5e-12));
        REQUIRE_THAT(res(2,1,3,1), WithinRel( df2(x) *   f1(y) *  df3(z), 5e-12));
        REQUIRE_THAT(res(2,1,3,2), WithinRel(  f2(x) *  df1(y) *  df3(z), 5e-12));
        REQUIRE_THAT(res(2,1,3,3), WithinRel(  f2(x) *   f1(y) * d2f3(z), 5e-12));

        REQUIRE_THAT(res(3,1,1,1), WithinRel(d2f1(x) *   f3(y) *   f2(z), 5e-12));
        REQUIRE_THAT(res(3,1,1,2), WithinRel( df1(x) *  df3(y) *   f2(z), 5e-12));
        REQUIRE_THAT(res(3,1,1,3), WithinRel( df1(x) *   f3(y) *  df2(z), 5e-12));
        REQUIRE_THAT(res(3,1,2,1), WithinRel( df1(x) *  df3(y) *   f2(z), 5e-12));
        REQUIRE_THAT(res(3,1,2,2), WithinRel(  f1(x) * d2f3(y) *   f2(z), 5e-12));
        REQUIRE_THAT(res(3,1,2,3), WithinRel(  f1(x) *  df3(y) *  df2(z), 5e-12));
        REQUIRE_THAT(res(3,1,3,1), WithinRel( df1(x) *   f3(y) *  df2(z), 5e-12));
        REQUIRE_THAT(res(3,1,3,2), WithinRel(  f1(x) *  df3(y) *  df2(z), 5e-12));
        REQUIRE_THAT(res(3,1,3,3), WithinRel(  f1(x) *   f3(y) * d2f2(z), 5e-12));

        REQUIRE_THAT(res(1,2,1,1), WithinRel(d2f1(x) *   f3(y) *   f2(z), 5e-12));
        REQUIRE_THAT(res(1,2,1,2), WithinRel( df1(x) *  df3(y) *   f2(z), 5e-12));
        REQUIRE_THAT(res(1,2,1,3), WithinRel( df1(x) *   f3(y) *  df2(z), 5e-12));
        REQUIRE_THAT(res(1,2,2,1), WithinRel( df1(x) *  df3(y) *   f2(z), 5e-12));
        REQUIRE_THAT(res(1,2,2,2), WithinRel(  f1(x) * d2f3(y) *   f2(z), 5e-12));
        REQUIRE_THAT(res(1,2,2,3), WithinRel(  f1(x) *  df3(y) *  df2(z), 5e-12));
        REQUIRE_THAT(res(1,2,3,1), WithinRel( df1(x) *   f3(y) *  df2(z), 5e-12));
        REQUIRE_THAT(res(1,2,3,2), WithinRel(  f1(x) *  df3(y) *  df2(z), 5e-12));
        REQUIRE_THAT(res(1,2,3,3), WithinRel(  f1(x) *   f3(y) * d2f2(z), 5e-12));

        REQUIRE_THAT(res(2,2,1,1), WithinRel(d2f2(x) *   f3(y) *   f1(z), 5e-12));
        REQUIRE_THAT(res(2,2,1,2), WithinRel( df2(x) *  df3(y) *   f1(z), 5e-12));
        REQUIRE_THAT(res(2,2,1,3), WithinRel( df2(x) *   f3(y) *  df1(z), 5e-12));
        REQUIRE_THAT(res(2,2,2,1), WithinRel( df2(x) *  df3(y) *   f1(z), 5e-12));
        REQUIRE_THAT(res(2,2,2,2), WithinRel(  f2(x) * d2f3(y) *   f1(z), 5e-12));
        REQUIRE_THAT(res(2,2,2,3), WithinRel(  f2(x) *  df3(y) *  df1(z), 5e-12));
        REQUIRE_THAT(res(2,2,3,1), WithinRel( df2(x) *   f3(y) *  df1(z), 5e-12));
        REQUIRE_THAT(res(2,2,3,2), WithinRel(  f2(x) *  df3(y) *  df1(z), 5e-12));
        REQUIRE_THAT(res(2,2,3,3), WithinRel(  f2(x) *   f3(y) * d2f1(z), 5e-12));

        REQUIRE_THAT(res(3,2,1,1), WithinRel(d2f3(x) *   f2(y) *   f1(z), 5e-12));
        REQUIRE_THAT(res(3,2,1,2), WithinRel( df3(x) *  df2(y) *   f1(z), 5e-12));
        REQUIRE_THAT(res(3,2,1,3), WithinRel( df3(x) *   f2(y) *  df1(z), 5e-12));
        REQUIRE_THAT(res(3,2,2,1), WithinRel( df3(x) *  df2(y) *   f1(z), 5e-12));
        REQUIRE_THAT(res(3,2,2,2), WithinRel(  f3(x) * d2f2(y) *   f1(z), 5e-12));
        REQUIRE_THAT(res(3,2,2,3), WithinRel(  f3(x) *  df2(y) *  df1(z), 5e-12));
        REQUIRE_THAT(res(3,2,3,1), WithinRel( df3(x) *   f2(y) *  df1(z), 5e-12));
        REQUIRE_THAT(res(3,2,3,2), WithinRel(  f3(x) *  df2(y) *  df1(z), 5e-12));
        REQUIRE_THAT(res(3,2,3,3), WithinRel(  f3(x) *   f2(y) * d2f1(z), 5e-12));

        REQUIRE_THAT(res(1,3,1,1), WithinRel(d2f2(x) *   f1(y) *   f3(z), 5e-12));
        REQUIRE_THAT(res(1,3,1,2), WithinRel( df2(x) *  df1(y) *   f3(z), 5e-12));
        REQUIRE_THAT(res(1,3,1,3), WithinRel( df2(x) *   f1(y) *  df3(z), 5e-12));
        REQUIRE_THAT(res(1,3,2,1), WithinRel( df2(x) *  df1(y) *   f3(z), 5e-12));
        REQUIRE_THAT(res(1,3,2,2), WithinRel(  f2(x) * d2f1(y) *   f3(z), 5e-12));
        REQUIRE_THAT(res(1,3,2,3), WithinRel(  f2(x) *  df1(y) *  df3(z), 5e-12));
        REQUIRE_THAT(res(1,3,3,1), WithinRel( df2(x) *   f1(y) *  df3(z), 5e-12));
        REQUIRE_THAT(res(1,3,3,2), WithinRel(  f2(x) *  df1(y) *  df3(z), 5e-12));
        REQUIRE_THAT(res(1,3,3,3), WithinRel(  f2(x) *   f1(y) * d2f3(z), 5e-12));

        REQUIRE_THAT(res(2,3,1,1), WithinRel(d2f3(x) *   f1(y) *   f2(z), 5e-12));
        REQUIRE_THAT(res(2,3,1,2), WithinRel( df3(x) *  df1(y) *   f2(z), 5e-12));
        REQUIRE_THAT(res(2,3,1,3), WithinRel( df3(x) *   f1(y) *  df2(z), 5e-12));
        REQUIRE_THAT(res(2,3,2,1), WithinRel( df3(x) *  df1(y) *   f2(z), 5e-12));
        REQUIRE_THAT(res(2,3,2,2), WithinRel(  f3(x) * d2f1(y) *   f2(z), 5e-12));
        REQUIRE_THAT(res(2,3,2,3), WithinRel(  f3(x) *  df1(y) *  df2(z), 5e-12));
        REQUIRE_THAT(res(2,3,3,1), WithinRel( df3(x) *   f1(y) *  df2(z), 5e-12));
        REQUIRE_THAT(res(2,3,3,2), WithinRel(  f3(x) *  df1(y) *  df2(z), 5e-12));
        REQUIRE_THAT(res(2,3,3,3), WithinRel(  f3(x) *   f1(y) * d2f2(z), 5e-12));

        REQUIRE_THAT(res(3,3,1,1), WithinRel(d2f3(x) *   f1(y) *   f2(z), 5e-12));
        REQUIRE_THAT(res(3,3,1,2), WithinRel( df3(x) *  df1(y) *   f2(z), 5e-12));
        REQUIRE_THAT(res(3,3,1,3), WithinRel( df3(x) *   f1(y) *  df2(z), 5e-12));
        REQUIRE_THAT(res(3,3,2,1), WithinRel( df3(x) *  df1(y) *   f2(z), 5e-12));
        REQUIRE_THAT(res(3,3,2,2), WithinRel(  f3(x) * d2f1(y) *   f2(z), 5e-12));
        REQUIRE_THAT(res(3,3,2,3), WithinRel(  f3(x) *  df1(y) *  df2(z), 5e-12));
        REQUIRE_THAT(res(3,3,3,1), WithinRel( df3(x) *   f1(y) *  df2(z), 5e-12));
        REQUIRE_THAT(res(3,3,3,2), WithinRel(  f3(x) *  df1(y) *  df2(z), 5e-12));
        REQUIRE_THAT(res(3,3,3,3), WithinRel(  f3(x) *   f1(y) * d2f2(z), 5e-12));
      }
}
