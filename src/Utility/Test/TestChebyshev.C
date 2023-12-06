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

#include "gtest/gtest.h"


namespace {

#include "refdata/cheby_strings.h"

}


TEST(TestChebyshevFunc, Value1D)
{
  ChebyshevFunc cheb(cheby_1, false);

  auto f = [](Real x) { return x*x*x + x*x + x - 1.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0}) {
    const Vec4 Xc({x,0.0, 0.0, 0.0}, &x);
    EXPECT_NEAR(cheb(Xc), f(x), 1e-12);
  }
}


TEST(TestChebyshevFunc, Value2D)
{
  ChebyshevFunc cheb(cheby_2, false);

  auto f = [](Real x) { return x*x*x + x*x + x - 1.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88}) {
      const double p[] {x,y};
      const Vec4 Xc({x, y, 0.0, 0.0}, p);
      EXPECT_NEAR(cheb(Xc), f(x) * f(y), 1e-12);
    }
}


TEST(TestChebyshevFunc, Value3D)
{
  ChebyshevFunc cheb(cheby_3, false);

  auto f = [](Real x) { return x*x*x + x*x + x - 1.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88})
      for (double z : {0.2, 0.4, 0.6, 0.77}) {
        const double p[] {x,y,z};
        const Vec4 Xc({x, y, z, 0.0}, p);
        EXPECT_NEAR(cheb(Xc), f(x) * f(y) * f(z), 1e-12);
      }
}


TEST(TestChebyshevVecFunc, Value2D)
{
  ChebyshevVecFunc cheb({cheby_12, cheby_21}, false, false);

  auto f1 = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto f2 = [](Real x) { return x*x*x - 2.0*x*x + 3.0*x - 1.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88}) {
      const double p[] {x,y};
      const Vec4 Xc({x, y, 0.0, 0.0}, p);
      const Vec3 res = cheb(Xc);
      EXPECT_NEAR(res[0], f1(x) * f2(y), 1e-12);
      EXPECT_NEAR(res[1], f1(y) * f2(x), 1e-12);
      EXPECT_DOUBLE_EQ(res[2], 0.0);
    }
}


TEST(TestChebyshevVecFunc, Value3D)
{
  ChebyshevVecFunc cheb({cheby_123, cheby_132, cheby_213}, false, false);

  auto f1 = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto f2 = [](Real x) { return x*x*x - 2.0*x*x + 3.0*x - 1.0; };
  auto f3 = [](Real x) { return 4.0*x*x*x - 7.0*x*x - 3.0*x + 5.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88})
      for (double z : {0.2, 0.4, 0.6, 0.77}) {
        const double p[] {x,y,z};
        const Vec4 Xc({x, y, z, 0.0}, p);
        const Vec3 res = cheb(Xc);
        EXPECT_NEAR(res[0], f1(x) * f2(y) * f3(z), 1e-12);
        EXPECT_NEAR(res[1], f1(x) * f2(z) * f3(y), 1e-12);
        EXPECT_NEAR(res[2], f2(x) * f1(y) * f3(z), 1e-12);
      }
}


TEST(TestChebyshevTensorFunc, Value2D)
{
  ChebyshevTensorFunc cheb({cheby_12, cheby_23, cheby_13, cheby_21}, false, false);

  auto f1 = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto f2 = [](Real x) { return x*x*x - 2.0*x*x + 3.0*x - 1.0; };
  auto f3 = [](Real x) { return 4.0*x*x*x - 7.0*x*x - 3.0*x + 5.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88}) {
      const double p[] {x,y};
      const Vec4 Xc({x, y, 0.0, 0.0}, p);
      const Tensor res = cheb(Xc);
      EXPECT_NEAR(res(1,1), f1(x) * f2(y), 1e-12);
      EXPECT_NEAR(res(2,1), f1(x) * f3(y), 1e-12);
      EXPECT_NEAR(res(1,2), f2(x) * f3(y), 1e-12);
      EXPECT_NEAR(res(2,2), f2(x) * f1(y), 1e-12);
    }
}


TEST(TestChebyshevTensorFunc, Value3D)
{
  ChebyshevTensorFunc cheb({cheby_123, cheby_132, cheby_213,
                            cheby_132, cheby_231, cheby_312,
                            cheby_132, cheby_321, cheby_312}, false, false);

  auto f1 = [](Real x) { return x*x*x + x*x + x - 1.0; };
  auto f2 = [](Real x) { return x*x*x - 2.0*x*x + 3.0*x - 1.0; };
  auto f3 = [](Real x) { return 4.0*x*x*x - 7.0*x*x - 3.0*x + 5.0; };

  for (double x : {0.0, 0.2, 0.7, 0.5, 0.66, 1.0})
    for (double y : {0.1, 0.5, 0.7, 0.88})
      for (double z : {0.2, 0.4, 0.6, 0.77}) {
        const double p[] {x,y,z};
        const Vec4 Xc({x, y, z, 0.0}, p);
        const Tensor res = cheb(Xc);
        EXPECT_NEAR(res(1,1), f1(x) * f2(y) * f3(z), 1e-12);
        EXPECT_NEAR(res(1,2), f1(x) * f3(y) * f2(z), 1e-12);
        EXPECT_NEAR(res(1,3), f2(x) * f1(y) * f3(z), 1e-12);
        EXPECT_NEAR(res(2,1), f1(x) * f3(y) * f2(z), 1e-12);
        EXPECT_NEAR(res(2,2), f2(x) * f3(y) * f1(z), 1e-12);
        EXPECT_NEAR(res(2,3), f3(x) * f1(y) * f2(z), 1e-12);
        EXPECT_NEAR(res(3,1), f1(x) * f3(y) * f2(z), 1e-12);
        EXPECT_NEAR(res(3,2), f3(x) * f2(y) * f1(z), 1e-12);
        EXPECT_NEAR(res(3,3), f3(x) * f1(y) * f2(z), 1e-12);
      }
}
