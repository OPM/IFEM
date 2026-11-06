//==============================================================================
//!
//! \file TestAndersonAcceleration.C
//!
//! \date Sep 18 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for Anderson acceleration of a fixed point iteration.
//!
//==============================================================================

#include "AndersonAcceleration.h"

#include "Catch2Support.h"

#include <cmath>
#include <functional>

namespace {

//! \brief Runs a fixed point iteration to a tolerance on its residual.
//! \param[in] G The fixed point map
//! \param x Starting point, holding the point reached on return
//! \param acc The acceleration to apply, of depth zero for none
//! \param[in] tol Tolerance on the norm of the residual
//! \param[in] maxIt Iteration limit
//! \return Number of iterations taken, or \a maxIt if the tolerance was
//! not reached
int iterate (const std::function<void(const Vector&,Vector&)>& G,
             Vector& x, AndersonAcceleration& acc,
             double tol = 1.0e-10, int maxIt = 10000)
{
  Vector gx(x.size());
  for (int it = 1; it <= maxIt; ++it) {
    G(x,gx);
    Vector f(gx);
    f -= x;
    if (f.norm2() < tol) {
      x = gx;
      return it;
    }
    acc.accelerate(x,gx);
    x = gx;
  }

  return maxIt;
}


//! \brief A linear fixed point map with a spread spectrum.
//! \details The map is x <- x + w*(b - A*x) with A diagonal, a relaxed
//! Richardson iteration whose convergence factors are 1-w*a_i. Its spectrum
//! is spread, which is what makes a two point extrapolation of little use.
void richardson (const Vector& x, Vector& gx)
{
  const size_t n = x.size();
  const double w = 0.1;
  for (size_t i = 1; i <= n; ++i) {
    const double a = 1.0 + 9.0*static_cast<double>(i-1)/(n-1); // 1 .. 10
    const double b = 1.0;
    gx(i) = x(i) + w*(b - a*x(i));
  }
}

}


TEST_CASE("TestAndersonAcceleration.DisabledIsPlainIteration")
{
  AndersonAcceleration none(0);
  Vector x(8), gx(8);
  richardson(x,gx);
  Vector expected(gx);

  // Accelerating with no depth must leave the step exactly as it was
  REQUIRE(none.accelerate(x,gx));
  REQUIRE(none.getSize() == 0);
  for (size_t i = 1; i <= gx.size(); ++i)
    REQUIRE_THAT(gx(i), WithinRel(expected(i), 1.0e-14));
}


TEST_CASE("TestAndersonAcceleration.ReachesTheSameFixedPoint")
{
  // The fixed point of the map above is x_i = 1/a_i
  const size_t n = 16;
  for (size_t depth : {size_t(0), size_t(2), size_t(5)}) {
    AndersonAcceleration acc(depth);
    Vector x(n);
    const int its = iterate(richardson,x,acc);
    REQUIRE(its < 10000);
    for (size_t i = 1; i <= n; ++i) {
      const double a = 1.0 + 9.0*static_cast<double>(i-1)/(n-1);
      REQUIRE_THAT(x(i), WithinAbs(1.0/a, 1.0e-8));
    }
  }
}


TEST_CASE("TestAndersonAcceleration.ConvergesInFewerIterations")
{
  const size_t n = 16;
  Vector x0(n), x2(n), x5(n);
  AndersonAcceleration none(0), two(2), five(5);

  const int plain = iterate(richardson,x0,none);
  const int acc2 = iterate(richardson,x2,two);
  const int acc5 = iterate(richardson,x5,five);

  // The plain iteration needs a couple of hundred steps on this map, and the
  // acceleration should cut that severalfold. A factor of two is asserted,
  // the measured factor being four, so that the test does not depend on
  // exactly how the least squares problem resolves
  REQUIRE(plain > 100);
  REQUIRE(acc2 < plain);
  REQUIRE(2*acc5 < plain);
}


TEST_CASE("TestAndersonAcceleration.HistoryIsBoundedByTheDepth")
{
  AndersonAcceleration acc(3);
  Vector x(8), gx(8);
  for (int it = 0; it < 10; ++it) {
    richardson(x,gx);
    REQUIRE(acc.accelerate(x,gx));
    REQUIRE(acc.getSize() <= 3);
    x = gx;
  }
  REQUIRE(acc.getSize() == 3);
}


TEST_CASE("TestAndersonAcceleration.DropsTheHistoryOnBlowUp")
{
  // A growth limit of one drops the history whenever the residual rises,
  // which the map below does on its second step
  AndersonAcceleration acc(5,1.0e-10,1.0);
  Vector x(4), gx(4);

  richardson(x,gx);
  REQUIRE(acc.accelerate(x,gx));
  x = gx;
  richardson(x,gx);
  REQUIRE(acc.accelerate(x,gx));
  REQUIRE(acc.getSize() == 1);

  // A step that lands far away raises the residual, so the history goes
  Vector far(4);
  for (size_t i = 1; i <= 4; ++i)
    far(i) = 1.0e3;
  REQUIRE(acc.accelerate(x,far));
  REQUIRE(acc.getSize() == 0);
}


TEST_CASE("TestAndersonAcceleration.ClearForgetsTheHistory")
{
  AndersonAcceleration acc(4);
  Vector x(8), gx(8);
  for (int it = 0; it < 4; ++it) {
    richardson(x,gx);
    acc.accelerate(x,gx);
    x = gx;
  }
  REQUIRE(acc.getSize() > 0);
  acc.clear();
  REQUIRE(acc.getSize() == 0);
}
