//==============================================================================
//!
//! \file TestBDF.C
//!
//! \date Oct 10 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for helper functions for BDF based time stepping.
//!
//==============================================================================

#include "BDF.h"

#include "Catch2Support.h"


TEST_CASE("TestBDF.BDF_1")
{
  TimeIntegration::BDF bdf(1);

  bdf.advanceStep();
  REQUIRE(bdf.getOrder() == 1);
  REQUIRE(bdf.getActualOrder() == 1);
  REQUIRE_THAT(bdf[0], WithinRel(1.0));
  REQUIRE_THAT(bdf[1], WithinRel(-1.0));
  bdf.advanceStep();
  REQUIRE(bdf.getOrder() == 1);
  REQUIRE_THAT(bdf[0], WithinRel(1.0));
  REQUIRE_THAT(bdf[1], WithinRel(-1.0));
}


TEST_CASE("TestBDF.BDF_2")
{
  TimeIntegration::BDF bdf(2);

  bdf.advanceStep();
  REQUIRE(bdf.getOrder() == 1);
  REQUIRE(bdf.getActualOrder() == 2);
  REQUIRE_THAT(bdf[0], WithinRel(1.0));
  REQUIRE_THAT(bdf[1], WithinRel(-1.0));
  bdf.advanceStep();
  REQUIRE(bdf.getOrder() == 2);
  REQUIRE_THAT(bdf[0], WithinRel(1.5));
  REQUIRE_THAT(bdf[1], WithinRel(-2.0));
  REQUIRE_THAT(bdf[2], WithinRel(0.5));
}


TEST_CASE("TestBDF.BDFD2_1")
{
  TimeIntegration::BDFD2 bdf(1);

  bdf.advanceStep(0.1, 0.1);
  REQUIRE(bdf.getOrder() == 1);
  REQUIRE(bdf.getActualOrder() == 1);
  REQUIRE_THAT(bdf[0], WithinRel(2.0));
  REQUIRE_THAT(bdf[1], WithinRel(-2.0));
  bdf.advanceStep(0.1, 0.1);
  REQUIRE(bdf.getOrder() == 1);
  REQUIRE_THAT(bdf[0], WithinRel(1.0));
  REQUIRE_THAT(bdf[1], WithinRel(-2.0));
  REQUIRE_THAT(bdf[2], WithinRel(1.0));
}


TEST_CASE("TestBDF.BDFD2_2")
{
  TimeIntegration::BDFD2 bdf(2);

  bdf.advanceStep(0.1, 0.1);
  REQUIRE(bdf.getOrder() == 1);
  REQUIRE(bdf.getActualOrder() == 2);
  REQUIRE_THAT(bdf[0], WithinRel(2.0));
  REQUIRE_THAT(bdf[1], WithinRel(-2.0));
  bdf.advanceStep(0.1, 0.1);
  REQUIRE(bdf.getOrder() == 1);
  REQUIRE_THAT(bdf[0], WithinRel(2.5));
  REQUIRE_THAT(bdf[1], WithinRel(-8.0));
  REQUIRE_THAT(bdf[2], WithinRel(5.5));
  bdf.advanceStep(0.1, 0.1);
  REQUIRE(bdf.getOrder() == 2);
  REQUIRE_THAT(bdf[0], WithinRel(2.0));
  REQUIRE_THAT(bdf[1], WithinRel(-5.0));
  REQUIRE_THAT(bdf[2], WithinRel(4.0));
  REQUIRE_THAT(bdf[3], WithinRel(-1.0));
}


namespace {

//! \brief Applies coefficients to the values of a function at given times.
template<class Function>
double apply (const TimeIntegration::BDF& bdf, const std::vector<double>& t,
              const Function& f)
{
  double sum = 0.0;
  for (size_t k = 0; k < t.size(); k++)
    sum += bdf[k]*f(t[k]);
  return sum;
}

}


TEST_CASE("TestBDF.BDF_2_Variable")
{
  const double tau = GENERATE(0.5, 2.0, 5.0);
  const double dt1 = 0.1, dt2 = tau*dt1, dt3 = 0.3*dt2;

  TimeIntegration::BDF bdf(2);
  bdf.advanceStep(dt1, 0.0);
  bdf.advanceStep(dt2, dt1);
  REQUIRE(bdf.getOrder() == 2);

  // Exact for the derivatives of polynomials of second order
  std::vector<double> t = { dt1 + dt2, dt1, 0.0 };
  REQUIRE_THAT(apply(bdf, t, [](double x) { return x; }) / dt2,
               WithinRel(1.0, 1e-12));
  REQUIRE_THAT(apply(bdf, t, [](double x) { return x*x; }) / dt2,
               WithinRel(2.0*t[0], 1e-12));

  // Extrapolation is exact for a linear function
  auto lin = [](double x) { return 1.0 + 3.0*x; };
  REQUIRE_THAT(bdf.extrapolate(std::vector<double>{lin(t[1]), lin(t[2])}),
               WithinRel(lin(t[0]), 1e-12));

  // And again on the next step with yet another step size
  bdf.advanceStep(dt3, dt2);
  t = { dt1 + dt2 + dt3, dt1 + dt2, dt1 };
  REQUIRE_THAT(apply(bdf, t, [](double x) { return x*x; }) / dt3,
               WithinRel(2.0*t[0], 1e-12));
  REQUIRE_THAT(bdf.extrapolate(std::vector<double>{lin(t[1]), lin(t[2])}),
               WithinRel(lin(t[0]), 1e-12));
}


TEST_CASE("TestBDF.BDFD2_1_Variable")
{
  const double tau = GENERATE(0.5, 2.0, 5.0);
  const double dt1 = 0.1, dt2 = tau*dt1;

  TimeIntegration::BDFD2 bdf(1);
  bdf.advanceStep(dt1, 0.0);
  bdf.advanceStep(dt2, dt1);

  // Exact for the second derivative of a polynomial of second order
  const std::vector<double> t = { dt1 + dt2, dt1, 0.0 };
  REQUIRE_THAT(apply(bdf, t, [](double x) { return x*x; }) / (dt2*dt2),
               WithinRel(2.0, 1e-12));
}


TEST_CASE("TestBDF.BDFD2_2_Variable")
{
  const double tau = GENERATE(0.5, 2.0, 5.0);
  const double dt1 = 0.1, dt2 = tau*dt1, dt3 = 0.7*dt2;

  TimeIntegration::BDFD2 bdf(2);
  bdf.advanceStep(dt1, 0.0);
  bdf.advanceStep(dt2, dt1);
  bdf.advanceStep(dt3, dt2);
  REQUIRE(bdf.getOrder() == 2);

  // Exact for the second derivatives of polynomials of third order
  const std::vector<double> t = { dt1 + dt2 + dt3, dt1 + dt2, dt1, 0.0 };
  REQUIRE_THAT(apply(bdf, t, [](double x) { return x*x; }) / (dt3*dt3),
               WithinRel(2.0, 1e-12));
  REQUIRE_THAT(apply(bdf, t, [](double x) { return x*x*x; }) / (dt3*dt3),
               WithinRel(6.0*t[0], 1e-12));
}
