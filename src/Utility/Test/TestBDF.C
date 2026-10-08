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
#include "Tensor.h"

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

  // And over a tensor, component by component
  auto linT = [&lin](double x) {
    Tensor T(2);
    T(1,1) = lin(x);  T(1,2) = 2.0*lin(x);
    T(2,1) = -lin(x); T(2,2) = 0.5*lin(x);
    return T;
  };
  const Tensor Tex = bdf.extrapolate(std::vector<Tensor>{linT(t[1]), linT(t[2])});
  const Tensor Tref = linT(t[0]);
  for (int i = 1; i <= 2; ++i)
    for (int j = 1; j <= 2; ++j)
      REQUIRE_THAT(Tex(i,j), WithinRel(Tref(i,j), 1e-12));

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


#ifdef HAS_CEREAL
TEST_CASE("TestBDF.Serialize")
{
  TimeIntegration::BDF bdf(2);
  bdf.advanceStep(0.1,0.1);
  bdf.advanceStep(0.1,0.1);
  bdf.advanceStep(0.2,0.1);

  std::map<std::string,std::string> data;
  REQUIRE(bdf.serialize(data,"Test"));

  // The restored scheme continues where the serialized one was
  TimeIntegration::BDF restored(2);
  REQUIRE(restored.deSerialize(data,"Test"));
  REQUIRE(restored.getOrder() == 2);
  REQUIRE_THAT(restored.getStepRatio(), WithinRel(bdf.getStepRatio()));
  for (int i = 0; i < 3; i++)
    REQUIRE_THAT(restored[i], WithinRel(bdf[i]));

  bdf.advanceStep(0.1,0.2);
  restored.advanceStep(0.1,0.2);
  for (int i = 0; i < 3; i++)
    REQUIRE_THAT(restored[i], WithinRel(bdf[i]));

  // Without a serialized state, the scheme is advanced to its full order
  TimeIntegration::BDF fresh(2);
  REQUIRE(fresh.deSerialize(data,"Other"));
  fresh.advanceStep();
  REQUIRE(fresh.getOrder() == 2);

  // A stationary scheme has no state, and is restored without one
  TimeIntegration::BDF stationary;
  std::map<std::string,std::string> none;
  REQUIRE(stationary.serialize(none,"Test"));
  REQUIRE(none.empty());
  REQUIRE(stationary.deSerialize(none,"Test"));
}


TEST_CASE("TestBDF.SerializeBDFD2")
{
  // The step size before the previous one is part of the state
  TimeIntegration::BDFD2 bdf(2);
  bdf.advanceStep(0.1,0.1);
  bdf.advanceStep(0.2,0.1);
  bdf.advanceStep(0.1,0.2);

  std::map<std::string,std::string> data;
  REQUIRE(bdf.serialize(data,"Test"));

  TimeIntegration::BDFD2 restored(2);
  REQUIRE(restored.deSerialize(data,"Test"));
  bdf.advanceStep(0.3,0.1);
  restored.advanceStep(0.3,0.1);
  for (int i = 0; i < 4; i++)
    REQUIRE_THAT(restored[i], WithinRel(bdf[i]));

  // Without a serialized state, the scheme continues at its full order,
  // past the start-up coefficients of its second step
  TimeIntegration::BDFD2 fresh(2);
  REQUIRE(fresh.deSerialize(data,"Other"));
  fresh.advanceStep();
  REQUIRE(fresh.getOrder() == 2);
  TimeIntegration::BDFD2 full(2);
  for (int i = 0; i < 3; i++)
    full.advanceStep();
  REQUIRE(fresh.getCoefs() == full.getCoefs());
}
#endif
