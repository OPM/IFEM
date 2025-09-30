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
