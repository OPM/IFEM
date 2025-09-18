//==============================================================================
//!
//! \file TestTimeIntUtils.C
//!
//! \date Oct 7 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for various helpers for time integration.
//!
//==============================================================================

#include "TimeIntUtils.h"

#include <catch2/catch_test_macros.hpp>

using namespace TimeIntegration;


TEST_CASE("TestTimeIntUtils.Order")
{
  REQUIRE(Order(EULER) == 1);
  REQUIRE(Order(BE) ==    1);
  REQUIRE(Order(HEUN) ==  2);
  REQUIRE(Order(BDF2) ==  2);
  REQUIRE(Order(RK3) ==   3);
  REQUIRE(Order(RK4) ==   4);
}


TEST_CASE("TestTimeIntUtils.Steps")
{
  REQUIRE(Steps(EULER) == 1);
  REQUIRE(Steps(BE) ==    1);
  REQUIRE(Steps(HEUN) ==  1);
  REQUIRE(Steps(BDF2) ==  2);
  REQUIRE(Steps(RK3) ==   1);
  REQUIRE(Steps(RK4) ==   1);
}
