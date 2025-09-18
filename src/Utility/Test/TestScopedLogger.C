//==============================================================================
//!
//! \file TestScopedLogger.C
//!
//! \date Oct 10 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for scoped logger.
//!
//==============================================================================

#include "ScopedLogger.h"

#include <catch2/catch_test_macros.hpp>

#include <sstream>


TEST_CASE("TestScopedLogger.General")
{
  std::stringstream str;
  auto&& MockFunction = [&str]() { ScopedLogger log("MockFunction",str); };

  MockFunction();
  char tmp[1024];
  str.getline(tmp, 1024);
  std::cout << tmp << std::endl;
#ifdef HAVE_MPI
  REQUIRE(!strcmp(tmp, "[0]: Entering \"MockFunction\""));
#else
  REQUIRE(!strcmp(tmp, "Entering \"MockFunction\""));
#endif
  str.getline(tmp, 1024);
  std::cout << tmp << std::endl;
#ifdef HAVE_MPI
  REQUIRE(!strcmp(tmp, "[0]: Exiting \"MockFunction\""));
#else
  REQUIRE(!strcmp(tmp, "Exiting \"MockFunction\""));
#endif
}
