//==============================================================================
//!
//! \file TestStringUtils.C
//!
//! \date Oct 7 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for some general string manipulation functions.
//!
//==============================================================================

#include "StringUtils.h"

#include <catch2/catch_test_macros.hpp>

TEST_CASE("TestStringUtils.ReplaceAll")
{
  std::string test1("abababab");
  replaceAll(test1, "ab", "ba");

  REQUIRE(test1 == "babababa");
}


TEST_CASE("TestStringUtils.SplitString")
{
  std::string test1("ab ba ba\tab");
  std::vector<std::string> s = splitString(test1);

  REQUIRE(s.size() == 4);

  REQUIRE(s[0] == "ab");
  REQUIRE(s[1] == "ba");
  REQUIRE(s[2] == "ba");
  REQUIRE(s[3] == "ab");
}
