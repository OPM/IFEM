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

#include "gtest/gtest.h"

TEST(TestStringUtils, ReplaceAll)
{
  std::string test1("abababab");
  replaceAll(test1, "ab", "ba");

  EXPECT_STREQ(test1.c_str(), "babababa");
}
