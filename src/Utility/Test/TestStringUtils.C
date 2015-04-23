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

TEST(TestStringUtils, SplitString)
{
  std::string test1("ab ba ba\tab");
  std::vector<std::string> s = splitString(test1);

  EXPECT_EQ(s.size(), 4u);

  EXPECT_STREQ(s[0].c_str(), "ab");
  EXPECT_STREQ(s[1].c_str(), "ba");
  EXPECT_STREQ(s[2].c_str(), "ba");
  EXPECT_STREQ(s[3].c_str(), "ab");
}
