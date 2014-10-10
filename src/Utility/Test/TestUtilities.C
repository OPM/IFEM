//==============================================================================
//!
//! \file TestUtilities.C
//!
//! \date Oct 10 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for various utility methods.
//!
//==============================================================================

#include "Utilities.h"
#include "tinyxml.h"

#include "gtest/gtest.h"

TEST(TestUtilities, ParseIntegers)
{
  const char* simple = "1";
  const char* ranged = "1:5";

  std::vector<int> values, values2;
  utl::parseIntegers(values, simple);
  utl::parseIntegers(values2, ranged);

  ASSERT_EQ(values.size(), 1U);
  ASSERT_EQ(values2.size(), 5U);
  ASSERT_EQ(values[0], 1);
  for (int i=0;i<5;++i)
    ASSERT_EQ(values2[i], i+1);
}

TEST(TestUtilities, ParseKnots)
{
  const char* simple = "0 0 0.5 1.0 1.0";
  char* meh = strdup(simple);
  strtok(meh, " ");
  std::vector<Real> xi;
  utl::parseKnots(xi);
  ASSERT_EQ(xi.size(), 4U);
  ASSERT_FLOAT_EQ(xi[0], 0.0);
  ASSERT_FLOAT_EQ(xi[1], 0.5);
  ASSERT_FLOAT_EQ(xi[2], 1.0);
  ASSERT_FLOAT_EQ(xi[3], 1.0);
}

TEST(TestUtilities, ReadLine)
{
  std::stringstream str;
  str << "Blah blah\n"
      << "# Commented line\n"
      << "Blah bluh\n";
  char* tmp = utl::readLine(str);
  ASSERT_STREQ(tmp, "Blah blah");
  tmp = utl::readLine(str);
  ASSERT_STREQ(tmp, "Blah bluh");
}

TEST(TestUtilities, IgnoreComments)
{
  std::stringstream str;
  str << "Blah blah\n"
      << "# Commented line\n"
      << "Blah bluh\n";
  char tmp[1024];
  str.getline(tmp,1024);
  ASSERT_STREQ(tmp, "Blah blah");
  utl::ignoreComments(str);
  str.getline(tmp,1024);
  ASSERT_STREQ(tmp, "Blah bluh");
}

TEST(TestUtilities, GetAttribute)
{
  TiXmlDocument doc;
  doc.LoadFile("src/Utility/Test/refdata/getattribute.xml");

  if (!doc.RootElement())
    ASSERT_TRUE(false);

  const std::vector<std::string> truebool = {"booltrue", "boolone",
                                             "boolon", "boolyes"};
  const std::vector<std::string> falsebool = {"boolfalse", "boolzero",
                                              "booloff", "boolno"};

  for (size_t i=0;i<truebool.size();++i) {
    const TiXmlElement* elem = doc.RootElement()->FirstChildElement(truebool[i]);
    bool b=false;
    ASSERT_TRUE(elem != NULL);
    ASSERT_TRUE(utl::getAttribute(elem, "bar", b));
    ASSERT_TRUE(b);
  }
  for (size_t i=0;i<falsebool.size();++i) {
    const TiXmlElement* elem = doc.RootElement()->FirstChildElement(falsebool[i]);
    bool b=false;
    ASSERT_TRUE(utl::getAttribute(elem, "bar", b));
    ASSERT_FALSE(b);
  }
  const TiXmlElement* elem = doc.RootElement()->FirstChildElement("intval");
  int val;
  ASSERT_TRUE(elem != NULL);
  ASSERT_TRUE(utl::getAttribute(elem, "bar", val));
  ASSERT_EQ(val, 1);
  size_t val2;
  ASSERT_TRUE(utl::getAttribute(elem, "bar", val2));
  ASSERT_EQ(val2, 1U);
  Real val3;
  elem = doc.RootElement()->FirstChildElement("realval");
  ASSERT_TRUE(elem != NULL);
  ASSERT_TRUE(utl::getAttribute(elem, "bar", val3));
  ASSERT_FLOAT_EQ(val3, 2.01);
}
