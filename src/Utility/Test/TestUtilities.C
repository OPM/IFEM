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
#include "Vec3.h"
#include "tinyxml2.h"

#include "gtest/gtest.h"


TEST(TestUtilities, ParseIntegers)
{
  const char* simple = "1";
  const char* ranged = "1:5";
  const char* multip = "1 3 5:8 10";
  std::vector<int> values1, values2, values3;
  utl::parseIntegers(values1, simple);
  utl::parseIntegers(values2, ranged);
  utl::parseIntegers(values3, multip);

  int i;
  ASSERT_EQ(values1.size(), 1U);
  ASSERT_EQ(values2.size(), 5U);
  ASSERT_EQ(values3.size(), 7U);
  EXPECT_EQ(values1.front(), 1);
  for (i = 0; i < 5; i++)
    EXPECT_EQ(values2[i], 1+i);
  for (i = 0; i < 3; i++)
    EXPECT_EQ(values3[i], 1+2*i);
  for (i = 2; i < 6; i++)
    EXPECT_EQ(values3[i], 3+i);
  EXPECT_EQ(values3[6], 10);
}

TEST(TestUtilities, ParseKnots)
{
  std::string simple("xi 0 0.5 0.9 1.0");
  strtok(const_cast<char*>(simple.c_str())," ");
  std::vector<double> xi;
  utl::parseKnots(xi);
  ASSERT_EQ(xi.size(), 4U);
  EXPECT_FLOAT_EQ(xi[0], 0.0);
  EXPECT_FLOAT_EQ(xi[1], 0.5);
  EXPECT_FLOAT_EQ(xi[2], 0.9);
  EXPECT_FLOAT_EQ(xi[3], 1.0);

  std::string graded5("xi C5 79 0.01 2.0");
  strtok(const_cast<char*>(graded5.c_str())," ");
  xi.clear();
  utl::parseKnots(xi);
  ASSERT_EQ(xi.size(), 79U);
  EXPECT_FLOAT_EQ(xi[39], 0.5);
  EXPECT_FLOAT_EQ(xi.front()+xi.back(), 1.0);

  for (size_t i = 0; i < xi.size(); i++)
    std::cout <<"xi["<< i <<"] = "<< xi[i] << std::endl;
}

TEST(TestUtilities, ReadLine)
{
  std::stringstream str;
  str << "Blah blah\n"
      << "# Commented line\n"
      << "Blah bluh\n";
  EXPECT_STREQ(utl::readLine(str), "Blah blah");
  EXPECT_STREQ(utl::readLine(str), "Blah bluh");
}

TEST(TestUtilities, IgnoreComments)
{
  std::stringstream str;
  str << "Blah blah\n"
      << "# Commented line\n"
      << "Blah bluh\n";
  char tmp[1024];
  str.getline(tmp,1024);
  EXPECT_STREQ(tmp, "Blah blah");
  utl::ignoreComments(str);
  str.getline(tmp,1024);
  EXPECT_STREQ(tmp, "Blah bluh");
}

TEST(TestUtilities, GetAttribute)
{
  const char* input = "<xml>"
    "  <booltrue bar=\"true\"/>"
    "  <boolone bar=\"1\"/>"
    "  <boolon bar=\"on\"/>"
    "  <boolyes bar=\"yes\"/>"
    "  <boolfalse bar=\"false\"/>"
    "  <boolzero bar=\"0\"/>"
    "  <booloff bar=\"off\"/>"
    "  <boolno bar=\"no\"/>"
    "  <intval bar=\"1\"/>"
    "  <realval bar=\"2.01\"/>"
    "  <vecval1 bar=\"1.2\"/>"
    "  <vecval2 bar=\"1.2 3.4\"/>"
    "  <vecval3 bar=\"1.2 3.4 5.6\"/>"
    "</xml>";

  tinyxml2::XMLDocument doc;
  doc.Parse(input);
  ASSERT_TRUE(doc.RootElement() != nullptr);

  const char* bTrue[4]  = { "booltrue" , "boolone" , "boolon" , "boolyes" };
  const char* bFalse[4] = { "boolfalse", "boolzero", "booloff", "boolno"  };

  int i;
  for (i = 0; i < 4; i++) {
    const tinyxml2::XMLElement* elem = doc.RootElement()->FirstChildElement(bTrue[i]);
    ASSERT_TRUE(elem != nullptr);
    bool b = false;
    EXPECT_TRUE(utl::getAttribute(elem, "bar", b));
    EXPECT_TRUE(b);
  }

  for (i = 0; i < 4; i++) {
    const tinyxml2::XMLElement* elem = doc.RootElement()->FirstChildElement(bFalse[i]);
    ASSERT_TRUE(elem != nullptr);
    bool b = false;
    EXPECT_TRUE(utl::getAttribute(elem, "bar", b));
    EXPECT_FALSE(b);
  }

  const tinyxml2::XMLElement* elem = doc.RootElement()->FirstChildElement("intval");
  ASSERT_TRUE(elem != nullptr);
  int val1 = 0;
  EXPECT_TRUE(utl::getAttribute(elem, "bar", val1));
  EXPECT_EQ(val1, 1);
  size_t val2 = 0;
  EXPECT_TRUE(utl::getAttribute(elem, "bar", val2));
  EXPECT_EQ(val2, 1U);

  elem = doc.RootElement()->FirstChildElement("realval");
  ASSERT_TRUE(elem != nullptr);
  double val3 = 0.0;
  EXPECT_TRUE(utl::getAttribute(elem, "bar", val3));
  EXPECT_FLOAT_EQ(val3, 2.01);

  elem = doc.RootElement()->FirstChildElement("vecval1");
  ASSERT_TRUE(elem != nullptr);
  Vec3 val4;
  EXPECT_TRUE(utl::getAttribute(elem, "bar", val4));
  EXPECT_FLOAT_EQ(val4.x, 1.2);
  EXPECT_FLOAT_EQ(val4.y, 0.0);
  EXPECT_FLOAT_EQ(val4.z, 0.0);
  EXPECT_TRUE(utl::getAttribute(elem, "bar", val4, 2));
  EXPECT_FLOAT_EQ(val4.x, 1.2);
  EXPECT_FLOAT_EQ(val4.y, 1.2);
  EXPECT_FLOAT_EQ(val4.z, 0.0);
  EXPECT_TRUE(utl::getAttribute(elem, "bar", val4, 3));
  EXPECT_FLOAT_EQ(val4.x, 1.2);
  EXPECT_FLOAT_EQ(val4.y, 1.2);
  EXPECT_FLOAT_EQ(val4.z, 1.2);

  elem = doc.RootElement()->FirstChildElement("vecval2");
  ASSERT_TRUE(elem != nullptr);
  EXPECT_TRUE(utl::getAttribute(elem, "bar", val4));
  EXPECT_FLOAT_EQ(val4.x, 1.2);
  EXPECT_FLOAT_EQ(val4.y, 3.4);
  EXPECT_FLOAT_EQ(val4.z, 0.0);
  EXPECT_TRUE(utl::getAttribute(elem, "bar", val4, 2));
  EXPECT_FLOAT_EQ(val4.x, 1.2);
  EXPECT_FLOAT_EQ(val4.y, 3.4);
  EXPECT_FLOAT_EQ(val4.z, 0.0);
  EXPECT_TRUE(utl::getAttribute(elem, "bar", val4, 3));
  EXPECT_FLOAT_EQ(val4.x, 1.2);
  EXPECT_FLOAT_EQ(val4.y, 3.4);
  EXPECT_FLOAT_EQ(val4.z, 3.4);

  elem = doc.RootElement()->FirstChildElement("vecval3");
  ASSERT_TRUE(elem != nullptr);
  EXPECT_TRUE(utl::getAttribute(elem, "bar", val4));
  EXPECT_FLOAT_EQ(val4.x, 1.2);
  EXPECT_FLOAT_EQ(val4.y, 3.4);
  EXPECT_FLOAT_EQ(val4.z, 5.6);
  EXPECT_TRUE(utl::getAttribute(elem, "bar", val4, 2));
  EXPECT_FLOAT_EQ(val4.x, 1.2);
  EXPECT_FLOAT_EQ(val4.y, 3.4);
  EXPECT_FLOAT_EQ(val4.z, 0.0);
  EXPECT_TRUE(utl::getAttribute(elem, "bar", val4, 3));
  EXPECT_FLOAT_EQ(val4.x, 1.2);
  EXPECT_FLOAT_EQ(val4.y, 3.4);
  EXPECT_FLOAT_EQ(val4.z, 5.6);
}

TEST(TestUtilities, getDirs)
{
  int cmp1 = utl::getDirs(1);
  int cmp2 = utl::getDirs(2);
  int cmp3 = utl::getDirs(3);
  int cmp4 = utl::getDirs(4);

  EXPECT_EQ(cmp1, 1);
  EXPECT_EQ(cmp2, 12);
  EXPECT_EQ(cmp3, 123);
  EXPECT_EQ(cmp4, 1234);
}
