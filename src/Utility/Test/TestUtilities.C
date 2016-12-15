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
  for (int i = 0; i < 5; i++)
    ASSERT_EQ(values2[i], i+1);
}

TEST(TestUtilities, ParseKnots)
{
  std::string simple("xi 0 0.5 0.9 1.0");
  strtok(const_cast<char*>(simple.c_str())," ");
  std::vector<double> xi;
  utl::parseKnots(xi);
  ASSERT_EQ(xi.size(), 4U);
  ASSERT_FLOAT_EQ(xi[0], 0.0);
  ASSERT_FLOAT_EQ(xi[1], 0.5);
  ASSERT_FLOAT_EQ(xi[2], 0.9);
  ASSERT_FLOAT_EQ(xi[3], 1.0);

  std::string graded5("xi C5 79 0.01 2.0");
  strtok(const_cast<char*>(graded5.c_str())," ");
  xi.clear();
  utl::parseKnots(xi);
  EXPECT_EQ(xi.size(), 79U);
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
  ASSERT_STREQ(utl::readLine(str), "Blah blah");
  ASSERT_STREQ(utl::readLine(str), "Blah bluh");
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
  ASSERT_TRUE(doc.RootElement());

  const char* bTrue[4]  = { "booltrue" , "boolone" , "boolon" , "boolyes" };
  const char* bFalse[4] = { "boolfalse", "boolzero", "booloff", "boolno"  };

  int i;
  for (i = 0; i < 4; i++) {
    const TiXmlElement* elem = doc.RootElement()->FirstChildElement(bTrue[i]);
    ASSERT_TRUE(elem != nullptr);
    bool b = false;
    ASSERT_TRUE(utl::getAttribute(elem, "bar", b));
    ASSERT_TRUE(b);
  }

  for (i = 0; i < 4; i++) {
    const TiXmlElement* elem = doc.RootElement()->FirstChildElement(bFalse[i]);
    ASSERT_TRUE(elem != nullptr);
    bool b = false;
    ASSERT_TRUE(utl::getAttribute(elem, "bar", b));
    ASSERT_FALSE(b);
  }

  const TiXmlElement* elem = doc.RootElement()->FirstChildElement("intval");
  ASSERT_TRUE(elem != nullptr);
  int val1 = 0;
  ASSERT_TRUE(utl::getAttribute(elem, "bar", val1));
  ASSERT_EQ(val1, 1);
  size_t val2 = 0;
  ASSERT_TRUE(utl::getAttribute(elem, "bar", val2));
  ASSERT_EQ(val2, 1U);

  elem = doc.RootElement()->FirstChildElement("realval");
  ASSERT_TRUE(elem != nullptr);
  double val3 = 0.0;
  ASSERT_TRUE(utl::getAttribute(elem, "bar", val3));
  ASSERT_FLOAT_EQ(val3, 2.01);

  elem = doc.RootElement()->FirstChildElement("vecval1");
  ASSERT_TRUE(elem != nullptr);
  Vec3 val4;
  ASSERT_TRUE(utl::getAttribute(elem, "bar", val4));
  ASSERT_FLOAT_EQ(val4.x, 1.2);
  ASSERT_FLOAT_EQ(val4.y, 0.0);
  ASSERT_FLOAT_EQ(val4.z, 0.0);
  ASSERT_TRUE(utl::getAttribute(elem, "bar", val4, 2));
  ASSERT_FLOAT_EQ(val4.x, 1.2);
  ASSERT_FLOAT_EQ(val4.y, 1.2);
  ASSERT_FLOAT_EQ(val4.z, 0.0);
  ASSERT_TRUE(utl::getAttribute(elem, "bar", val4, 3));
  ASSERT_FLOAT_EQ(val4.x, 1.2);
  ASSERT_FLOAT_EQ(val4.y, 1.2);
  ASSERT_FLOAT_EQ(val4.z, 1.2);

  elem = doc.RootElement()->FirstChildElement("vecval2");
  ASSERT_TRUE(elem != nullptr);
  ASSERT_TRUE(utl::getAttribute(elem, "bar", val4));
  ASSERT_FLOAT_EQ(val4.x, 1.2);
  ASSERT_FLOAT_EQ(val4.y, 3.4);
  ASSERT_FLOAT_EQ(val4.z, 0.0);
  ASSERT_TRUE(utl::getAttribute(elem, "bar", val4, 2));
  ASSERT_FLOAT_EQ(val4.x, 1.2);
  ASSERT_FLOAT_EQ(val4.y, 3.4);
  ASSERT_FLOAT_EQ(val4.z, 0.0);
  ASSERT_TRUE(utl::getAttribute(elem, "bar", val4, 3));
  ASSERT_FLOAT_EQ(val4.x, 1.2);
  ASSERT_FLOAT_EQ(val4.y, 3.4);
  ASSERT_FLOAT_EQ(val4.z, 3.4);

  elem = doc.RootElement()->FirstChildElement("vecval3");
  ASSERT_TRUE(elem != nullptr);
  ASSERT_TRUE(utl::getAttribute(elem, "bar", val4));
  ASSERT_FLOAT_EQ(val4.x, 1.2);
  ASSERT_FLOAT_EQ(val4.y, 3.4);
  ASSERT_FLOAT_EQ(val4.z, 5.6);
  ASSERT_TRUE(utl::getAttribute(elem, "bar", val4, 2));
  ASSERT_FLOAT_EQ(val4.x, 1.2);
  ASSERT_FLOAT_EQ(val4.y, 3.4);
  ASSERT_FLOAT_EQ(val4.z, 0.0);
  ASSERT_TRUE(utl::getAttribute(elem, "bar", val4, 3));
  ASSERT_FLOAT_EQ(val4.x, 1.2);
  ASSERT_FLOAT_EQ(val4.y, 3.4);
  ASSERT_FLOAT_EQ(val4.z, 5.6);
}
