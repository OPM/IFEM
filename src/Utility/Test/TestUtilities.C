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

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <sstream>
#include <string>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;


TEST_CASE("TestUtilities.ParseIntegers")
{
  const char* simple = "1";
  const char* ranged = "1:5";
  const char* multip = "1 3 5:8 10";
  std::vector<int> values1, values2, values3;
  utl::parseIntegers(values1, simple);
  utl::parseIntegers(values2, ranged);
  utl::parseIntegers(values3, multip);

  int i;
  REQUIRE(values1.size() == 1);
  REQUIRE(values2.size() == 5);
  REQUIRE(values3.size() == 7);
  REQUIRE(values1.front() == 1);
  for (i = 0; i < 5; i++)
    REQUIRE(values2[i] == 1+i);
  for (i = 0; i < 3; i++)
    REQUIRE(values3[i] == 1+2*i);
  for (i = 2; i < 6; i++)
    REQUIRE(values3[i] == 3+i);
  REQUIRE(values3[6] == 10);
}


TEST_CASE("TestUtilities.ParseKnots")
{
  std::string simple("xi 0 0.5 0.9 1.0");
  strtok(const_cast<char*>(simple.c_str())," ");
  std::vector<double> xi;
  utl::parseKnots(xi);
  REQUIRE(xi.size() == 4);
  REQUIRE_THAT(xi[0], WithinAbs(0.0, 1e-14));
  REQUIRE_THAT(xi[1], WithinRel(0.5));
  REQUIRE_THAT(xi[2], WithinRel(0.9));
  REQUIRE_THAT(xi[3], WithinRel(1.0));

  std::string graded5("xi C5 79 0.01 2.0");
  strtok(const_cast<char*>(graded5.c_str())," ");
  xi.clear();
  utl::parseKnots(xi);
  REQUIRE(xi.size() == 79);
  REQUIRE_THAT(xi[39], WithinRel(0.5));
  REQUIRE_THAT(xi.front()+xi.back(), WithinRel(1.0));

  for (size_t i = 0; i < xi.size(); i++)
    std::cout <<"xi["<< i <<"] = "<< xi[i] << std::endl;
}


TEST_CASE("TestUtilities.ReadLine")
{
  using namespace std::string_literals;
  std::stringstream str;
  str << "Blah blah\n"
      << "# Commented line\n"
      << "Blah bluh\n";
  REQUIRE(utl::readLine(str) == "Blah blah"s);
  REQUIRE(utl::readLine(str) == "Blah bluh"s);
}


TEST_CASE("TestUtilities.IgnoreComments")
{
  using namespace std::string_literals;
  std::stringstream str;
  str << "Blah blah\n"
      << "# Commented line\n"
      << "Blah bluh\n";
  char tmp[1024];
  str.getline(tmp,1024);
  REQUIRE(tmp == "Blah blah"s);
  utl::ignoreComments(str);
  str.getline(tmp,1024);
  REQUIRE(tmp == "Blah bluh"s);
}


TEST_CASE("TestUtilities.GetAttribute")
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
  REQUIRE(doc.RootElement() != nullptr);

  const char* bTrue[4]  = { "booltrue" , "boolone" , "boolon" , "boolyes" };
  const char* bFalse[4] = { "boolfalse", "boolzero", "booloff", "boolno"  };

  int i;
  for (i = 0; i < 4; i++) {
    const tinyxml2::XMLElement* elem = doc.RootElement()->FirstChildElement(bTrue[i]);
    REQUIRE(elem != nullptr);
    bool b = false;
    REQUIRE(utl::getAttribute(elem, "bar", b));
    REQUIRE(b);
  }

  for (i = 0; i < 4; i++) {
    const tinyxml2::XMLElement* elem = doc.RootElement()->FirstChildElement(bFalse[i]);
    REQUIRE(elem != nullptr);
    bool b = false;
    REQUIRE(utl::getAttribute(elem, "bar", b));
    REQUIRE(!b);
  }

  const tinyxml2::XMLElement* elem = doc.RootElement()->FirstChildElement("intval");
  REQUIRE(elem != nullptr);
  int val1 = 0;
  REQUIRE(utl::getAttribute(elem, "bar", val1));
  REQUIRE(val1 == 1);
  size_t val2 = 0;
  REQUIRE(utl::getAttribute(elem, "bar", val2));
  REQUIRE(val2 == 1);

  elem = doc.RootElement()->FirstChildElement("realval");
  REQUIRE(elem != nullptr);
  double val3 = 0.0;
  REQUIRE(utl::getAttribute(elem, "bar", val3));
  REQUIRE_THAT(val3, WithinRel(2.01));

  elem = doc.RootElement()->FirstChildElement("vecval1");
  REQUIRE(elem != nullptr);
  Vec3 val4;
  REQUIRE(utl::getAttribute(elem, "bar", val4));
  REQUIRE_THAT(val4.x, WithinRel(1.2));
  REQUIRE_THAT(val4.y, WithinAbs(0.0, 1e-14));
  REQUIRE_THAT(val4.z, WithinAbs(0.0, 1e-14));
  REQUIRE(utl::getAttribute(elem, "bar", val4, 2));
  REQUIRE_THAT(val4.x, WithinRel(1.2));
  REQUIRE_THAT(val4.y, WithinRel(1.2));
  REQUIRE_THAT(val4.z, WithinAbs(0.0, 1e-14));
  REQUIRE(utl::getAttribute(elem, "bar", val4, 3));
  REQUIRE_THAT(val4.x, WithinRel(1.2));
  REQUIRE_THAT(val4.y, WithinRel(1.2));
  REQUIRE_THAT(val4.z, WithinRel(1.2));

  elem = doc.RootElement()->FirstChildElement("vecval2");
  REQUIRE(elem != nullptr);
  REQUIRE(utl::getAttribute(elem, "bar", val4));
  REQUIRE_THAT(val4.x, WithinRel(1.2));
  REQUIRE_THAT(val4.y, WithinRel(3.4));
  REQUIRE_THAT(val4.z, WithinAbs(0.0, 1e-14));
  REQUIRE(utl::getAttribute(elem, "bar", val4, 2));
  REQUIRE_THAT(val4.x, WithinRel(1.2));
  REQUIRE_THAT(val4.y, WithinRel(3.4));
  REQUIRE_THAT(val4.z, WithinAbs(0.0, 1e-14));
  REQUIRE(utl::getAttribute(elem, "bar", val4, 3));
  REQUIRE_THAT(val4.x, WithinRel(1.2));
  REQUIRE_THAT(val4.y, WithinRel(3.4));
  REQUIRE_THAT(val4.z, WithinRel(3.4));

  elem = doc.RootElement()->FirstChildElement("vecval3");
  REQUIRE(elem != nullptr);
  REQUIRE(utl::getAttribute(elem, "bar", val4));
  REQUIRE_THAT(val4.x, WithinRel(1.2));
  REQUIRE_THAT(val4.y, WithinRel(3.4));
  REQUIRE_THAT(val4.z, WithinRel(5.6));
  REQUIRE(utl::getAttribute(elem, "bar", val4, 2));
  REQUIRE_THAT(val4.x, WithinRel(1.2));
  REQUIRE_THAT(val4.y, WithinRel(3.4));
  REQUIRE_THAT(val4.z, WithinAbs(0.0, 1e-14));
  REQUIRE(utl::getAttribute(elem, "bar", val4, 3));
  REQUIRE_THAT(val4.x, WithinRel(1.2));
  REQUIRE_THAT(val4.y, WithinRel(3.4));
  REQUIRE_THAT(val4.z, WithinRel(5.6));
}


TEST_CASE("TestUtilities.GetDirs")
{
  int cmp1 = utl::getDirs(1);
  int cmp2 = utl::getDirs(2);
  int cmp3 = utl::getDirs(3);
  int cmp4 = utl::getDirs(4);

  REQUIRE(cmp1 == 1);
  REQUIRE(cmp2 == 12);
  REQUIRE(cmp3 == 123);
  REQUIRE(cmp4 == 1234);
}
