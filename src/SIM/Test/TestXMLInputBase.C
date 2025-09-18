//==============================================================================
//!
//! \file TestXMLInputBase.C
//!
//! \date Nov 28 2023
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for base class for XML input parsing.
//!
//==============================================================================

#include "XMLInputBase.h"

#include <catch2/catch_test_macros.hpp>

#include <string>
#include <vector>
#include <tinyxml2.h>

namespace {

class TestXMLInputBase : public XMLInputBase
{
public:
  bool parse(const tinyxml2::XMLElement* elem)
  {
    bool child = false;
    while (elem) {
      strings.push_back(elem->Value());
      const tinyxml2::XMLAttribute* attribute = elem->FirstAttribute();
      while (attribute)
      {
        strings.push_back(attribute->Name());
        strings.push_back(attribute->Value());
        attribute = attribute->Next();
      }
      if (elem->GetText())
        strings.push_back(elem->GetText());
      if (child)
        elem = elem->NextSiblingElement();
      else
        elem = elem->FirstChildElement();
      child = true;
    }
    return true;
  }

  std::vector<std::string> strings;
};

}


TEST_CASE("TestXMLInputBase.IncludeFiles")
{
  TestXMLInputBase x;
  x.readXML("src/SIM/Test/with_include.xml");

  const std::vector<std::string> ref = {
    "boundaryconditions",
    "dirichlet", "set", "foo", "comp", "12",
    "dirichlet", "set", "bar", "comp", "12", "type", "expression", "a*b*c",
    "neumann", "set", "foobar", "type", "constant", "1.0",
    "someothertag", "is_here"
  };

  REQUIRE(x.strings == ref);
}
