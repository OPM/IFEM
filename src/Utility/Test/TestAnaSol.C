//==============================================================================
//!
//! \file TestAnaSol.C
//!
//! \date May 13 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Tests for parsing of analytical solutions.
//!
//==============================================================================

#include "AnaSol.h"
#include "Function.h"
#include "Vec3Oper.h"
#include "tinyxml.h"

#include "gtest/gtest.h"


TEST(TestAnaSol, ParseDerivatives)
{
  const char* input = "<anasol>"
    "  <primary>"
    "    x^3+y^2+z"
    "    <derivative dir=\"1\">3*x^2</derivative>"
    "    <derivative dir=\"2\">2*y</derivative>"
    "    <derivative dir=\"3\">1</derivative>"
    "    <derivative d1=\"1\" d2=\"1\">6*x</derivative>"
    "    <derivative d1=\"2\" d2=\"2\">2</derivative>"
    "  </primary>"
    "  <secondary>"
    "    3*x^2|2*y|1"
    "    <derivative dir=\"1\">6*x|0|0</derivative>"
    "    <derivative dir=\"2\">0|2|0</derivative>"
    "    <derivative d1=\"1\" d2=\"1\">6|0|0</derivative>"
    "  </secondary>"
    "</anasol>";

  TiXmlDocument doc;
  doc.Parse(input,nullptr,TIXML_ENCODING_UTF8);
  const TiXmlElement* tag = doc.RootElement();
  ASSERT_TRUE(tag != nullptr);
  ASSERT_EQ(strcmp(tag->Value(),"anasol"),0);

  AnaSol mySol(tag,true);
  RealFunc* v = mySol.getScalarSol();
  ASSERT_TRUE(v != nullptr);

  Vec3 X(1.0,2.0,3.0);
  EXPECT_FLOAT_EQ((*v)(X),8.0);
  EXPECT_FLOAT_EQ(v->deriv(X,1),3.0);
  EXPECT_FLOAT_EQ(v->deriv(X,2),4.0);
  EXPECT_FLOAT_EQ(v->deriv(X,3),1.0);
  EXPECT_FLOAT_EQ(v->dderiv(X,1,1),6.0);
  EXPECT_FLOAT_EQ(v->dderiv(X,2,2),2.0);
  EXPECT_FLOAT_EQ(v->dderiv(X,3,3),0.0);
  EXPECT_FLOAT_EQ(v->dderiv(X,1,2),0.0);
  EXPECT_FLOAT_EQ(v->dderiv(X,3,2),0.0);

  VecFunc* v2 = mySol.getScalarSecSol();
  ASSERT_TRUE(v2 != nullptr);
  EXPECT_TRUE((*v2)(X) == Vec3(3.0,4.0,1.0));
  EXPECT_TRUE(v2->deriv(X,1) == Vec3(6.0,0.0,0.0));
  EXPECT_TRUE(v2->deriv(X,2) == Vec3(0.0,2.0,0.0));
  EXPECT_TRUE(v2->dderiv(X,1,1) == Vec3(6.0,0.0,0.0));
}
