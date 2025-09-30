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
#include "Tensor.h"
#include "TensorFunction.h"
#include "Vec3Oper.h"
#include "tinyxml2.h"

#include "Catch2Support.h"


TEST_CASE("TestAnaSol.ParseDerivatives")
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

  tinyxml2::XMLDocument doc;
  doc.Parse(input);
  const tinyxml2::XMLElement* tag = doc.RootElement();
  REQUIRE(tag != nullptr);
  REQUIRE(strcmp(tag->Value(),"anasol") == 0);

  AnaSol mySol(tag,true);
  RealFunc* v = mySol.getScalarSol();
  REQUIRE(v != nullptr);

  Vec3 X(1.0,2.0,3.0);
  REQUIRE_THAT((*v)(X), WithinRel(8.0));
  REQUIRE_THAT(v->deriv(X,1), WithinRel(3.0));
  REQUIRE_THAT(v->deriv(X,2), WithinRel(4.0));
  REQUIRE_THAT(v->deriv(X,3), WithinRel(1.0));
  REQUIRE_THAT(v->dderiv(X,1,1), WithinRel(6.0));
  REQUIRE_THAT(v->dderiv(X,2,2), WithinRel(2.0));
  REQUIRE_THAT(v->dderiv(X,3,3), WithinRel(0.0));
  REQUIRE_THAT(v->dderiv(X,1,2), WithinRel(0.0));
  REQUIRE_THAT(v->dderiv(X,3,2), WithinRel(0.0));

  VecFunc* v2 = mySol.getScalarSecSol();
  REQUIRE(v2 != nullptr);
  REQUIRE((*v2)(X) == Vec3(3.0,4.0,1.0));
  REQUIRE(v2->deriv(X,1) == Vec3(6.0,0.0,0.0));
  REQUIRE(v2->deriv(X,2) == Vec3(0.0,2.0,0.0));
  REQUIRE(v2->dderiv(X,1,1) == Vec3(6.0,0.0,0.0));
}


TEST_CASE("TestAnaSol.ParseAD")
{
  const char* input = "<anasol type=\"expression\" autodiff=\"true\">"
    "  <primary>"
    "    x^3+y^2+z"
    "  </primary>"
    "</anasol>";

  tinyxml2::XMLDocument doc;
  doc.Parse(input);
  const tinyxml2::XMLElement* tag = doc.RootElement();
  REQUIRE(tag != nullptr);
  REQUIRE(strcmp(tag->Value(),"anasol") == 0);

  AnaSol mySol(tag,true);
  RealFunc* v = mySol.getScalarSol();
  REQUIRE(v != nullptr);

  mySol.setupSecondarySolutions();

  Vec3 X(1.0,2.0,3.0);
  const Vec3 grad1 = v->gradient(X);
  const SymmTensor hess = v->hessian(X);
  VecFunc* v2 = mySol.getScalarSecSol();
  REQUIRE(v2 != nullptr);
  const Vec3 grad2 = (*v2)(X);
  REQUIRE_THAT((*v)(X), WithinRel(8.0));
  REQUIRE_THAT(v->deriv(X,1), WithinRel(3.0));
  REQUIRE_THAT(grad1[0], WithinRel(3.0));
  REQUIRE_THAT(grad2[0], WithinRel(3.0));
  REQUIRE_THAT(v->deriv(X,2), WithinRel(4.0));
  REQUIRE_THAT(grad1[1], WithinRel(4.0));
  REQUIRE_THAT(grad2[1], WithinRel(4.0));
  REQUIRE_THAT(v->deriv(X,3), WithinRel(1.0));
  REQUIRE_THAT(grad1[2], WithinRel(1.0));
  REQUIRE_THAT(grad2[2], WithinRel(1.0));

  REQUIRE_THAT(v->dderiv(X,1,1), WithinRel(6.0));
  REQUIRE_THAT(hess(1,1), WithinRel(6.0));
  REQUIRE_THAT(v->dderiv(X,2,2), WithinRel(2.0));
  REQUIRE_THAT(hess(2,2), WithinRel(2.0));
  REQUIRE_THAT(v->dderiv(X,3,3), WithinRel(0.0));
  REQUIRE_THAT(hess(3,3), WithinRel(0.0));
  REQUIRE_THAT(v->dderiv(X,1,2), WithinRel(0.0));
  REQUIRE_THAT(hess(1,2), WithinRel(0.0));
  REQUIRE_THAT(v->dderiv(X,3,2), WithinRel(0.0));
  REQUIRE_THAT(hess(3,2), WithinRel(0.0));
}


TEST_CASE("TestAnaSol.ParseFD")
{
  const char* input = "<anasol type=\"expression\">"
    "  <primary>"
    "    x^3+y^2+z"
    "  </primary>"
    "</anasol>";

  tinyxml2::XMLDocument doc;
  doc.Parse(input);
  const tinyxml2::XMLElement* tag = doc.RootElement();
  REQUIRE(tag != nullptr);
  REQUIRE(strcmp(tag->Value(),"anasol") == 0);

  AnaSol mySol(tag,true);
  RealFunc* v = mySol.getScalarSol();
  REQUIRE(v != nullptr);

  mySol.setupSecondarySolutions();

  Vec3 X(1.0,2.0,3.0);
  const Vec3 grad1 = v->gradient(X);
  VecFunc* v2 = mySol.getScalarSecSol();
  REQUIRE(v2 != nullptr);
  const Vec3 grad2 = (*v2)(X);
  REQUIRE_THAT((*v)(X), WithinRel(8.0));
  REQUIRE_THAT(v->deriv(X,1), WithinRel(3.0, 1e-6));
  REQUIRE_THAT(grad1[0], WithinRel(3.0, 1e-6));
  REQUIRE_THAT(grad2[0], WithinRel(3.0, 1e-6));
  REQUIRE_THAT(v->deriv(X,2), WithinRel(4.0, 1e-6));
  REQUIRE_THAT(grad1[1], WithinRel(4.0, 1e-6));
  REQUIRE_THAT(grad2[1], WithinRel(4.0, 1e-6));
  REQUIRE_THAT(v->deriv(X,3), WithinRel(1.0, 1e-6));
  REQUIRE_THAT(grad1[2], WithinRel(1.0, 1e-6));
  REQUIRE_THAT(grad2[2], WithinRel(1.0, 1e-6));
}

TEST_CASE("TestAnaSol.ParseADStress")
{
  const char* input = "<anasol type=\"expression\" symmetric=\"true\" autodiff=\"true\">"
    "  <primary>"
    "    x^3+y^2+z | x^2+y^3+z"
    "  </primary>"
    "</anasol>";

  tinyxml2::XMLDocument doc;
  doc.Parse(input);
  const tinyxml2::XMLElement* tag = doc.RootElement();
  REQUIRE(tag != nullptr);
  REQUIRE(strcmp(tag->Value(),"anasol") == 0);

  AnaSol mySol(tag,false);
  REQUIRE(mySol.getScalarSol() == nullptr);
  REQUIRE(mySol.getVectorSol() != nullptr);

  mySol.setupSecondarySolutions();

  REQUIRE(mySol.getVectorSecSol() == nullptr);
  const STensorFunc* v = mySol.getStressSol();
  REQUIRE(v != nullptr);

  Vec3 X(1.0,2.0,3.0);
  const SymmTensor grad1 = (*v)(X);
  REQUIRE_THAT(grad1(1,1), WithinRel(3*X.x*X.x));
  REQUIRE_THAT(grad1(1,2), WithinRel(2*X.y));
  REQUIRE_THAT(grad1(2,1), WithinRel(2*X.y));
  REQUIRE_THAT(grad1(2,2), WithinRel(3*X.y*X.y));
}
