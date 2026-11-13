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


namespace
{
  //! \brief Convenience class for parsing anasol tags.
  class AnaSolParser
  {
  public:
    //! \brief The constructor parses the specified XML input.
    explicit AnaSolParser(const char* xinp)
    {
      doc.Parse(xinp);
      const tinyxml2::XMLElement* tag = doc.RootElement();
      REQUIRE(tag != nullptr);
      REQUIRE(strcmp(tag->Value(),"anasol") == 0);
    }

    //! \brief Returns the root element.
    operator const tinyxml2::XMLElement*() { return doc.RootElement(); }

  private:
    tinyxml2::XMLDocument doc; //!< XML document holding the parsed input
  };
}


TEST_CASE("TestAnaSol.ParseDerivatives")
{
  AnaSolParser input
    ("<anasol>"
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
     "</anasol>");

  AnaSol mySol(input,true);

  RealFunc* v = mySol.getScalarSol();
  VecFunc* v2 = mySol.getScalarSecSol();
  REQUIRE(v != nullptr);
  REQUIRE(v2 != nullptr);

  const Vec3 X(1.0,2.0,3.0);
  CHECK_THAT((*v)(X), WithinRel(8.0));
  CHECK_THAT(v->deriv(X,1), WithinRel(3.0));
  CHECK_THAT(v->deriv(X,2), WithinRel(4.0));
  CHECK_THAT(v->deriv(X,3), WithinRel(1.0));
  CHECK_THAT(v->dderiv(X,1,1), WithinRel(6.0));
  CHECK_THAT(v->dderiv(X,2,2), WithinRel(2.0));
  CHECK_THAT(v->dderiv(X,3,3), WithinRel(0.0));
  CHECK_THAT(v->dderiv(X,1,2), WithinRel(0.0));
  CHECK_THAT(v->dderiv(X,3,2), WithinRel(0.0));

  CHECK((*v2)(X) == Vec3(3.0,4.0,1.0));
  CHECK(v2->deriv(X,1) == Vec3(6.0,0.0,0.0));
  CHECK(v2->deriv(X,2) == Vec3(0.0,2.0,0.0));
  CHECK(v2->dderiv(X,1,1) == Vec3(6.0,0.0,0.0));
}


TEST_CASE("TestAnaSol.ParseAD")
{
  AnaSolParser input
    ("<anasol type=\"expression\" autodiff=\"true\">"
     "  <primary>x^3+y^2+z</primary>"
     "</anasol>");

  AnaSol mySol(input,true);
  mySol.setupSecondarySolutions();

  RealFunc* v = mySol.getScalarSol();
  VecFunc* v2 = mySol.getScalarSecSol();
  REQUIRE(v != nullptr);
  REQUIRE(v2 != nullptr);

  const Vec3 X(1.0,2.0,3.0);
  const Vec3 grad1 = v->gradient(X);
  const SymmTensor hess = v->hessian(X);
  const Vec3 grad2 = (*v2)(X);
  CHECK_THAT((*v)(X), WithinRel(8.0));
  CHECK_THAT(v->deriv(X,1), WithinRel(3.0));
  CHECK_THAT(grad1[0], WithinRel(3.0));
  CHECK_THAT(grad2[0], WithinRel(3.0));
  CHECK_THAT(v->deriv(X,2), WithinRel(4.0));
  CHECK_THAT(grad1[1], WithinRel(4.0));
  CHECK_THAT(grad2[1], WithinRel(4.0));
  CHECK_THAT(v->deriv(X,3), WithinRel(1.0));
  CHECK_THAT(grad1[2], WithinRel(1.0));
  CHECK_THAT(grad2[2], WithinRel(1.0));

  CHECK_THAT(v->dderiv(X,1,1), WithinRel(6.0));
  CHECK_THAT(hess(1,1), WithinRel(6.0));
  CHECK_THAT(v->dderiv(X,2,2), WithinRel(2.0));
  CHECK_THAT(hess(2,2), WithinRel(2.0));
  CHECK_THAT(v->dderiv(X,3,3), WithinRel(0.0));
  CHECK_THAT(hess(3,3), WithinRel(0.0));
  CHECK_THAT(v->dderiv(X,1,2), WithinRel(0.0));
  CHECK_THAT(hess(1,2), WithinRel(0.0));
  CHECK_THAT(v->dderiv(X,3,2), WithinRel(0.0));
  CHECK_THAT(hess(3,2), WithinRel(0.0));
}


TEST_CASE("TestAnaSol.ParseFD")
{
  AnaSolParser input
    ("<anasol type=\"expression\">"
     "  <primary>x^3+y^2+z</primary>"
     "</anasol>");

  AnaSol mySol(input,true);
  mySol.setupSecondarySolutions();

  RealFunc* v = mySol.getScalarSol();
  VecFunc* v2 = mySol.getScalarSecSol();
  REQUIRE(v != nullptr);
  REQUIRE(v2 != nullptr);

  const Vec3 X(1.0,2.0,3.0);
  const Vec3 grad1 = v->gradient(X);
  const Vec3 grad2 = (*v2)(X);
  CHECK_THAT((*v)(X), WithinRel(8.0));
  CHECK_THAT(v->deriv(X,1), WithinRel(3.0, 1e-6));
  CHECK_THAT(grad1[0], WithinRel(3.0, 1e-6));
  CHECK_THAT(grad2[0], WithinRel(3.0, 1e-6));
  CHECK_THAT(v->deriv(X,2), WithinRel(4.0, 1e-6));
  CHECK_THAT(grad1[1], WithinRel(4.0, 1e-6));
  CHECK_THAT(grad2[1], WithinRel(4.0, 1e-6));
  CHECK_THAT(v->deriv(X,3), WithinRel(1.0, 1e-6));
  CHECK_THAT(grad1[2], WithinRel(1.0, 1e-6));
  CHECK_THAT(grad2[2], WithinRel(1.0, 1e-6));
}


TEST_CASE("TestAnaSol.ParseADStress")
{
  AnaSolParser input
    ("<anasol type=\"expression\" symmetric=\"true\" autodiff=\"true\">"
     "  <primary>x^3+y^2+z | x^2+y^3+z</primary>"
     "</anasol>");

  AnaSol mySol(input,false);
  mySol.setupSecondarySolutions();

  REQUIRE(mySol.getScalarSol() == nullptr);
  REQUIRE(mySol.getVectorSol() != nullptr);
  REQUIRE(mySol.getVectorSecSol() == nullptr);
  STensorFunc* v = mySol.getStressSol();
  REQUIRE(v != nullptr);

  const Vec3 X(1.0,2.0,3.0);
  const SymmTensor grad1 = (*v)(X);
  CHECK_THAT(grad1(1,1), WithinRel(3*X.x*X.x));
  CHECK_THAT(grad1(1,2), WithinRel(2*X.y));
  CHECK_THAT(grad1(2,1), WithinRel(2*X.y));
  CHECK_THAT(grad1(2,2), WithinRel(3*X.y*X.y));
}


TEST_CASE("TestAnaSol.ParseMultiple")
{
  AnaSolParser input
    ("<anasol>"
     "  <primary>x^3+y^2+z</primary>"
     "  <primary>x/4+y/2+z</primary>"
     "  <secondary>3*x^2 | 2*y | 1</secondary>"
     "  <secondary type=\"constant\">0.25 0.5 0.0</secondary>"
     "</anasol>");

  AnaSol mySol(input,true);

  RealFunc* f1 = mySol.getScalarSol(0);
  RealFunc* f2 = mySol.getScalarSol(1);
  VecFunc* vf1 = mySol.getScalarSecSol(0);
  VecFunc* vf2 = mySol.getScalarSecSol(1);
  REQUIRE(f1 != nullptr);
  REQUIRE(f2 != nullptr);
  REQUIRE(vf1 != nullptr);
  REQUIRE(vf2 != nullptr);

  const Vec3 X(1.0,2.0,3.0);
  CHECK_THAT((*f1)(X), WithinRel(8.0));
  CHECK_THAT((*f2)(X), WithinRel(4.25));
  CHECK((*vf1)(X).equal(Vec3(3.0,4.0,1.0)));
  CHECK((*vf2)(X).equal(Vec3(0.25,0.5,0.0)));
}
