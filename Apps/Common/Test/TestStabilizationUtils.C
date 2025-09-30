//==============================================================================
//!
//! \file TestStabilizationUtils.C
//!
//! \date Oct 13 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for stabilization utilities.
//!
//==============================================================================

#include "StabilizationUtils.h"
#include "CoordinateMapping.h"

#include "Catch2Support.h"


TEST_CASE("TestStabilizationUtils.GetElementSize")
{
  std::vector<Vec3> XC;
  XC.push_back(Vec3(0.0, 0.0, 0.0));
  XC.push_back(Vec3(1.0, 0.0, 0.0));
  XC.push_back(Vec3(0.0, 0.45, 0.0));
  XC.push_back(Vec3(0.43, 1.0, 0.0));

  double result1 = StabilizationUtils::getElementSize(XC, 2);
  XC.clear();
  XC.push_back(Vec3(0.0, 0.0, 0.0));
  XC.push_back(Vec3(0.1, 0.0, 0.0));
  XC.push_back(Vec3(0.0, 0.45, 0.0));
  XC.push_back(Vec3(0.2, 1.0, 0.0));
  XC.push_back(Vec3(0.0, 0.0, 0.32));
  XC.push_back(Vec3(0.1, 0.0, 0.33));
  XC.push_back(Vec3(0.0, 0.45, 0.34));
  XC.push_back(Vec3(0.2, 1.0, 0.35));
  double result2 = StabilizationUtils::getElementSize(XC, 3);

  REQUIRE_THAT(result1, WithinRel(0.45));
  REQUIRE_THAT(result2, WithinRel(0.1));
}

TEST_CASE("TestStabilizationUtils.GetTauPt")
{
  Vector U(3);
  U.fill(1.0);
  Matrix J(3, 3);
  J.fill(2.0);
  Vector du(3);
  du.fill(3.0);

  Matrix G;
  utl::getGmat(J, &du[0], G);
  REQUIRE_THAT(StabilizationUtils::getTauPt(0.1, 0.3, U, G, 4.0, 5.0),
               WithinRel(0.032328788, 1e-7));
}


TEST_CASE("TestStabilizationUtils.GetTauNSPt")
{
  Vector U(3);
  U.fill(1.0);
  Matrix J(3, 3);
  J.fill(2.0);
  Vector du(3);
  du.fill(3.0);

  Matrix G;
  utl::getGmat(J, &du[0], G);

  double tauM, tauC;
  std::tie(tauM,tauC) = StabilizationUtils::getTauNSPt(0.1, 0.3, U, G);
  REQUIRE_THAT(tauM, WithinRel(0.016634906, 1e-7));
  REQUIRE_THAT(tauC, WithinRel(1.87858, 1e-7));
}


TEST_CASE("TestStabilizationUtils.GetTauNSALEPt")
{
  Vector U(3);
  U.fill(1.0);
  Matrix J(3, 3);
  J.fill(2.0);
  Vector du(3);
  du.fill(3.0);

  Matrix G;
  utl::getGmat(J, &du[0], G);

  double tauM, tauC;
  std::tie(tauM,tauC) = StabilizationUtils::getTauNSALEPt(0.1, 0.3, U, G);
  REQUIRE_THAT(tauM, WithinRel(0.016634906, 1e-7));
  REQUIRE_THAT(tauC, WithinRel(0.049904719, 1e-7));
}
