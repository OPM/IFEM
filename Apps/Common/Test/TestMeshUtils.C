//==============================================================================
//!
//! \file TestMeshUtils.C
//!
//! \date Feb 16 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for various mesh quality indicators
//!
//==============================================================================

#include "MeshUtils.h"
#include "SIM2D.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;


TEST_CASE("TestMeshUtils.Aspect2D")
{
  SIM2D model;
  REQUIRE(model.createDefaultModel());
  REQUIRE(model.createFEMmodel());

  std::vector<double> aspect;
  MeshUtils::computeAspectRatios(aspect, model);
  REQUIRE_THAT(aspect.front(), WithinRel(1.0));

  Vector tmp(2*model.getNoNodes());
  tmp(2) = -1.0;
  MeshUtils::computeAspectRatios(aspect, model, tmp);
  REQUIRE_THAT(aspect.front(), WithinRel(2.0));
}


TEST_CASE("TestMeshUtils.Skewness2D")
{
  SIM2D model;
  REQUIRE(model.createDefaultModel());
  REQUIRE(model.createFEMmodel());

  std::vector<double> skewness;
  MeshUtils::computeMeshSkewness(skewness, model);
  REQUIRE_THAT(skewness.front(), WithinAbs(0.0, 1e-13));

  Vector tmp(2*model.getNoNodes());
  tmp(1) = -1.0;
  MeshUtils::computeMeshSkewness(skewness, model, tmp);
  REQUIRE_THAT(skewness.front(), WithinRel(0.5));
}
