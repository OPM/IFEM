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

#include "gtest/gtest.h"


TEST(TestMeshUtils, Aspect2D)
{
  SIM2D model;
  ASSERT_TRUE(model.createDefaultModel());
  ASSERT_TRUE(model.createFEMmodel());

  std::vector<double> aspect;
  MeshUtils::computeAspectRatios(aspect, model);
  ASSERT_FLOAT_EQ(aspect.front(), 1.0);

  Vector tmp(2*model.getNoNodes());
  tmp(2) = -1.0;
  MeshUtils::computeAspectRatios(aspect, model, tmp);
  ASSERT_FLOAT_EQ(aspect.front(), 2.0);
}


TEST(TestMeshUtils, Skewness2D)
{
  SIM2D model;
  ASSERT_TRUE(model.createDefaultModel());
  ASSERT_TRUE(model.createFEMmodel());

  std::vector<double> skewness;
  MeshUtils::computeMeshSkewness(skewness, model);
  ASSERT_FLOAT_EQ(skewness.front(), 0.0);

  Vector tmp(2*model.getNoNodes());
  tmp(1) = -1.0;
  MeshUtils::computeMeshSkewness(skewness, model, tmp);
  ASSERT_FLOAT_EQ(skewness.front(), 0.5);
}
