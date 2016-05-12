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

#include "IntegrandBase.h"
#include "MeshUtils.h"
#include "SIM2D.h"

#include "gtest/gtest.h"

class DummyIntegrand : public IntegrandBase {
};

class TestSIM : public SIM2D {
  public:
    TestSIM() : SIM2D(new DummyIntegrand, 2)
    {
      myModel = this->createDefaultGeometry(nullptr);
    }
};


TEST(TestMeshUtils, Aspect2D)
{
  TestSIM model;
  model.preprocess();

  Vector aspect;
  MeshUtils::computeAspectRatios(aspect, model);
  ASSERT_FLOAT_EQ(aspect.front(), 1.0);

  Vector tmp(model.getNoDOFs());
  tmp(2) = -1.0;
  MeshUtils::computeAspectRatios(aspect, model, tmp);
  ASSERT_FLOAT_EQ(aspect.front(), 2.0);
}


TEST(TestMeshUtils, Skewness2D)
{
  TestSIM model;
  model.preprocess();

  Vector skewness;
  MeshUtils::computeMeshSkewness(skewness, model);
  ASSERT_FLOAT_EQ(skewness.front(), 0.0);

  Vector tmp(model.getNoDOFs());
  tmp(1) = -1.0;
  MeshUtils::computeMeshSkewness(skewness, model, tmp);
  ASSERT_FLOAT_EQ(skewness.front(), 0.5);
}
