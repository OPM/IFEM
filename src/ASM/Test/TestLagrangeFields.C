//==============================================================================
//!
//! \file TestSplineFields.C
//!
//! \date Oct 6 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for vector spline fields.
//!
//==============================================================================

#include "ASMSquare.h"
#include "ASMCube.h"
#include "Test/FieldTests.h"


TEST_CASE("TestLagrangeFields.Value2D")
{
  Fields2DTests<ASMSquareLag>::Value();
  Fields2DTests<ASMSquareLag>::ValueQuad();
}


TEST_CASE("TestLagrangeFields.Grad2D")
{
  Fields2DTests<ASMSquareLag>::Grad();
}


TEST_CASE("TestLagrangeFields.Value3D")
{
  Fields3DTests<ASMCubeLag>::Value();
  Fields3DTests<ASMCubeLag>::ValueQuad();
}


TEST_CASE("TestLagrangeFields.Grad3D")
{
  Fields3DTests<ASMCubeLag>::Grad();
}
