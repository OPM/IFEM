//==============================================================================
//!
//! \file TestLagrangeField.C
//!
//! \date Nov 18 2025
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for scalar lagrange fields.
//!
//==============================================================================

#include "ASMSquare.h"
#include "ASMCube.h"
#include "FieldTests.h"


TEST_CASE("TestLagrangeField.Value2D")
{
  Field2DTests<ASMSquareLag>::Value();
  Field2DTests<ASMSquareLag>::ValueQuad();
}


TEST_CASE("TestLagrangeField.Grad2D")
{
  Field2DTests<ASMSquareLag>::Grad();
}


TEST_CASE("TestLagrangeField.Value3D")
{
  Field3DTests<ASMCubeLag>::Value();
  Field3DTests<ASMCubeLag>::ValueQuad();
}


TEST_CASE("TestLagrangeField.Grad3D")
{
  Field3DTests<ASMCubeLag>::Grad();
}
