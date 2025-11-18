//==============================================================================
//!
//! \file TestSplineField.C
//!
//! \date Oct 6 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for scalar spline fields. 
//!
//==============================================================================

#include "ASMSquare.h"
#include "ASMCube.h"
#include "FieldTests.h"


TEST_CASE("TestSplineField.Value2D")
{
  Field2DTests<ASMSquare>::Value();
  Field2DTests<ASMSquare>::ValueQuad();
}


TEST_CASE("TestSplineField.Grad2D")
{
  Field2DTests<ASMSquare>::Grad();
}


TEST_CASE("TestSplineField.GradSepGeom2D")
{
  Field2DTests<ASMSquare>::GradSepGeom();
}


TEST_CASE("TestSplineField.Hessian2D")
{
  Field2DTests<ASMSquare>::Hessian();
}


TEST_CASE("TestSplineField.HessianSepGeom2D")
{
  Field2DTests<ASMSquare>::HessianSepGeom();
}


TEST_CASE("TestSplineField.Value3D")
{
  Field3DTests<ASMCube>::Value();
  Field3DTests<ASMCube>::ValueQuad();
}


TEST_CASE("TestSplineField.Grad3D")
{
  Field3DTests<ASMCube>::Grad();
}


TEST_CASE("TestSplineField.GradSepGeom3D")
{
  Field3DTests<ASMCube>::GradSepGeom();
}


TEST_CASE("TestSplineField.Hessian3D")
{
  Field3DTests<ASMCube>::Hessian();
}


TEST_CASE("TestSplineField.HessianSepGeom3D")
{
  Field3DTests<ASMCube>::HessianSepGeom();
}
