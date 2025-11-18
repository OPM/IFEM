//==============================================================================
//!
//! \file TestLRSplineField.C
//!
//! \date Mar 9 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for scalar LR spline fields.
//!
//==============================================================================

#include "ASMuCube.h"
#include "ASMuSquare.h"
#include "Test/FieldTests.h"


TEST_CASE("TestLRSplineField.Value2D")
{
  Field2DTests<ASMuSquare>::Value();
  Field2DTests<ASMuSquare>::ValueQuad();
}


TEST_CASE("TestLRSplineField.Grad2D")
{
  Field2DTests<ASMuSquare>::Grad();
}


TEST_CASE("TestLRSplineField.GradSepGeom2D")
{
  Field2DTests<ASMuSquare>::GradSepGeom();
}


TEST_CASE("TestLRSplineField.Hessian2D")
{
  Field2DTests<ASMuSquare>::Hessian();
}


TEST_CASE("TestLRSplineField.HessianSepGeom2D")
{
  Field2DTests<ASMuSquare>::HessianSepGeom();
}


TEST_CASE("TestLRSplineField.Value3D")
{
  Field3DTests<ASMuCube>::Value();
  Field3DTests<ASMuCube>::ValueQuad();
}


TEST_CASE("TestLRSplineField.Grad3D")
{
  Field3DTests<ASMuCube>::Grad();
}


TEST_CASE("TestLRSplineField.GradSepGeom3D")
{
  Field3DTests<ASMuCube>::GradSepGeom();
}


TEST_CASE("TestLRSplineField.Hessian3D")
{
  Field3DTests<ASMuCube>::Hessian();
}


TEST_CASE("TestLRSplineField.HessianSepGeom3D")
{
  Field3DTests<ASMuCube>::HessianSepGeom();
}
