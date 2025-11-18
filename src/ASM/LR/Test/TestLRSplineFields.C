//==============================================================================
//!
//! \file TestLRSplineFields.C
//!
//! \date Mar 9 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for vector LR spline fields.
//!
//==============================================================================

#include "ASMuCube.h"
#include "ASMuSquare.h"
#include "Test/FieldTests.h"


TEST_CASE("TestLRSplineFields.Value2D")
{
  Fields2DTests<ASMuSquare>::Value();
  Fields2DTests<ASMuSquare>::ValueQuad();
}


TEST_CASE("TestLRSplineFields.Grad2D")
{
  Fields2DTests<ASMuSquare>::Grad();
}


TEST_CASE("TestLRSplineFields.GradSepGeom2D")
{
  Fields2DTests<ASMuSquare>::Grad();
}


TEST_CASE("TestLRSplineFields.Hessian2D")
{
  Fields2DTests<ASMuSquare>::Hessian();
}


TEST_CASE("TestLRSplineFields.HessianSepGeom2D")
{
  Fields2DTests<ASMuSquare>::Hessian();
}


TEST_CASE("TestLRSplineFields.Value2Dmx")
{
  Fields2DTests<ASMmxuSquare>::Valuemx();
}


TEST_CASE("TestLRSplineFields.Grad2Dmx")
{
  Fields2DTests<ASMmxuSquare>::Gradmx();
}


TEST_CASE("TestLRSplineFields.Hessian2Dmx")
{
  Fields2DTests<ASMmxuSquare>::Hessianmx();
}


TEST_CASE("TestLRSplineFields.Value3D")
{
  Fields3DTests<ASMuCube>::Value();
  Fields3DTests<ASMuCube>::ValueQuad();
}


TEST_CASE("TestLRSplineFields.Grad3D")
{
  Fields3DTests<ASMuCube>::Grad();
}


TEST_CASE("TestLRSplineFields.GradSepGeom3D")
{
  Fields3DTests<ASMuCube>::GradSepGeom();
}


TEST_CASE("TestLRSplineFields.Hessian3D")
{
  Fields3DTests<ASMuCube>::Hessian();
}


TEST_CASE("TestLRSplineFields.HessianSepGeom3D")
{
  Fields3DTests<ASMuCube>::HessianSepGeom();
}


TEST_CASE("TestLRSplineFields.Value3Dmx")
{
  Fields3DTests<ASMmxuCube>::Valuemx();
}


TEST_CASE("TestLRSplineFields.Grad3Dmx")
{
  Fields3DTests<ASMmxuCube>::Gradmx();
}


TEST_CASE("TestLRSplineFields.Hessian3Dmx")
{
  Fields3DTests<ASMmxuCube>::Hessianmx();
}
