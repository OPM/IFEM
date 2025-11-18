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


TEST_CASE("TestSplineFields.Value2D")
{
  Fields2DTests<ASMSquare>::Value();
  Fields2DTests<ASMSquare>::ValueQuad();
}


TEST_CASE("TestSplineFields.Grad2D")
{
  Fields2DTests<ASMSquare>::Grad();
}


TEST_CASE("TestSplineFields.GradSepGeom2D")
{
  Fields2DTests<ASMSquare>::Grad();
}


TEST_CASE("TestSplineFields.Hessian2D")
{
  Fields2DTests<ASMSquare>::Hessian();
}


TEST_CASE("TestSplineFields.HessianSepGeom2D")
{
  Fields2DTests<ASMSquare>::Hessian();
}


TEST_CASE("TestSplineFields.Value2Dmx")
{
  Fields2DTests<ASMmxSquare>::Valuemx();
}


TEST_CASE("TestSplineFields.Grad2Dmx")
{
  Fields2DTests<ASMmxSquare>::Gradmx();
}


TEST_CASE("TestSplineFields.Hessian2Dmx")
{
  Fields2DTests<ASMmxSquare>::Hessianmx();
}


TEST_CASE("TestSplineFields.Value3D")
{
  Fields3DTests<ASMCube>::Value();
  Fields3DTests<ASMCube>::ValueQuad();
}


TEST_CASE("TestSplineFields.Grad3D")
{
  Fields3DTests<ASMCube>::Grad();
}


TEST_CASE("TestSplineFields.GradSepGeom3D")
{
  Fields3DTests<ASMCube>::GradSepGeom();
}


TEST_CASE("TestSplineFields.Hessian3D")
{
  Fields3DTests<ASMCube>::Hessian();
}


TEST_CASE("TestSplineFields.HessianSepGeom3D")
{
  Fields3DTests<ASMCube>::HessianSepGeom();
}


TEST_CASE("TestSplineFields.Value3Dmx")
{
  Fields3DTests<ASMmxCube>::Valuemx();
}


TEST_CASE("TestSplineFields.Grad3Dmx")
{
  Fields3DTests<ASMmxCube>::Gradmx();
}


TEST_CASE("TestSplineFields.Hessian3Dmx")
{
  Fields3DTests<ASMmxCube>::Hessianmx();
}
