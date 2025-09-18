//==============================================================================
//!
//! \file TestMatrixFallback.C
//!
//! \date Jan 28 2025
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for the fallback code for vector, matrix and matrix3d.
//!
//==============================================================================

#undef USE_MKL
#undef USE_ACCELERATE
#undef USE_CBLAS
#undef HAS_BLAS

#include "../MatrixTests.h"

#include <catch2/catch_template_test_macros.hpp>


TEMPLATE_TEST_CASE("TestVectorFallback.Add", "", float, double)
{
  vectorAddTest<TestType>();
}


TEMPLATE_TEST_CASE("TestVectorFallback.Dot", "", float, double)
{
  vectorDotTest<TestType>();
}


TEMPLATE_TEST_CASE("TestVectorFallback.Multiply", "", float, double)
{
  vectorMultiplyTest<TestType>();
}


TEMPLATE_TEST_CASE("TestVectorFallback.Norm", "", float, double)
{
  vectorNormTest<TestType>();
}


TEMPLATE_TEST_CASE("TestMatrixFallback.Multiply", "", float, double)
{
  multiplyTest<TestType>();
}


TEMPLATE_TEST_CASE("TestMatrixFallback.Norm", "", float, double)
{
  normTest<TestType>();
}


TEMPLATE_TEST_CASE("TestMatrixFallback.OuterProduct", "", float, double)
{
  outerProductTest<TestType>();
}


TEMPLATE_TEST_CASE("TestMatrix3DFallback.Multiply", "", float, double)
{
  matrix3DMultiplyTest<TestType>();
}
