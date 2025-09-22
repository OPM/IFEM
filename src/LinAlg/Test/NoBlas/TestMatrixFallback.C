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


TEST(TestVectorFallback, Add)
{
  vectorAddTest<double>();
  vectorAddTest<float>();
}


TEST(TestVectorFallback, Dot)
{
  vectorDotTest<double>();
  vectorDotTest<float>();
}


TEST(TestVectorFallback, Multiply)
{
  vectorMultiplyTest<double>();
  vectorMultiplyTest<float>();
}


TEST(TestVectorFallback, Norm)
{
  vectorNormTest<double>();
  vectorNormTest<float>();
}


TEST(TestMatrixFallback, Multiply)
{
  multiplyTest<double>();
  multiplyTest<float>();
}


TEST(TestMatrixFallback, Norm)
{
  normTest<double>();
  normTest<float>();
}


TEST(TestMatrixFallback, OuterProduct)
{
  outerProductTest<double>();
  outerProductTest<float>();
}


TEST(TestMatrix3DFallback, Multiply)
{
  matrix3DMultiplyTest<double>();
  matrix3DMultiplyTest<float>();
}
