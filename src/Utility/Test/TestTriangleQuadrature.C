//==============================================================================
//!
//! \file TestTriangleQuadrature.C
//!
//! \date Jul 31 2019
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for triangle quadratures.
//!
//==============================================================================

#include "TriangleQuadrature.h"

#include "gtest/gtest.h"
#include <numeric>


TEST(TestTriangleQuadrature, Sanity)
{
  // Quadratures that are intended to not exist
  for (int i : {-1,2,5}) {
    EXPECT_EQ(TriangleQuadrature::getCoord(i), nullptr);
    EXPECT_EQ(TriangleQuadrature::getWeight(i), nullptr);
  }

  for (int i : {1,3,4}) {
    const double* c = TriangleQuadrature::getCoord(i);
    const double* w = TriangleQuadrature::getWeight(i);
    double sum = std::accumulate(w, w+i, 0.0);
    EXPECT_FLOAT_EQ(sum, 1.0);
    for (int j = 0; j < i; ++j)
      EXPECT_TRUE(c[j] >= 0.0 && c[j] <= 1.0);
  }
}
