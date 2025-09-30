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

#include "Catch2Support.h"

#include <numeric>


TEST_CASE("TestTriangleQuadrature.Sanity")
{
  // Quadratures that are intended to not exist
  for (int i : {-1,2,5}) {
    REQUIRE(TriangleQuadrature::getCoord(i) == nullptr);
    REQUIRE(TriangleQuadrature::getWeight(i) == nullptr);
  }

  for (int i : {1,3,4}) {
    const double* c = TriangleQuadrature::getCoord(i);
    const double* w = TriangleQuadrature::getWeight(i);
    double sum = std::accumulate(w, w+i, 0.0);
    REQUIRE_THAT(sum, WithinRel(1.0));
    for (int j = 0; j < i; ++j){
      REQUIRE(c[j] >= 0.0);
      REQUIRE(c[j] <= 1.0);
    }
  }
}
