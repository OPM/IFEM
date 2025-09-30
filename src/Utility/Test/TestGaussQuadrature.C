//==============================================================================
//!
//! \file TestGaussQuadrature.C
//!
//! \date Dec 31 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for Gaussian quadratures.
//!
//==============================================================================

#include "GaussQuadrature.h"

#include <cmath>

#include "Catch2Support.h"


TEST_CASE("TestGaussQuadrature.Integrate")
{
  auto param = GENERATE(1,2,3,4,5,6,7,8,9,10);

  SECTION("points =  " + std::to_string(param)) {
    const double* xi = GaussQuadrature::getCoord(param);
    const double* wi = GaussQuadrature::getWeight(param);

    double res = 0.0;
    for (int i = 0; i < param; ++i)
      for (int p = 0; p < param; ++p)
        res += wi[i]*pow(xi[i], p); // 1+x+x^2+x^3+..

    double exa = 0.0;
    for (int p = 0; p < param; ++p)
      exa += 1.0/(p+1) * (pow(1.0, p+1) - pow(-1.0, p+1));

    REQUIRE_THAT(res, WithinRel(exa, param == 10 ? 1e-7 : 1e-14));
  }
}
