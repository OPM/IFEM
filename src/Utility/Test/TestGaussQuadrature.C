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

#include "gtest/gtest.h"
#include <cmath>


class TestGaussQuadrature : public testing::Test,
                            public testing::WithParamInterface<int>
{
};


TEST_P(TestGaussQuadrature, Integrate)
{
  const double* xi = GaussQuadrature::getCoord(GetParam());
  const double* wi = GaussQuadrature::getWeight(GetParam());

  double res = 0.0;
  for (int i = 0; i < GetParam(); ++i)
    for (int p = 0; p < GetParam(); ++p)
      res += wi[i]*pow(xi[i], p); // 1+x+x^2+x^3+..

  double exa = 0.0;
  for (int p = 0; p < GetParam(); ++p)
    exa += 1.0/(p+1) * (pow(1.0, p+1) - pow(-1.0, p+1));

  EXPECT_FLOAT_EQ(res, exa);
}


INSTANTIATE_TEST_CASE_P(TestGaussQuadrature,
                        TestGaussQuadrature,
                        testing::Values(1, 2, 3, 4, 5, 6, 7, 8, 9, 10));
