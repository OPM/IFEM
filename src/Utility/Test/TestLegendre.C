//==============================================================================
//!
//! \file TestLegendre.C
//!
//! \date Oct 10 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for various utility methods for Spectral elements.
//!
//==============================================================================

#include "Legendre.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;


#ifdef HAS_BLAS
TEST_CASE("TestLegendre.GL")
{
  std::vector<Real> w, p;
  REQUIRE(Legendre::GL(w, p, 3));
  REQUIRE_THAT(w[0], WithinRel(0.5555555555555561));
  REQUIRE_THAT(w[1], WithinRel(0.8888888888888888));
  REQUIRE_THAT(w[2], WithinRel(0.5555555555555561));
  REQUIRE_THAT(p[0], WithinRel(-0.7745966692414832));
  REQUIRE_THAT(p[1], WithinAbs(0.0, 1e-14));
  REQUIRE_THAT(p[2], WithinRel(0.7745966692414832));

  REQUIRE(Legendre::GL(w, p, 4));
  REQUIRE_THAT(w[0], WithinRel(0.3478548451374518));
  REQUIRE_THAT(w[1], WithinRel(0.6521451548625461));
  REQUIRE_THAT(w[2], WithinRel(0.6521451548625461));
  REQUIRE_THAT(w[3], WithinRel(0.3478548451374518));
  REQUIRE_THAT(p[0], WithinRel(-0.8611363115940535));
  REQUIRE_THAT(p[1], WithinRel(-0.3399810435848561));
  REQUIRE_THAT(p[2], WithinRel(0.3399810435848561));
  REQUIRE_THAT(p[3], WithinRel(0.8611363115940535));
}


TEST_CASE("TestLegendre.GLL")
{
  std::vector<Real> w, p;
  REQUIRE(Legendre::GLL(w, p, 3));
  REQUIRE_THAT(w[0], WithinRel(0.3333333333333333));
  REQUIRE_THAT(w[1], WithinRel(1.333333333333333));
  REQUIRE_THAT(w[2], WithinRel(0.3333333333333333));
  REQUIRE_THAT(p[0], WithinRel(-1.0));
  REQUIRE_THAT(p[1], WithinAbs(0.0, 1e-14));
  REQUIRE_THAT(p[2], WithinRel(1.0));

  REQUIRE(Legendre::GLL(w, p, 4));
  REQUIRE_THAT(w[0], WithinRel(0.1666666666666667));
  REQUIRE_THAT(w[1], WithinRel(0.8333333333333335));
  REQUIRE_THAT(w[2], WithinRel(0.8333333333333335));
  REQUIRE_THAT(w[3], WithinRel(0.1666666666666667));
  REQUIRE_THAT(p[0], WithinRel(-1.0));
  REQUIRE_THAT(p[1], WithinRel(-0.447213595499958));
  REQUIRE_THAT(p[2], WithinRel(0.447213595499958));
  REQUIRE_THAT(p[3], WithinRel(1.0));
}
#endif


TEST_CASE("TestLegendre.Eval")
{
  Real result1, result2;

  Legendre::eval(3, 0.3, result1);
  Legendre::eval(5, 0.7, result2);
  REQUIRE_THAT(result1, WithinRel(-0.3825));
  REQUIRE_THAT(result2, WithinRel(-0.3651987499999999));
}


TEST_CASE("TestLegendre.Derivative")
{
  Real result1, result2;

  Legendre::derEval(3, 0.3, result1);
  Legendre::derEval(5, 0.7, result2);
  REQUIRE_THAT(result1, WithinRel(-0.8250000000000001));
  REQUIRE_THAT(result2, WithinRel(-1.5335625));
}


#ifdef HAS_BLAS
TEST_CASE("TestLegendre.BasisDerivatives")
{
  utl::matrix<Real> result1, result2;

  Legendre::basisDerivatives(3, result1);
  const Real ref3[3][3] = {{-1.5,  2.0, -0.5},
                           {-0.5,  0.0,  0.5},
                           { 0.5, -2.0,  1.5}};

  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      REQUIRE_THAT(result1(i+1, j+1), WithinRel(ref3[i][j]));

  Legendre::basisDerivatives(10, result2);

  const Real ref10[10][10] =
    {{-22.5,                30.43814502928187,   -12.17794670742983,   6.943788485133953,
      -4.599354761103132,    3.294643033749184,   -2.452884175442686,  1.829563931903247,
      -1.275954836092663,    0.5},
     {-5.074064702978062,    0.0,                  7.185502869705838, -3.351663862746774,
       2.078207994036418,   -1.444948448751456,    1.059154463645442, -0.783239293137909,
       0.5437537382357057,  -0.2127027580091889},
     { 1.203351992852207,   -4.259297354965217,    0.0, 4.368674557010181,
      -2.104350179413155,    1.334915483878251,   -0.9366032131394465,
       0.6767970871960859,  -0.4642749589081571,   0.1807865854892499},
     {-0.528369376820273,    1.529902638181603,   -3.364125868297819, 0.0,
       3.387318101202445,   -1.64649408398706,     1.046189365502494,
      -0.7212373127216044,   0.4834623263339481,  -0.1866457893937361},
     { 0.3120472556084113,  -0.8458135734064247,   1.444850315601661,
      -3.020217958199347,    0.0,                  3.025188487751975, -1.468055509389994,
       0.9165551803364357,  -0.5880821430451693,   0.2235279447424539},
     {-0.2235279447424539,   0.5880821430451693,  -0.9165551803364357,
       1.468055509389994,   -3.025188487751975,    0.0, 3.020217958199347,
      -1.444850315601661,    0.8458135734064247,  -0.3120472556084113},
     { 0.1866457893937361,  -0.4834623263339481,   0.7212373127216044,
      -1.046189365502494,    1.64649408398706,    -3.387318101202445,
       0.0,                  3.364125868297819,   -1.529902638181603, 0.528369376820273},
     {-0.1807865854892499,   0.4642749589081571,  -0.6767970871960859,
       0.9366032131394465,  -1.334915483878251,    2.104350179413155,
      -4.368674557010181,    0.0,                  4.259297354965217,
      -1.203351992852207},
     { 0.2127027580091889,  -0.5437537382357057,   0.783239293137909,
      -1.059154463645442,    1.444948448751456,   -2.078207994036418,
       3.351663862746774,   -7.185502869705838,    0.0, 5.074064702978062},
     {-0.5,                  1.275954836092663,   -1.829563931903247,
       2.452884175442686,   -3.294643033749184,    4.599354761103132,
      -6.943788485133953,   12.17794670742983,   -30.43814502928187, 22.5}};

  for (size_t i = 0; i < 10; ++i)
    for (size_t j = 0; j < 10; ++j)
      REQUIRE_THAT(result2(i+1, j+1), WithinRel(ref10[i][j]));
}
#endif
