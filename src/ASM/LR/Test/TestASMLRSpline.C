//==============================================================================
//!
//! \file TestASMLRSpline.C
//!
//! \date Oct 31 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for the base class of LR B-spline patches.
//!
//==============================================================================

#include "LR/ASMLRSpline.h"

#include "LRSpline/Basisfunction.h"

#include "Catch2Support.h"

#include <array>


TEST_CASE("TestASMLRSpline.SortDifferentOrders")
{
  // Two functions of order 3 in u and 2 in v, which differ only in the last
  // knot in u. They are sorted by all knots in u first, also when the order
  // in v is lower, which would make the v knots decide.
  const std::array<double,4> ua{0.0, 0.0, 0.0, 1.0}, ub{0.0, 0.0, 0.0, 0.5};
  const std::array<double,3> va{0.0, 0.0, 1.0}, vb{0.0, 1.0, 1.0};
  const std::array<double,2> cp{0.0, 0.0};
  LR::Basisfunction a(ua.begin(), va.begin(), cp.begin(), 2, 3, 2);
  LR::Basisfunction b(ub.begin(), vb.begin(), cp.begin(), 2, 3, 2);

  std::vector<LR::Basisfunction*> functions{&a, &b};
  ASMLRSpline::Sort(0, 1, 0, functions);
  REQUIRE(functions.size() == 2);
  CHECK(functions[0] == &b);
  CHECK(functions[1] == &a);

  // With the directions swapped, the v knots are compared first
  ASMLRSpline::Sort(0, 1, 4, functions);
  CHECK(functions[0] == &a);
  CHECK(functions[1] == &b);
}
