//==============================================================================
//!
//! \file TestFiniteElement.C
//!
//! \date Nov 17 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Unit tests for FiniteElement.
//!
//==============================================================================

#include "FiniteElement.h"

#include "Catch2Support.h"


TEST_CASE("TestFiniteElement.Print")
{
  FiniteElement   fe1(4);
  MxFiniteElement fe2({3,2,4});

  std::cout << fe1 << std::endl;
  std::cout << fe2 << std::endl;

  REQUIRE(fe1.N.size() == 4);
  REQUIRE(fe2.basis(1).size() == 3);
  REQUIRE(fe2.basis(2).size() == 2);
  REQUIRE(fe2.basis(3).size() == 4);
}
