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

#include "gtest/gtest.h"


TEST(TestFiniteElement, Print)
{
  FiniteElement   fe1(4);
  MxFiniteElement fe2({3,2,4});

  std::cout << fe1 << std::endl;
  std::cout << fe2 << std::endl;

  ASSERT_EQ((int)fe1.N.size(),4);
  ASSERT_EQ((int)fe2.basis(1).size(),3);
  ASSERT_EQ((int)fe2.basis(2).size(),2);
  ASSERT_EQ((int)fe2.basis(3).size(),4);
}
