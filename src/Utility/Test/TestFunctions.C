//==============================================================================
//!
//! \file TestFunctions.C
//!
//! \date Apr 28 2018
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Tests for parsing of functions.
//!
//==============================================================================

#include "Functions.h"
#include <cstdlib>
#include <cmath>

#include "gtest/gtest.h"


TEST(TestScalarFunc, ParseDerivative)
{
  const char* func1 = "sin(1.5*t)*t";
  const char* func2 = "sin(1.5*t)*t:1.5*cos(1.5*t)*t+sin(1.5*t)";

  ScalarFunc* f1 = utl::parseTimeFunc(func1,"expression");
  ScalarFunc* f2 = utl::parseTimeFunc(func2,"expression");

  ASSERT_TRUE(f1 != nullptr);
  ASSERT_TRUE(f2 != nullptr);

  double t = 0.0;
  for (int i = 0; i < 20; i++)
  {
    t += 0.314*(double)random()/(double)RAND_MAX;
    std::cout <<"f("<< t <<") = "<< (*f1)(t)
	      <<"  f'("<< t <<") = "<< f1->deriv(t) << std::endl;
    EXPECT_FLOAT_EQ((*f1)(t),sin(1.5*t)*t);
    EXPECT_FLOAT_EQ((*f2)(t),sin(1.5*t)*t);
    EXPECT_FLOAT_EQ(f1->deriv(t),1.5*cos(1.5*t)*t+sin(1.5*t));
    EXPECT_FLOAT_EQ(f2->deriv(t),1.5*cos(1.5*t)*t+sin(1.5*t));
  }
}
