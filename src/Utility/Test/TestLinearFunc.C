//==============================================================================
//!
//! \file TestLinearFunc.C
//!
//! \date Jul 21 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Tests for evaluation of piece-wise linear functions.
//!
//==============================================================================

#include "Functions.h"
#include "Vec3Oper.h"

#include "gtest/gtest.h"

TEST(TestLinearFunc, Evaluate)
{
  LinearFunc f1("src/Utility/Test/refdata/func.dat",2);
  EXPECT_EQ(f1(0.15), 2.5);
  EXPECT_EQ(f1(0.4),  2.5);

  LinVecFunc f2("src/Utility/Test/refdata/func.dat",3);
  Vec3 v1 = f2(0.15);
  Vec3 v2 = f2(0.2);
  Vec3 v3 = f2(0.0);
  std::cout <<"v1 = "<< v1 << std::endl;
  std::cout <<"v2 = "<< v2 << std::endl;
  std::cout <<"v3 = "<< v3 << std::endl;
  EXPECT_EQ(v1.y, 1.2);
  EXPECT_EQ(v2.x, 1.6);
  EXPECT_EQ(v3.z, 1.3);

  std::cout <<"Parsing time function: ";
  ScalarFunc* f3 = utl::parseTimeFunc("src/Utility/Test/refdata/func.dat",
                                      "PiecewiseLinear");
  std::cout <<"Parsing time function: ";
  ScalarFunc* f4 = utl::parseTimeFunc("src/Utility/Test/refdata/func.dat 4 2.0",
                                      "PiecewiseLinear");
  std::cout <<"Parsing vector-valued time function: ";
  VecTimeFunc* f5 = utl::parseVecTimeFunc("src/Utility/Test/refdata/func.dat 3",
                                          "PiecewiseLinear");
  ASSERT_FALSE(f3 == nullptr);
  EXPECT_NEAR((*f3)(0.12),2.2,1.0e-15);
  ASSERT_FALSE(f4 == nullptr);
  EXPECT_NEAR((*f4)(0.35),4.3,1.0e-15);
  EXPECT_EQ((*f5)(0.15),v1);
  EXPECT_EQ((*f5)(0.20),v2);
  EXPECT_EQ((*f5)(0.00),v3);
  delete f3;
  delete f4;
  delete f5;
}
