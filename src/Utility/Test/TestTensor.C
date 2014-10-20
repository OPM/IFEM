//==============================================================================
//!
//! \file TestTensor.C
//!
//! \date Oct 14 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for second-order tensors with some basic operations.
//!
//==============================================================================

#include "Tensor.h"

#include "gtest/gtest.h"

TEST(TestTensor, Shift)
{
  const double data[9] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };

  Tensor T2(2), t2(2);
  T2 = std::vector<double>(data,data+4);
  Tensor T3(3), t3(3);
  T3 = std::vector<double>(data,data+9);

  t2=T2;
  t2.shift(1);
  ASSERT_FLOAT_EQ(t2(1,1), 3.0);
  ASSERT_FLOAT_EQ(t2(1,2), 1.0);
  ASSERT_FLOAT_EQ(t2(2,1), 4.0);
  ASSERT_FLOAT_EQ(t2(2,2), 2.0);
  t2=T2;
  t2.shift(-1);
  ASSERT_FLOAT_EQ(t2(1,1), 3.0);
  ASSERT_FLOAT_EQ(t2(1,2), 1.0);
  ASSERT_FLOAT_EQ(t2(2,1), 4.0);
  ASSERT_FLOAT_EQ(t2(2,2), 2.0);
  t2=T2;
  t2.shift(2);
  ASSERT_FLOAT_EQ(t2(1,1), 1.0);
  ASSERT_FLOAT_EQ(t2(1,2), 3.0);
  ASSERT_FLOAT_EQ(t2(2,1), 2.0);
  ASSERT_FLOAT_EQ(t2(2,2), 4.0);
  t2=T2;
  t2.shift(-2);
  ASSERT_FLOAT_EQ(t2(1,1), 1.0);
  ASSERT_FLOAT_EQ(t2(1,2), 3.0);
  ASSERT_FLOAT_EQ(t2(2,1), 2.0);
  ASSERT_FLOAT_EQ(t2(2,2), 4.0);

  t3=T3;
  t3.shift(1);
  ASSERT_FLOAT_EQ(t3(1,1), 7.0);
  ASSERT_FLOAT_EQ(t3(1,2), 1.0);
  ASSERT_FLOAT_EQ(t3(1,3), 4.0);
  ASSERT_FLOAT_EQ(t3(2,1), 8.0);
  ASSERT_FLOAT_EQ(t3(2,2), 2.0);
  ASSERT_FLOAT_EQ(t3(2,3), 5.0);
  ASSERT_FLOAT_EQ(t3(3,1), 9.0);
  ASSERT_FLOAT_EQ(t3(3,2), 3.0);
  ASSERT_FLOAT_EQ(t3(3,3), 6.0);
  t3=T3;
  t3.shift(-1);
  ASSERT_FLOAT_EQ(t3(1,1), 4.0);
  ASSERT_FLOAT_EQ(t3(1,2), 7.0);
  ASSERT_FLOAT_EQ(t3(1,3), 1.0);
  ASSERT_FLOAT_EQ(t3(2,1), 5.0);
  ASSERT_FLOAT_EQ(t3(2,2), 8.0);
  ASSERT_FLOAT_EQ(t3(2,3), 2.0);
  ASSERT_FLOAT_EQ(t3(3,1), 6.0);
  ASSERT_FLOAT_EQ(t3(3,2), 9.0);
  ASSERT_FLOAT_EQ(t3(3,3), 3.0);
  t3=T3;
  t3.shift(2);
  ASSERT_FLOAT_EQ(t3(1,1), 4.0);
  ASSERT_FLOAT_EQ(t3(1,2), 7.0);
  ASSERT_FLOAT_EQ(t3(1,3), 1.0);
  ASSERT_FLOAT_EQ(t3(2,1), 5.0);
  ASSERT_FLOAT_EQ(t3(2,2), 8.0);
  ASSERT_FLOAT_EQ(t3(2,3), 2.0);
  ASSERT_FLOAT_EQ(t3(3,1), 6.0);
  ASSERT_FLOAT_EQ(t3(3,2), 9.0);
  ASSERT_FLOAT_EQ(t3(3,3), 3.0);
  t3=T3;
  t3.shift(-2);
  ASSERT_FLOAT_EQ(t3(1,1), 7.0);
  ASSERT_FLOAT_EQ(t3(1,2), 1.0);
  ASSERT_FLOAT_EQ(t3(1,3), 4.0);
  ASSERT_FLOAT_EQ(t3(2,1), 8.0);
  ASSERT_FLOAT_EQ(t3(2,2), 2.0);
  ASSERT_FLOAT_EQ(t3(2,3), 5.0);
  ASSERT_FLOAT_EQ(t3(3,1), 9.0);
  ASSERT_FLOAT_EQ(t3(3,2), 3.0);
  ASSERT_FLOAT_EQ(t3(3,3), 6.0);
}
