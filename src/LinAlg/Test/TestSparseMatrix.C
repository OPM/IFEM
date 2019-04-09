//==============================================================================
//!
//! \file TestSparseMatrix.C
//!
//! \date Apr 9 2019
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for sparse matrices
//!
//==============================================================================

#include "SparseMatrix.h"

#include "gtest/gtest.h"


TEST(TestSparseMatrix, CalcCSR)
{
  SparseMatrix Mat1(2,3);
  Mat1(1,1) = 1.0;
  Mat1(2,1) = 2.0;
  Mat1(2,3) = 3.0;
  IntVec IA1, JA1;
  SparseMatrix::calcCSR(IA1, JA1, 2, Mat1.getValues());

  EXPECT_EQ(IA1.size(), 3U);
  EXPECT_EQ(IA1[0], 0);
  EXPECT_EQ(IA1[1], 1);

  EXPECT_EQ(JA1.size(), 3U);
  EXPECT_EQ(JA1[0], 0);
  EXPECT_EQ(JA1[1], 0);
  EXPECT_EQ(JA1[2], 2);
}
