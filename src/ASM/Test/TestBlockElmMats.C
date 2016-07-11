//==============================================================================
//!
//! \file TestBlockElmMats.C
//!
//! \date Jul 11 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Unit tests for block element matrices
//!
//==============================================================================

#include "BlockElmMats.h"

#include "gtest/gtest.h"


TEST(TestBlockElmMats, 1Basis2BlocksDiag)
{
  BlockElmMats mats(2);

  mats.resize(3, 3);
  ASSERT_TRUE(mats.redim(1, 2, 1));
  ASSERT_TRUE(mats.redim(2, 2, 1));
  mats.redimNewtonMat();

  mats.A[1].fill(1);
  mats.A[2].fill(2);

  const Matrix& N = mats.getNewtonMatrix();
  for (size_t b = 1; b <= 2; ++b)
    for (size_t i = 1; i <= 2; ++i)
      for (size_t j = 1; j <= 2 ; ++j)
        ASSERT_FLOAT_EQ(N(b+2*(i-1), b + 2*(j-1)), b);

  // check off-diagonal blocks
  for (size_t i = 1; i <= 2; ++i)
    for (size_t j = 1; j <= 2 ; ++j) {
      ASSERT_FLOAT_EQ(N(2*(i-1)+1, 2*j), 0.0);
      ASSERT_FLOAT_EQ(N(2*j, 2*(i-1)+1), 0.0);
    }
}


TEST(TestBlockElmMats, 1Basis2BlocksSymmetric)
{
  BlockElmMats mats(2);

  mats.resize(4, 3);
  ASSERT_TRUE(mats.redim(1, 2, 1));
  ASSERT_TRUE(mats.redim(2, 2, 1));
  ASSERT_TRUE(mats.redimOffDiag(3, 1));
  mats.redimNewtonMat();

  mats.A[1].fill(1);
  mats.A[2].fill(2);
  mats.A[3].fill(3);

  const Matrix& N = mats.getNewtonMatrix();

  // diagonal blocks
  for (size_t b = 1; b <= 2; ++b)
    for (size_t i = 1; i <= 2; ++i)
      for (size_t j = 1; j <= 2 ; ++j)
        ASSERT_FLOAT_EQ(N(b+2*(i-1), b + 2*(j-1)), b);

  // check off-diagonal blocks and symmetry
  for (size_t i = 1; i <= 2; ++i)
    for (size_t j = 1; j <= 2 ; ++j) {
      ASSERT_FLOAT_EQ(N(2*(i-1)+1, 2*j), 3);
      ASSERT_FLOAT_EQ(N(2*j, 2*(i-1)+1),
                      N(2*(i-1)+1, 2*j));
    }
}


TEST(TestBlockElmMats, 1Basis2BlocksSkewSymmetric)
{
  BlockElmMats mats(2);

  mats.resize(4, 3);
  ASSERT_TRUE(mats.redim(1, 2, 1));
  ASSERT_TRUE(mats.redim(2, 2, 1));
  ASSERT_TRUE(mats.redimOffDiag(3, -1));
  mats.redimNewtonMat();

  mats.A[1].fill(1);
  mats.A[2].fill(2);
  mats.A[3].fill(3);

  const Matrix& N = mats.getNewtonMatrix();

  // diagonal blocks
  for (size_t b = 1; b <= 2; ++b)
    for (size_t i = 1; i <= 2; ++i)
      for (size_t j = 1; j <= 2 ; ++j)
        ASSERT_FLOAT_EQ(N(b+2*(i-1), b + 2*(j-1)), b);

  // check off-diagonal blocks and symmetry
  for (size_t i = 1; i <= 2; ++i)
    for (size_t j = 1; j <= 2 ; ++j) {
      ASSERT_FLOAT_EQ(N(2*(i-1)+1, 2*j), 3);
      ASSERT_FLOAT_EQ(N(2*j, 2*(i-1)+1),
                      -N(2*(i-1)+1, 2*j));
    }
}


TEST(TestBlockElmMats, 1Basis2BlocksFull)
{
  BlockElmMats mats(2);

  mats.resize(5, 3);
  ASSERT_TRUE(mats.redim(1, 2, 1));
  ASSERT_TRUE(mats.redim(2, 2, 1));
  ASSERT_TRUE(mats.redimOffDiag(3, 0));
  ASSERT_TRUE(mats.redimOffDiag(4, 0));
  mats.redimNewtonMat();

  mats.A[1].fill(1);
  mats.A[2].fill(2);
  mats.A[3].fill(3);
  mats.A[4].fill(4);

  const Matrix& N = mats.getNewtonMatrix();

  // diagonal blocks
  for (size_t b = 1; b <= 2; ++b)
    for (size_t i = 1; i <= 2; ++i)
      for (size_t j = 1; j <= 2 ; ++j)
        ASSERT_FLOAT_EQ(N(b+2*(i-1), b + 2*(j-1)), b);

  // check off-diagonal blocks
  for (size_t i = 1; i <= 2; ++i)
    for (size_t j = 1; j <= 2 ; ++j) {
      ASSERT_FLOAT_EQ(N(2*(i-1)+1, 2*j), 3);
      ASSERT_FLOAT_EQ(N(2*j, 2*(i-1)+1), 4);
    }
}


TEST(TestBlockElmMats, 2Basis2BlocksDiag)
{
  BlockElmMats mats(2, 2);

  mats.resize(3, 3);
  ASSERT_TRUE(mats.redim(1, 2, 2, 1));
  ASSERT_TRUE(mats.redim(2, 2, 1, 2));
  mats.redimNewtonMat();

  mats.A[1].fill(1);
  mats.A[2].fill(2);

  const Matrix& N = mats.getNewtonMatrix();

  // check first basis blocks
  for (size_t i = 1; i <= 4; ++i)
    for (size_t j = 1; j <= 6; ++j)
      ASSERT_FLOAT_EQ(N(i,j), (j <= 4 ? 1.0 : 0.0));

  // check second basis blocks
  for (size_t i = 5; i <= 6; ++i)
    for (size_t j = 1; j <= 6; ++j)
      ASSERT_FLOAT_EQ(N(i,j), (j <= 4 ? 0.0 : 2.0));
}


TEST(TestBlockElmMats, 2Basis2BlocksSymmetric)
{
  BlockElmMats mats(2, 2);

  mats.resize(4, 3);
  ASSERT_TRUE(mats.redim(1, 2, 2, 1));
  ASSERT_TRUE(mats.redim(2, 2, 1, 2));
  ASSERT_TRUE(mats.redimOffDiag(3, 1));
  mats.redimNewtonMat();

  mats.A[1].fill(1);
  mats.A[2].fill(2);
  mats.A[3].fill(3);

  const Matrix& N = mats.getNewtonMatrix();

  // check first basis blocks
  for (size_t i = 1; i <= 4; ++i)
    for (size_t j = 1; j <= 6; ++j)
      ASSERT_FLOAT_EQ(N(i,j), (j <= 4 ? 1.0 : 3.0));

  // check second basis blocks
  for (size_t i = 5; i <= 6; ++i)
    for (size_t j = 1; j <= 6; ++j)
      ASSERT_FLOAT_EQ(N(i,j), (j <= 4 ? 3.0 : 2.0));
}


TEST(TestBlockElmMats, 2Basis2BlocksSkewSymmetric)
{
  BlockElmMats mats(2, 2);

  mats.resize(4, 3);
  ASSERT_TRUE(mats.redim(1, 2, 2, 1));
  ASSERT_TRUE(mats.redim(2, 2, 1, 2));
  ASSERT_TRUE(mats.redimOffDiag(3, -1));
  mats.redimNewtonMat();

  mats.A[1].fill(1);
  mats.A[2].fill(2);
  mats.A[3].fill(3);

  const Matrix& N = mats.getNewtonMatrix();

  // check first basis blocks
  for (size_t i = 1; i <= 4; ++i)
    for (size_t j = 1; j <= 6; ++j)
      ASSERT_FLOAT_EQ(N(i,j), (j <= 4 ? 1.0 : 3.0));

  // check second basis blocks
  for (size_t i = 5; i <= 6; ++i)
    for (size_t j = 1; j <= 6; ++j)
      ASSERT_FLOAT_EQ(N(i,j), (j <= 4 ? -3.0 : 2.0));
}


TEST(TestBlockElmMats, 2Basis2BlocksFull)
{
  BlockElmMats mats(2, 2);

  mats.resize(5, 3);
  ASSERT_TRUE(mats.redim(1, 2, 2, 1));
  ASSERT_TRUE(mats.redim(2, 2, 1, 2));
  ASSERT_TRUE(mats.redimOffDiag(3, 0));
  ASSERT_TRUE(mats.redimOffDiag(4, 0));
  mats.redimNewtonMat();

  mats.A[1].fill(1);
  mats.A[2].fill(2);
  mats.A[3].fill(3);
  mats.A[4].fill(4);

  const Matrix& N = mats.getNewtonMatrix();

  // check first basis blocks
  for (size_t i = 1; i <= 4; ++i)
    for (size_t j = 1; j <= 6; ++j)
      ASSERT_FLOAT_EQ(N(i,j), (j <= 4 ? 1.0 : 3.0));

  // check second basis blocks
  for (size_t i = 5; i <= 6; ++i)
    for (size_t j = 1; j <= 6; ++j)
      ASSERT_FLOAT_EQ(N(i,j), (j <= 4 ? 4.0 : 2.0));
}
