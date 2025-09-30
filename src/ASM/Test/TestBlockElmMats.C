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

#include "Catch2Support.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;


TEST_CASE("TestBlockElmMats.1Basis2BlocksDiag")
{
  BlockElmMats mats(2, 1);

  mats.resize(3, 3);
  REQUIRE(mats.redim(1, 2, 1));
  REQUIRE(mats.redim(2, 2, 1));
  mats.finalize();

  mats.A[1].fill(1);
  mats.A[2].fill(2);

  const Matrix& N = mats.getNewtonMatrix();
  for (size_t b = 1; b <= 2; ++b)
    for (size_t i = 1; i <= 2; ++i)
      for (size_t j = 1; j <= 2 ; ++j)
        REQUIRE_THAT(N(b+2*(i-1), b + 2*(j-1)),
                     WithinRel(static_cast<Real>(b)));

  // check off-diagonal blocks
  for (size_t i = 1; i <= 2; ++i)
    for (size_t j = 1; j <= 2 ; ++j) {
      REQUIRE_THAT(N(2*(i-1)+1, 2*j), WithinAbs(0.0, 1e-14));
      REQUIRE_THAT(N(2*j, 2*(i-1)+1), WithinAbs(0.0, 1e-14));
    }
}


TEST_CASE("TestBlockElmMats.1Basis2BlocksSymmetric")
{
  BlockElmMats mats(2, 1);

  mats.resize(4, 3);
  REQUIRE(mats.redim(1, 2, 1));
  REQUIRE(mats.redim(2, 2, 1));
  REQUIRE(mats.redimOffDiag(3, 1));
  mats.finalize();

  mats.A[1].fill(1);
  mats.A[2].fill(2);
  mats.A[3].fill(3);

  const Matrix& N = mats.getNewtonMatrix();

  // diagonal blocks
  for (size_t b = 1; b <= 2; ++b)
    for (size_t i = 1; i <= 2; ++i)
      for (size_t j = 1; j <= 2 ; ++j)
        REQUIRE_THAT(N(b+2*(i-1), b + 2*(j-1)),
                     WithinRel(static_cast<Real>(b)));

  // check off-diagonal blocks and symmetry
  for (size_t i = 1; i <= 2; ++i)
    for (size_t j = 1; j <= 2 ; ++j) {
      REQUIRE_THAT(N(2*(i-1)+1, 2*j), WithinRel(3.0));
      REQUIRE_THAT(N(2*j, 2*(i-1)+1), WithinRel(N(2*(i-1)+1, 2*j)));
    }
}


TEST_CASE("TestBlockElmMats.1Basis2BlocksSkewSymmetric")
{
  BlockElmMats mats(2, 1);

  mats.resize(4, 3);
  REQUIRE(mats.redim(1, 2, 1));
  REQUIRE(mats.redim(2, 2, 1));
  REQUIRE(mats.redimOffDiag(3, -1));
  mats.finalize();

  mats.A[1].fill(1);
  mats.A[2].fill(2);
  mats.A[3].fill(3);

  const Matrix& N = mats.getNewtonMatrix();

  // diagonal blocks
  for (size_t b = 1; b <= 2; ++b)
    for (size_t i = 1; i <= 2; ++i)
      for (size_t j = 1; j <= 2 ; ++j)
        REQUIRE_THAT(N(b+2*(i-1), b + 2*(j-1)),
                     WithinRel(static_cast<Real>(b)));

  // check off-diagonal blocks and symmetry
  for (size_t i = 1; i <= 2; ++i)
    for (size_t j = 1; j <= 2 ; ++j) {
      REQUIRE_THAT(N(2*(i-1)+1, 2*j), WithinRel(3.0));
      REQUIRE_THAT(N(2*j, 2*(i-1)+1),
                   WithinRel(-N(2*(i-1)+1, 2*j)));
    }
}


TEST_CASE("TestBlockElmMats.1Basis2BlocksFull")
{
  BlockElmMats mats(2, 1);

  mats.resize(5, 3);
  REQUIRE(mats.redim(1, 2, 1));
  REQUIRE(mats.redim(2, 2, 1));
  REQUIRE(mats.redimOffDiag(3, 0));
  REQUIRE(mats.redimOffDiag(4, 0));
  mats.finalize();

  mats.A[1].fill(1);
  mats.A[2].fill(2);
  mats.A[3].fill(3);
  mats.A[4].fill(4);

  const Matrix& N = mats.getNewtonMatrix();

  // diagonal blocks
  for (size_t b = 1; b <= 2; ++b)
    for (size_t i = 1; i <= 2; ++i)
      for (size_t j = 1; j <= 2 ; ++j)
        REQUIRE_THAT(N(b+2*(i-1), b + 2*(j-1)),
                     WithinRel(static_cast<Real>(b)));

  // check off-diagonal blocks
  for (size_t i = 1; i <= 2; ++i)
    for (size_t j = 1; j <= 2 ; ++j) {
      REQUIRE_THAT(N(2*(i-1)+1, 2*j), WithinRel(3.0));
      REQUIRE_THAT(N(2*j, 2*(i-1)+1), WithinRel(4.0));
    }
}


TEST_CASE("TestBlockElmMats.2Basis2BlocksDiag")
{
  BlockElmMats mats(2, 2);

  mats.resize(3, 3);
  REQUIRE(mats.redim(1, 2, 2, 1));
  REQUIRE(mats.redim(2, 2, 1, 2));
  mats.finalize();

  mats.A[1].fill(1);
  mats.A[2].fill(2);

  const Matrix& N = mats.getNewtonMatrix();

  // check first basis blocks
  for (size_t i = 1; i <= 4; ++i)
    for (size_t j = 1; j <= 6; ++j)
      if (j <= 4)
        REQUIRE_THAT(N(i,j), WithinRel(1.0));
      else
        REQUIRE_THAT(N(i,j), WithinAbs(0.0, 1e-14));

  // check second basis blocks
  for (size_t i = 5; i <= 6; ++i)
    for (size_t j = 1; j <= 6; ++j)
      if (j <= 4)
        REQUIRE_THAT(N(i,j), WithinAbs(0.0, 1e-14));
      else
        REQUIRE_THAT(N(i,j), WithinRel(2.0));
}


TEST_CASE("TestBlockElmMats.2Basis2BlocksSymmetric")
{
  BlockElmMats mats(2, 2);

  mats.resize(4, 3);
  REQUIRE(mats.redim(1, 2, 2, 1));
  REQUIRE(mats.redim(2, 2, 1, 2));
  REQUIRE(mats.redimOffDiag(3, 1));
  mats.finalize();

  mats.A[1].fill(1);
  mats.A[2].fill(2);
  mats.A[3].fill(3);

  const Matrix& N = mats.getNewtonMatrix();

  // check first basis blocks
  for (size_t i = 1; i <= 4; ++i)
    for (size_t j = 1; j <= 6; ++j)
      REQUIRE_THAT(N(i,j), WithinRel(j <= 4 ? 1.0 : 3.0));

  // check second basis blocks
  for (size_t i = 5; i <= 6; ++i)
    for (size_t j = 1; j <= 6; ++j)
      REQUIRE_THAT(N(i,j), WithinRel(j <= 4 ? 3.0 : 2.0));
}


TEST_CASE("TestBlockElmMats.2Basis2BlocksSkewSymmetric")
{
  BlockElmMats mats(2, 2);

  mats.resize(4, 3);
  REQUIRE(mats.redim(1, 2, 2, 1));
  REQUIRE(mats.redim(2, 2, 1, 2));
  REQUIRE(mats.redimOffDiag(3, -1));
  mats.finalize();

  mats.A[1].fill(1);
  mats.A[2].fill(2);
  mats.A[3].fill(3);

  const Matrix& N = mats.getNewtonMatrix();

  // check first basis blocks
  for (size_t i = 1; i <= 4; ++i)
    for (size_t j = 1; j <= 6; ++j)
      REQUIRE_THAT(N(i,j), WithinRel(j <= 4 ? 1.0 : 3.0));

  // check second basis blocks
  for (size_t i = 5; i <= 6; ++i)
    for (size_t j = 1; j <= 6; ++j)
      REQUIRE_THAT(N(i,j), WithinRel(j <= 4 ? -3.0 : 2.0));
}


TEST_CASE("TestBlockElmMats.2Basis2BlocksFull")
{
  BlockElmMats mats(2, 2);

  mats.resize(5, 3);
  REQUIRE(mats.redim(1, 2, 2, 1));
  REQUIRE(mats.redim(2, 2, 1, 2));
  REQUIRE(mats.redimOffDiag(3, 0));
  REQUIRE(mats.redimOffDiag(4, 0));
  mats.finalize();

  mats.A[1].fill(1);
  mats.A[2].fill(2);
  mats.A[3].fill(3);
  mats.A[4].fill(4);

  const Matrix& N = mats.getNewtonMatrix();

  // check first basis blocks
  for (size_t i = 1; i <= 4; ++i)
    for (size_t j = 1; j <= 6; ++j)
      REQUIRE_THAT(N(i,j), WithinRel(j <= 4 ? 1.0 : 3.0));

  // check second basis blocks
  for (size_t i = 5; i <= 6; ++i)
    for (size_t j = 1; j <= 6; ++j)
      REQUIRE_THAT(N(i,j), WithinRel(j <= 4 ? 4.0 : 2.0));
}
