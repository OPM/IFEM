//==============================================================================
//!
//! \file TestTensor4.C
//!
//! \date Oct 29 2015
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Tests for fourth-order tensors with some basic operations.
//!
//==============================================================================

#include "Tensor4.h"

#include "gtest/gtest.h"


TEST(TestTensor4, Constructor)
{
  unsigned short int i, j, k, l, n;
  for (n = 1; n <= 3; n++)
  {
    SymmTensor4 I(n), J(n,true);
    Tensor4 Is(n), Js(n,1.2,true);
    for (i = 1; i <= n; i++)
      for (j = 1; j <= n; j++)
        for (k = 1; k <= n; k++)
          for (l = 1; l <= n; l++)
          {
            if (i == j && j == k && k == l)
            {
              ASSERT_FLOAT_EQ( I(i,j,k,l), 1.0);
              ASSERT_FLOAT_EQ(Is(i,j,k,l), 1.0);
            }
            else
            {
              ASSERT_FLOAT_EQ( I(i,j,k,l), 0.0);
              ASSERT_FLOAT_EQ(Is(i,j,k,l), 0.0);
            }

            if (i == j && k == l)
            {
              ASSERT_FLOAT_EQ( J(i,j,k,l), 1.0);
              ASSERT_FLOAT_EQ(Js(i,j,k,l), 1.2);
            }
            else
            {
              ASSERT_FLOAT_EQ( J(i,j,k,l), 0.0);
              ASSERT_FLOAT_EQ(Js(i,j,k,l), 0.0);
            }
          }
  }
}
