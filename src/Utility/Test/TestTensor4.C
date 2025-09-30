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

#include "Catch2Support.h"


TEST_CASE("TestTensor4.Constructor")
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
              REQUIRE_THAT( I(i,j,k,l), WithinRel(1.0));
              REQUIRE_THAT(Is(i,j,k,l), WithinRel(1.0));
            }
            else
            {
              REQUIRE_THAT( I(i,j,k,l), WithinAbs(0.0, 1e-14));
              REQUIRE_THAT(Is(i,j,k,l), WithinAbs(0.0, 1e-14));
            }

            if (i == j && k == l)
            {
              REQUIRE_THAT( J(i,j,k,l), WithinRel(1.0));
              REQUIRE_THAT(Js(i,j,k,l), WithinRel(1.2));
            }
            else
            {
              REQUIRE_THAT( J(i,j,k,l), WithinAbs(0.0, 1e-14));
              REQUIRE_THAT(Js(i,j,k,l), WithinAbs(0.0, 1e-14));
            }
          }
  }
}
