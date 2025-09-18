// $Id$
//==============================================================================
//!
//! \file TestDiagMatrix.C
//!
//! \date Aug 30 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Unit tests for diagonal system matrices.
//!
//==============================================================================

#include "DiagMatrix.h"
#include "SAM.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <numeric>
#include <algorithm>

using Catch::Matchers::WithinRel;


/*!
  \brief A simple SAM class for diagonal systems.
*/

class SAMdiag : public SAM
{
public:
  //! \brief The constructor initializes the arrays for a diagonal system.
  SAMdiag(int n)
  {
    nmmnpc = nel = nnod = ndof = neq = n;
    mmnpc  = new int[n];
    mpmnpc = new int[n+1];
    madof  = new int[n+1];
    msc    = new int[n];
    std::iota(mpmnpc,mpmnpc+n+1,1);
    std::iota(mmnpc ,mmnpc +n  ,1);
    std::iota(madof ,madof +n+1,1);
    std::fill(msc   ,msc   +n  ,1);
    REQUIRE(this->initSystemEquations());
  }

  //! \brief Empty destructor.
  virtual ~SAMdiag() {}
};


TEST_CASE("TestDiagMatrix.AssembleAndSolve")
{
  const int n = 6;

  SAMdiag sam(n);
  StdVector b(n), x;
  DiagMatrix A;

  A.initAssembly(sam,false);
  A.init();
  b.init();

  Matrix eM(1,1);
  for (int e = 1; e <= n; e++)
  {
    eM(1,1) = (double)e;
    REQUIRE(A.assemble(eM,sam,e));
    REQUIRE(sam.assembleSystem(b,{(double)2*e*e},e));
  }
  REQUIRE(A.solve(b,x));

  std::cout <<"A = "<< A;
  std::cout <<"b = "<< b;
  std::cout <<"x = "<< x;

  for (int i = 1; i <= n; i++)
    REQUIRE_THAT(x(i), WithinRel(static_cast<double>(2*i)));
}
