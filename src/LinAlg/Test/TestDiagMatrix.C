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

#include "gtest/gtest.h"
#include <numeric>
#include <algorithm>


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
    EXPECT_TRUE(this->initSystemEquations());
  }

  //! \brief Empty destructor.
  virtual ~SAMdiag() {}
};


TEST(TestDiagMatrix, AssembleAndSolve)
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
    EXPECT_TRUE(A.assemble(eM,sam,e));
    EXPECT_TRUE(sam.assembleSystem(b,{(double)2*e*e},e));
  }
  ASSERT_TRUE(A.solve(b,x));

  std::cout <<"A = "<< A;
  std::cout <<"b = "<< b;
  std::cout <<"x = "<< x;

  for (int i = 1; i <= n; i++)
    EXPECT_FLOAT_EQ(x(i),(double)2*i);
}
