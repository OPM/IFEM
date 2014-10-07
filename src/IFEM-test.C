//==============================================================================
//!
//! \file IFEM-test.C
//!
//! \date Oct 07 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Main testing application function.
//!
//==============================================================================

#include "gtest/gtest.h"

#include <cstdio>
#include <cstdlib>

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);

  int ret = RUN_ALL_TESTS();

  return ret;
}
