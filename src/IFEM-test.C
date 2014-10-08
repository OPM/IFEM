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

#include "IFEM.h"
#include "Profiler.h"

#include <cstdio>
#include <cstdlib>

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  IFEM::Init(argc, argv);
  Profiler prof(argv[0]);

  int ret = RUN_ALL_TESTS();

  return ret;
}
