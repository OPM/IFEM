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


/*!
  \brief Main program for the IFEM unit tests.
*/

int main (int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  IFEM::Init(argc, argv);
  Profiler prof(argv[0], false);

  return RUN_ALL_TESTS();
}
