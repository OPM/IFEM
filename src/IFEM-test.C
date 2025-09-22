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

#include <algorithm>

#ifdef HAVE_MPI
#include <mpi.h>
#endif


/*!
  \brief Main program for the IFEM unit tests.
*/

int main (int argc, char** argv)
{
  const bool list_tests =
        std::any_of(argv+1, argv + argc,
                    [](const char* arg)
                    { return !strcmp(arg, "--gtest_list_tests"); });

  testing::InitGoogleTest(&argc, argv);
  IFEM::Init(argc, argv, nullptr, list_tests);
  Profiler prof(argv[0], false);

#ifdef HAVE_MPI
  if (list_tests) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0)
      return 0;
  }
#endif

  return RUN_ALL_TESTS();
}
