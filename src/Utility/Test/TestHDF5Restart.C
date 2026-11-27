// $Id$
//==============================================================================
//!
//! \file TestHDF5Restart.C
//!
//! \date Oct 9 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for the restart output stride.
//!
//==============================================================================

#include "HDF5Restart.h"
#include "ProcessAdm.h"
#include "TimeStep.h"

#include "Catch2Support.h"

#include <vector>


namespace {

//! \brief Returns the steps dumped from step \a first to \a last.
std::vector<int> dumped (HDF5Restart& hdf, int first, int last)
{
  TimeStep tp;
  std::vector<int> steps;
  for (tp.step = first; tp.step <= last; ++tp.step)
    if (hdf.dumpStep(tp))
      steps.push_back(tp.step);
  return steps;
}

}


TEST_CASE("TestHDF5Restart.Stride")
{
  ProcessAdm adm;
  HDF5Restart hdf("dummy", adm, 50);
  REQUIRE(dumped(hdf,0,160) == std::vector<int>{50, 100, 150});
}


TEST_CASE("TestHDF5Restart.StrideAfterRestart")
{
  // A simulation restarted from step 201 continues at the multiples
  ProcessAdm adm;
  HDF5Restart hdf("dummy", adm, 50, 2, 201);
  REQUIRE(dumped(hdf,201,310) == std::vector<int>{250, 300});

  // and the restored step is not dumped again
  HDF5Restart hdf2("dummy", adm, 50, 2, 200);
  REQUIRE(dumped(hdf2,200,260) == std::vector<int>{250});
}
