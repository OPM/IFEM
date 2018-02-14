//==============================================================================
//!
//! \file TestASMs2D.C
//!
//! \date Feb 14 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for structured 2D spline FE models.
//!
//==============================================================================

#include "ASMs2D.h"
#include "SIM2D.h"

#include "gtest/gtest.h"
#include <numeric>


class TestASMs2D : public testing::Test,
                   public testing::WithParamInterface<int>
{
};


TEST_P(TestASMs2D, Connect)
{
  SIM2D sim(1);
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.createFEMmodel());
}

const std::vector<int> orientations2D = {0,1};
INSTANTIATE_TEST_CASE_P(TestASMs2D, TestASMs2D, testing::ValuesIn(orientations2D));
