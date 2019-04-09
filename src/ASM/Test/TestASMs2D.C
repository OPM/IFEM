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


class ASMSquare : public ASMs2D
{
public:
  ASMSquare()
  {
    std::stringstream geo("200 1 0 0\n2 0\n2 2\n0 0 1 1\n2 2\n0 0 1 1\n0 0\n1 0\n0 1\n1 1\n");
    EXPECT_TRUE(this->read(geo));
  }
  virtual ~ASMSquare() {}
};


TEST(TestASMs2D, ElementNeighbours)
{
  ASMSquare pch1;
  ASSERT_TRUE(pch1.uniformRefine(0,1));
  ASSERT_TRUE(pch1.uniformRefine(1,1));
  ASSERT_TRUE(pch1.generateFEMTopology());
  ASMbase::NeighArray neigh(4);
  pch1.getNeighbours(neigh);
  const std::array<std::vector<int>,4> ref = {{{-1,  1, -1,  2},
                                               { 0, -1, -1,  3},
                                               {-1,  3,  0, -1},
                                               { 2, -1,  1, -1}}};
  ASSERT_EQ(neigh.size(), 4U);
  for (size_t n = 0; n < neigh.size(); ++n) {
    ASSERT_EQ(neigh[n].size(), ref[n].size());
    for (size_t i = 0; i < neigh[n].size(); ++i)
      EXPECT_EQ(neigh[n][i], ref[n][i]);
  }
}


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
