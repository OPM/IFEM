//==============================================================================
//!
//! \file TestASMs3D.C
//!
//! \date Feb 14 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for structured 3D spline FE models.
//!
//==============================================================================

#include "ASMs3D.h"
#include "SIM3D.h"

#include "gtest/gtest.h"
#include <numeric>


class ASMCube : public ASMs3D
{
public:
  ASMCube()
  {
    std::stringstream geo("700 1 0 0\n3 0\n2 2\n0 0 1 1\n2 2\n0 0 1 1\n2 2\n0 0 1 1\n0 0 0\n1 0 0\n0 1 0\n1 1 0\n0 0 1\n1 0 0\n0 1 1\n1 1 1\n");
    EXPECT_TRUE(this->read(geo));
  }
  virtual ~ASMCube() {}
};


TEST(TestASMs3D, ElementNeighbours)
{
  ASMCube pch1;
  ASSERT_TRUE(pch1.uniformRefine(0,1));
  ASSERT_TRUE(pch1.uniformRefine(1,1));
  ASSERT_TRUE(pch1.uniformRefine(2,1));
  ASSERT_TRUE(pch1.generateFEMTopology());
  ASMbase::NeighArray neigh(8);
  pch1.getNeighbours(neigh);
  const std::array<std::vector<int>,8> ref = {{{-1,  1, -1,  2, -1,  4},
                                               { 0, -1, -1,  3, -1,  5},
                                               {-1,  3,  0, -1, -1,  6},
                                               { 2, -1,  1, -1, -1,  7},
                                               {-1,  5, -1,  6,  0, -1},
                                               { 4, -1, -1,  7,  1, -1},
                                               {-1,  7,  4, -1,  2, -1},
                                               { 6, -1,  5, -1,  3, -1}}};
  ASSERT_EQ(neigh.size(), 8U);
  for (size_t n = 0; n < neigh.size(); ++n) {
    ASSERT_EQ(neigh[n].size(), ref[n].size());
    for (size_t i = 0; i < neigh[n].size(); ++i)
      EXPECT_EQ(neigh[n][i], ref[n][i]);
  }
}


class TestASMs3D : public testing::Test,
                   public testing::WithParamInterface<int>
{
};


TEST_P(TestASMs3D, Connect)
{
  if (GetParam() > 7)
    return;

  SIM3D sim(3);
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.createFEMmodel());
}


TEST_P(TestASMs3D, ConnectUneven)
{
  SIM3D sim(1);
  std::stringstream str;
  str << "src/ASM/Test/refdata/3d_uneven";
  if (GetParam() > 0)
    str << "_" << (GetParam()-1)/3 << (GetParam()-1) % 3;
  str << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.createFEMmodel());
}


const std::vector<int> orientations3D = {0,1,2,3,5,6,7,8,9,10,11,12};
INSTANTIATE_TEST_CASE_P(TestASMs3D, TestASMs3D, testing::ValuesIn(orientations3D));
