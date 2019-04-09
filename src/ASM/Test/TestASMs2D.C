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


TEST(TestASMs2D, ElementConnectivities)
{
  ASMSquare pch1;
  ASSERT_TRUE(pch1.uniformRefine(0,1));
  ASSERT_TRUE(pch1.uniformRefine(1,1));
  ASSERT_TRUE(pch1.generateFEMTopology());
  IntMat neigh(4);
  pch1.getElmConnectivities(neigh);
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


INSTANTIATE_TEST_CASE_P(TestASMs2D, TestASMs2D, testing::Values(0,1));


class ASMdegenerate2D : public ASMs2D
{
public:
  ASMdegenerate2D(int iedge)
  {
    // Create a degenerated patch where edge "iedge" is collapsed into a vertex
    std::stringstream geo; // -- Xi ----    --- Eta ----
    geo <<"200 1 0 0 2 0"<<" 2 2 0 0 1 1"<<" 2 2 0 0 1 1"<<" 0 0 ";
    switch (iedge) {
    case 1:
      geo <<"3 0 0 0 3 1\n"; // P3=P1
      break;
    case 2:
      geo <<"3 0 0 1 3 0\n"; // P4=P2
      break;
    case 3:
      geo <<"0 0 0 1 3 1\n"; // P2=P1
      break;
    case 4:
      geo <<"3 0 0 1 0 1\n"; // P4=P3
      break;
    default:
      geo <<"3 0 0 1 3 1\n";
      break;
    }
    EXPECT_TRUE(this->read(geo));
  }
  virtual ~ASMdegenerate2D() {}
};


TEST(TestASMs2D, Collapse)
{
  for (int iedge = 1; iedge <= 4; iedge++)
  {
    ASMdegenerate2D pch(iedge);
    ASSERT_TRUE(pch.uniformRefine(0,2));
    ASSERT_TRUE(pch.uniformRefine(1,1));
    ASSERT_TRUE(pch.generateFEMTopology());
    std::cout <<"Degenerating E"<< iedge << std::endl;
#ifdef SP_DEBUG
    pch.write(std::cout);
#endif
    EXPECT_TRUE(pch.collapseEdge(iedge));
  }
}
