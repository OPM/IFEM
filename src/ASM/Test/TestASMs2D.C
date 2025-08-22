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

#include "ASMSquare.h"
#include "SIM2D.h"

#include "gtest/gtest.h"


TEST(TestASMs2D, ElementConnectivities)
{
  ASMSquare pch1;
  ASMbase::resetNumbering();
  ASSERT_TRUE(pch1.uniformRefine(0,1));
  ASSERT_TRUE(pch1.uniformRefine(1,1));
  ASSERT_TRUE(pch1.generateFEMTopology());
  const size_t nel = pch1.getNoElms();
  pch1.shiftElemNumbers(nel);
  IntMat neighGlb(2*nel), neighLoc(nel);
  pch1.getElmConnectivities(neighGlb);
  pch1.getElmConnectivities(neighLoc, ASM::GEOMETRY_BASIS);
  const std::array<std::vector<int>,4> ref = {{{-1,  1, -1,  2},
                                               { 0, -1, -1,  3},
                                               {-1,  3,  0, -1},
                                               { 2, -1,  1, -1}}};
  ASSERT_EQ(neighLoc.size(), nel);
  ASSERT_EQ(neighGlb.size(), 2*nel);
  for (size_t n = 0; n < neighLoc.size(); ++n) {
    ASSERT_EQ(neighLoc[n].size(), ref[n].size());
    ASSERT_EQ(neighGlb[n+nel].size(), ref[n].size());
    for (size_t i = 0; i < neighLoc[n].size(); ++i) {
      EXPECT_EQ(neighLoc[n][i], ref[n][i]);
      EXPECT_EQ(neighGlb[n+nel][i], ref[n][i] > -1 ? ref[n][i] + nel : -1);
    }
  }
}


TEST(TestASMs2D, BoundaryElements)
{
  ASMSquare pch1;
  ASMbase::resetNumbering();
  ASSERT_TRUE(pch1.uniformRefine(0,1));
  ASSERT_TRUE(pch1.uniformRefine(1,1));
  ASSERT_TRUE(pch1.generateFEMTopology());

  const std::array<std::array<int,2>,4> ref = {{{{0, 2}}, {{1, 3}}, {{0, 1}}, {{2, 3}}}};

  std::array<IntVec,4> n;
  for (size_t i = 1; i <= 4; ++i) {
    pch1.getBoundaryElms(i, n[i-1]);
    ASSERT_EQ(n[i-1].size(), ref[i-1].size());
    for (size_t j = 0; j < ref[i-1].size(); ++j)
      EXPECT_EQ(n[i-1][j], ref[i-1][j]);
  }
}


TEST(TestASMs2D, Write)
{
  ASMbase::resetNumbering();
  ASMSquare pch1;
  EXPECT_TRUE(pch1.generateFEMTopology());

  std::stringstream str;
  EXPECT_TRUE(pch1.write(str, 1));
  EXPECT_EQ(str.str(), ASMSquare::square);

  EXPECT_FALSE(pch1.write(str, 2));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::GEOMETRY_BASIS));
  EXPECT_EQ(str.str(), ASMSquare::square);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::PROJECTION_BASIS));
  EXPECT_EQ(str.str(), ASMSquare::square);

  EXPECT_FALSE(pch1.write(str, ASM::PROJECTION_BASIS_2));
  EXPECT_FALSE(pch1.write(str, ASM::REFINEMENT_BASIS));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::INTEGRATION_BASIS));
  EXPECT_EQ(str.str(), ASMSquare::square);
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


INSTANTIATE_TEST_SUITE_P(TestASMs2D, TestASMs2D, testing::Values(0,1));


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
    ASMbase::resetNumbering();
    ASMdegenerate2D pch(iedge);
    ASSERT_TRUE(pch.uniformRefine(0,2));
    ASSERT_TRUE(pch.uniformRefine(1,1));
    ASSERT_TRUE(pch.generateFEMTopology());
    std::cout <<"Degenerating E"<< iedge << std::endl;
#ifdef SP_DEBUG
    pch.write(std::cout,0);
#endif
    EXPECT_TRUE(pch.collapseEdge(iedge));
  }
}


TEST(TestASMs2D, ElmNodes)
{
  ASMbase::resetNumbering();

  ASMSquare pch1;
  pch1.createProjectionBasis(true);
  pch1.raiseOrder(1,1);
  pch1.uniformRefine(0,1);
  pch1.uniformRefine(1,1);
  pch1.createProjectionBasis(false);
  ASSERT_TRUE(pch1.uniformRefine(0,1));
  ASSERT_TRUE(pch1.uniformRefine(1,1));
  ASSERT_TRUE(pch1.generateFEMTopology());

  const IntMat mnpc = pch1.getElmNodes(1);

  const auto ref = std::array{
      std::array{0,1,3,4},
      std::array{1,2,4,5},
      std::array{3,4,6,7},
      std::array{4,5,7,8},
  };
  ASSERT_EQ(mnpc.size(), ref.size());
  for (size_t i = 0; i < mnpc.size(); ++i) {
    EXPECT_EQ(mnpc[i].size(), ref[i].size());
    for (size_t j = 0; j < mnpc[i].size(); ++j)
      EXPECT_EQ(mnpc[i][j], ref[i][j]);
  }

  const auto ref_proj = std::array{
      std::array{0,1,2,4,5,6,8,9,10},
      std::array{1,2,3,5,6,7,9,10,11},
      std::array{4,5,6,8,9,10,12,13,14},
      std::array{5,6,7,9,10,11,13,14,15},
  };
  const IntMat mnpc_proj = pch1.getElmNodes(ASM::PROJECTION_BASIS);
  ASSERT_EQ(mnpc_proj.size(), ref_proj.size());
  for (size_t i = 0; i < mnpc_proj.size(); ++i) {
    EXPECT_EQ(mnpc_proj[i].size(), ref_proj[i].size());
    for (size_t j = 0; j < mnpc_proj[i].size(); ++j)
      EXPECT_EQ(mnpc_proj[i][j], ref_proj[i][j]);
  }
}
