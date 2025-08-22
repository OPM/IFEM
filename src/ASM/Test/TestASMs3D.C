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

#include "ASMCube.h"
#include "IntegrandBase.h"
#include "GlobalIntegral.h"
#include "SIM3D.h"
#include <array>

#include "gtest/gtest.h"


TEST(TestASMs3D, ElementConnectivities)
{
  ASMCube pch1;
  ASMbase::resetNumbering();
  ASSERT_TRUE(pch1.uniformRefine(0,1));
  ASSERT_TRUE(pch1.uniformRefine(1,1));
  ASSERT_TRUE(pch1.uniformRefine(2,1));
  ASSERT_TRUE(pch1.generateFEMTopology());
  const size_t nel = pch1.getNoElms();
  pch1.shiftElemNumbers(nel);
  IntMat neighGlb(2*nel), neighLoc(nel);
  pch1.getElmConnectivities(neighGlb);
  pch1.getElmConnectivities(neighLoc, ASM::GEOMETRY_BASIS);
  const std::array<std::vector<int>,8> ref = {{{-1,  1, -1,  2, -1,  4},
                                               { 0, -1, -1,  3, -1,  5},
                                               {-1,  3,  0, -1, -1,  6},
                                               { 2, -1,  1, -1, -1,  7},
                                               {-1,  5, -1,  6,  0, -1},
                                               { 4, -1, -1,  7,  1, -1},
                                               {-1,  7,  4, -1,  2, -1},
                                               { 6, -1,  5, -1,  3, -1}}};
  ASSERT_EQ(neighGlb.size(), 2*nel);
  ASSERT_EQ(neighLoc.size(), nel);
  for (size_t n = 0; n < neighLoc.size(); ++n) {
    ASSERT_EQ(neighLoc[n].size(), ref[n].size());
    ASSERT_EQ(neighGlb[n+nel].size(), ref[n].size());
    for (size_t i = 0; i < neighLoc[n].size(); ++i) {
      EXPECT_EQ(neighLoc[n][i], ref[n][i]);
      EXPECT_EQ(neighGlb[n+nel][i], ref[n][i] > -1 ? ref[n][i] + nel : -1);
    }
  }
}


TEST(TestASMs3D, BoundaryElements)
{
  ASMCube pch1;
  ASMbase::resetNumbering();
  ASSERT_TRUE(pch1.uniformRefine(0,1));
  ASSERT_TRUE(pch1.uniformRefine(1,1));
  ASSERT_TRUE(pch1.uniformRefine(2,1));
  ASSERT_TRUE(pch1.generateFEMTopology());

  const std::array<std::array<int,4>,6> ref = {{{{0, 2, 4, 6}},
                                                {{1, 3, 5, 7}},
                                                {{0, 1, 4, 5}},
                                                {{2, 3, 6, 7}},
                                                {{0, 1, 2, 3}},
                                                {{4, 5, 6, 7}}}};

  std::array<IntVec,6> n;
  for (size_t i = 1; i <= 6; ++i) {
    pch1.getBoundaryElms(i, n[i-1]);
    ASSERT_EQ(n[i-1].size(), ref[i-1].size());
    for (size_t j = 0; j < ref[i-1].size(); ++j)
      EXPECT_EQ(n[i-1][j], ref[i-1][j]);
  }
}


TEST(TestASMs3D, Write)
{
  ASMbase::resetNumbering();
  ASMCube pch1;
  EXPECT_TRUE(pch1.generateFEMTopology());

  std::stringstream str;
  EXPECT_TRUE(pch1.write(str, 1));
  EXPECT_EQ(str.str(), ASMCube::cube);

  EXPECT_FALSE(pch1.write(str, 2));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::GEOMETRY_BASIS));
  EXPECT_EQ(str.str(), ASMCube::cube);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::PROJECTION_BASIS));
  EXPECT_EQ(str.str(), ASMCube::cube);

  EXPECT_FALSE(pch1.write(str, ASM::PROJECTION_BASIS_2));
  EXPECT_FALSE(pch1.write(str, ASM::REFINEMENT_BASIS));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::INTEGRATION_BASIS));
  EXPECT_EQ(str.str(), ASMCube::cube);
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


INSTANTIATE_TEST_SUITE_P(TestASMs3D, TestASMs3D,
                         testing::Values(0,1,2,3,5,6,7,8,9,10,11,12));


class ASMdegenerate3D : public ASMs3D
{
public:
  ASMdegenerate3D(int iface, int iedge)
  {
    std::stringstream geo; // --- Xi ---    --- Eta ----    --- Zeta ---
    geo <<"700 1 0 0 3 0"<<" 2 2 0 0 1 1"<<" 2 2 0 0 1 1"<<" 2 2 0 0 1 1\n";
    auto&& defaultGeo = [&geo]() {
      geo <<"0 0 0  3 0 0  0 1 0  3 1 0\n"
          <<"0 0 2  3 0 2  0 1 2  3 1 2\n";
    };

    // Create a degenerated patch where face "iface"
    // is collapsed into the edge "iedge"
    switch (iface) {
    case 1:
      switch (iedge) {
      case 0:
        geo <<"0 0 0  3 0 0  0 0 0  3 1 0\n"
            <<"0 0 0  3 0 2  0 0 0  3 1 2\n";
        break;
      case 5:
        geo <<"0 0 0  3 0 0  0 1 0  3 1 0\n"
            <<"0 0 0  3 0 2  0 1 0  3 1 2\n";
        break;
      case 7:
        geo <<"0 0 2  3 0 0  0 1 2  3 1 0\n"
            <<"0 0 2  3 0 2  0 1 2  3 1 2\n";
        break;
      case 9:
        geo <<"0 0 0  3 0 0  0 0 0  3 1 0\n"
            <<"0 0 2  3 0 2  0 0 2  3 1 2\n";
        break;
      case 11:
        geo <<"0 1 0  3 0 0  0 1 0  3 1 0\n"
            <<"0 1 2  3 0 2  0 1 2  3 1 0\n";
        break;
      default:
        defaultGeo();
      }
      break;

    case 2:
      switch (iedge) {
      case 0:
        geo <<"0 0 0  3 0 0  0 1 0  3 0 0\n"
            <<"0 0 2  3 0 0  0 1 2  3 0 0\n";
        break;
      case 6:
        geo <<"0 0 0  3 0 0  0 1 0  3 1 0\n"
            <<"0 0 2  3 0 0  0 1 2  3 1 0\n";
        break;
      case 8:
        geo <<"0 0 0  3 0 2  0 1 0  3 1 2\n"
            <<"0 0 2  3 0 2  0 1 2  3 1 2\n";
        break;
      case 10:
        geo <<"0 0 0  3 0 0  0 1 0  3 0 0\n"
            <<"0 0 2  3 0 2  0 1 2  3 0 2\n";
        break;
      case 12:
        geo <<"0 0 0  3 1 0  0 1 0  3 1 0\n"
            <<"0 0 2  3 1 2  0 1 2  3 1 2\n";
        break;
      default:
        defaultGeo();
      }
      break;

    case 3:
      switch (iedge) {
      case 0:
        geo <<"0 0 0  0 0 0  0 1 0  3 1 0\n"
            <<"0 0 0  0 0 0  0 1 2  3 1 2\n";
        break;
      case 1:
        geo <<"0 0 0  3 0 0  0 1 0  3 1 0\n"
            <<"0 0 0  3 0 0  0 1 2  3 1 2\n";
        break;
      case 3:
        geo <<"0 0 2  3 0 2  0 1 0  3 1 0\n"
            <<"0 0 2  3 0 2  0 1 2  3 1 2\n";
        break;
      case 9:
        geo <<"0 0 0  0 0 0  0 1 0  3 1 0\n"
            <<"0 0 2  0 0 2  0 1 2  3 1 2\n";
        break;
      case 10:
        geo <<"3 0 0  3 0 0  0 1 0  3 1 0\n"
            <<"3 0 2  3 0 2  0 1 2  3 1 2\n";
        break;
      default:
        defaultGeo();
      }
      break;

    case 4:
      switch (iedge) {
      case 0:
        geo <<"0 0 0  3 0 0  0 1 0  0 1 0\n"
            <<"0 0 2  3 0 2  0 1 0  0 1 0\n";
        break;
      case 2:
        geo <<"0 0 0  3 0 0  0 1 0  3 1 0\n"
            <<"0 0 2  3 0 2  0 1 0  3 1 0\n";
        break;
      case 4:
        geo <<"0 0 0  3 0 0  0 1 2  3 1 2\n"
            <<"0 0 2  3 0 2  0 1 2  3 1 2\n";
        break;
      case 11:
        geo <<"0 0 0  3 0 0  0 1 0  0 1 0\n"
            <<"0 0 2  3 0 2  0 1 2  0 1 2\n";
        break;
      case 12:
        geo <<"0 0 0  3 0 0  3 1 0  3 1 0\n"
            <<"0 0 2  3 0 2  3 1 2  3 1 2\n";
        break;
      default:
        defaultGeo();
      }
      break;

    case 5:
      switch (iedge) {
      case 0:
        geo <<"0 0 0  0 0 0  0 0 0  0 0 0\n"
            <<"0 0 2  3 0 2  0 1 2  3 1 2\n";
        break;
      case 1:
        geo <<"0 0 0  3 0 0  0 0 0  3 0 0\n"
            <<"0 0 2  3 0 2  0 1 2  3 1 2\n";
        break;
      case 2:
        geo <<"0 1 0  3 1 0  0 1 0  3 1 0\n"
            <<"0 0 2  3 0 2  0 1 2  3 1 2\n";
        break;
      case 5:
        geo <<"0 0 0  0 0 0  0 1 0  0 1 0\n"
            <<"0 0 2  3 0 2  0 1 2  3 1 2\n";
        break;
      case 6:
        geo <<"3 0 0  3 0 0  3 1 0  3 1 0\n"
            <<"0 0 2  3 0 2  0 1 2  3 1 2\n";
        break;
      default:
        defaultGeo();
      }
      break;

    case 6:
      switch (iedge) {
      case 0:
        geo <<"0 0 0  3 0 0  0 1 0  3 1 0\n"
            <<"1 0 2  1 0 2  1 0 2  1 0 2\n";
        break;
      case 3:
        geo <<"0 0 0  3 0 0  0 1 0  3 1 0\n"
            <<"0 0 2  3 0 2  0 0 2  3 0 2\n";
        break;
      case 4:
        geo <<"0 0 0  3 0 0  0 1 0  3 1 0\n"
            <<"0 1 2  3 1 2  0 1 2  3 1 2\n";
        break;
      case 7:
        geo <<"0 0 0  3 0 0  0 1 0  3 1 0\n"
            <<"0 0 2  0 0 2  0 1 2  0 1 2\n";
        break;
      case 8:
        geo <<"0 0 0  3 0 0  0 1 0  3 1 0\n"
            <<"3 0 2  3 0 2  3 1 2  3 1 2\n";
        break;
      default:
        defaultGeo();
      }
      break;

    default:
      defaultGeo();
      break;
    }
    EXPECT_TRUE(this->read(geo));
  }
  virtual ~ASMdegenerate3D() {}
};


TEST(TestASMs3D, Collapse)
{
  // Face-to-edge topology, see
  // https://github.com/OPM/IFEM/blob/master/doc/sim-input.pdf Figures 2 and 3.
  std::array<std::array<int,4>,6> faceTop = {{
      {{ 5, 7, 9,11 }},
      {{ 6, 8,10,12 }},
      {{ 1, 3, 9,10 }},
      {{ 2, 4,11,12 }},
      {{ 1, 2, 5, 6 }},
      {{ 3, 4, 7, 8 }} }};

  for (int iface = 1; iface <= 6; iface++)
    for (int iedge = 0; iedge <= 12; iedge++)
    {
      ASMbase::resetNumbering();
      ASMdegenerate3D pch(iface,iedge);
      ASSERT_TRUE(pch.uniformRefine(0,4));
      ASSERT_TRUE(pch.uniformRefine(1,1));
      ASSERT_TRUE(pch.uniformRefine(2,3));
      ASSERT_TRUE(pch.generateFEMTopology());
      std::cout <<"Degenerating F"<< iface <<" onto E"<< iedge << std::endl;
#ifdef SP_DEBUG
      pch.write(std::cout,0);
#endif
      const std::array<int,4>& face = faceTop[iface-1];
      if (iedge == 0 || std::find(face.begin(),face.end(),iedge) != face.end())
        EXPECT_TRUE(pch.collapseFace(iface,iedge));
      else
        EXPECT_FALSE(pch.collapseFace(iface,iedge));
    }
}

class NoProblem : public IntegrandBase
{
public:
  NoProblem() : IntegrandBase(3) {}
  virtual ~NoProblem() {}
protected:
  virtual bool evalBou(LocalIntegral&, const FiniteElement&,
                       const Vec3&, const Vec3&) const { return true; }
};


TEST(TestASMs3D, FaceIntegrate)
{
  GlobalIntegral dummy;
  NoProblem prb;
  ASMCube patch;
  ASMbase::resetNumbering();

  ASSERT_TRUE(patch.raiseOrder(2,1,0));
  ASSERT_TRUE(patch.generateFEMTopology());
  for (int lIndex = 1; lIndex <= 6; lIndex++)
  {
    patch.generateThreadGroups(lIndex,false,false);
    ASSERT_TRUE(patch.integrate(prb,lIndex,dummy,TimeDomain()));
  }

  patch.clear(true);
  ASMbase::resetNumbering();

  ASSERT_TRUE(patch.uniformRefine(0,3));
  ASSERT_TRUE(patch.uniformRefine(1,2));
  ASSERT_TRUE(patch.uniformRefine(2,1));
  ASSERT_TRUE(patch.generateFEMTopology());
  for (int lIndex = 1; lIndex <= 6; lIndex++)
  {
    patch.generateThreadGroups(lIndex,false,false);
    ASSERT_TRUE(patch.integrate(prb,lIndex,dummy,TimeDomain()));
  }
}


TEST(TestASMs3D, ElmNodes)
{
  ASMbase::resetNumbering();

  ASMCube pch1;
  pch1.createProjectionBasis(true);
  pch1.raiseOrder(1,1,1);
  pch1.uniformRefine(0,1);
  pch1.uniformRefine(1,1);
  pch1.uniformRefine(2,1);
  pch1.createProjectionBasis(false);
  ASSERT_TRUE(pch1.uniformRefine(0,1));
  ASSERT_TRUE(pch1.uniformRefine(1,1));
  ASSERT_TRUE(pch1.uniformRefine(2,1));
  ASSERT_TRUE(pch1.generateFEMTopology());

  const IntMat mnpc = pch1.getElmNodes(1);

  const auto ref = std::array{
      std::array{0,1,3,4,9,10,12,13},
      std::array{1,2,4,5,10,11,13,14},
      std::array{3,4,6,7,12,13,15,16},
      std::array{4,5,7,8,13,14,16,17},
      std::array{9,10,12,13,18,19,21,22},
      std::array{10,11,13,14,19,20,22,23},
      std::array{12,13,15,16,21,22,24,25},
      std::array{13,14,16,17,22,23,25,26},
  };
  ASSERT_EQ(mnpc.size(), ref.size());
  for (size_t i = 0; i < mnpc.size(); ++i) {
    EXPECT_EQ(mnpc[i].size(), ref[i].size());
    for (size_t j = 0; j < mnpc[i].size(); ++j)
      EXPECT_EQ(mnpc[i][j], ref[i][j]);
  }

  const auto ref_proj = std::array{
      std::array{0,1,2,4,5,6,8,9,10,16,17,18,20,21,22,24,25,26,32,33,34,36,37,38,40,41,42},
      std::array{1,2,3,5,6,7,9,10,11,17,18,19,21,22,23,25,26,27,33,34,35,37,38,39,41,42,43},
      std::array{4,5,6,8,9,10,12,13,14,20,21,22,24,25,26,28,29,30,36,37,38,40,41,42,44,45,46},
      std::array{5,6,7,9,10,11,13,14,15,21,22,23,25,26,27,29,30,31,37,38,39,41,42,43,45,46,47},
      std::array{16,17,18,20,21,22,24,25,26,32,33,34,36,37,38,40,41,42,48,49,50,52,53,54,56,57,58},
      std::array{17,18,19,21,22,23,25,26,27,33,34,35,37,38,39,41,42,43,49,50,51,53,54,55,57,58,59},
      std::array{20,21,22,24,25,26,28,29,30,36,37,38,40,41,42,44,45,46,52,53,54,56,57,58,60,61,62},
      std::array{21,22,23,25,26,27,29,30,31,37,38,39,41,42,43,45,46,47,53,54,55,57,58,59,61,62,63},
  };
  const IntMat mnpc_proj = pch1.getElmNodes(ASM::PROJECTION_BASIS);
  ASSERT_EQ(mnpc_proj.size(), ref_proj.size());
  for (size_t i = 0; i < mnpc_proj.size(); ++i) {
    EXPECT_EQ(mnpc_proj[i].size(), ref_proj[i].size());
    for (size_t j = 0; j < mnpc_proj[i].size(); ++j)
      EXPECT_EQ(mnpc_proj[i][j], ref_proj[i][j]);
  }
}
