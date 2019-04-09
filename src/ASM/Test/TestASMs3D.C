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
#include <array>

#include "gtest/gtest.h"


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


TEST(TestASMs3D, ElementConnectivities)
{
  ASMCube pch1;
  ASSERT_TRUE(pch1.uniformRefine(0,1));
  ASSERT_TRUE(pch1.uniformRefine(1,1));
  ASSERT_TRUE(pch1.uniformRefine(2,1));
  ASSERT_TRUE(pch1.generateFEMTopology());
  IntMat neigh(8);
  pch1.getElmConnectivities(neigh);
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


INSTANTIATE_TEST_CASE_P(TestASMs3D, TestASMs3D,
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
      ASMdegenerate3D pch(iface,iedge);
      ASSERT_TRUE(pch.uniformRefine(0,4));
      ASSERT_TRUE(pch.uniformRefine(1,1));
      ASSERT_TRUE(pch.uniformRefine(2,3));
      ASSERT_TRUE(pch.generateFEMTopology());
      std::cout <<"Degenerating F"<< iface <<" onto E"<< iedge << std::endl;
#ifdef SP_DEBUG
      pch.write(std::cout);
#endif
      const std::array<int,4>& face = faceTop[iface-1];
      if (iedge == 0 || std::find(face.begin(),face.end(),iedge) != face.end())
        EXPECT_TRUE(pch.collapseFace(iface,iedge));
      else
        EXPECT_FALSE(pch.collapseFace(iface,iedge));
    }
}
