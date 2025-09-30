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

#include "Catch2Support.h"


TEST_CASE("TestASMs2D.ElementConnectivities")
{
  ASMSquare pch1;
  ASMbase::resetNumbering();
  REQUIRE(pch1.uniformRefine(0,1));
  REQUIRE(pch1.uniformRefine(1,1));
  REQUIRE(pch1.generateFEMTopology());
  const size_t nel = pch1.getNoElms();
  pch1.shiftElemNumbers(nel);
  IntMat neighGlb(2*nel), neighLoc(nel);
  pch1.getElmConnectivities(neighGlb);
  pch1.getElmConnectivities(neighLoc, ASM::GEOMETRY_BASIS);
  const std::array<std::vector<int>,4> ref = {{{-1,  1, -1,  2},
                                               { 0, -1, -1,  3},
                                               {-1,  3,  0, -1},
                                               { 2, -1,  1, -1}}};
  REQUIRE(neighLoc.size() == nel);
  REQUIRE(neighGlb.size() == 2*nel);
  for (size_t n = 0; n < neighLoc.size(); ++n) {
    REQUIRE(neighLoc[n].size() == ref[n].size());
    REQUIRE(neighGlb[n+nel].size() == ref[n].size());
    for (size_t i = 0; i < neighLoc[n].size(); ++i) {
      REQUIRE(neighLoc[n][i] == ref[n][i]);
      REQUIRE(neighGlb[n+nel][i] == (ref[n][i] > -1 ?
                                      ref[n][i] + static_cast<int>(nel) : -1));
    }
  }
}


TEST_CASE("TestASMs2D.BoundaryElements")
{
  ASMSquare pch1;
  ASMbase::resetNumbering();
  REQUIRE(pch1.uniformRefine(0,1));
  REQUIRE(pch1.uniformRefine(1,1));
  REQUIRE(pch1.generateFEMTopology());

  const std::array<std::array<int,2>,4> ref = {{{{0, 2}}, {{1, 3}}, {{0, 1}}, {{2, 3}}}};

  std::array<IntVec,4> n;
  for (size_t i = 1; i <= 4; ++i) {
    pch1.getBoundaryElms(i, n[i-1]);
    REQUIRE(n[i-1].size() == ref[i-1].size());
    for (size_t j = 0; j < ref[i-1].size(); ++j)
      REQUIRE(n[i-1][j] == ref[i-1][j]);
  }
}


TEST_CASE("TestASMs2D.Write")
{
  ASMbase::resetNumbering();
  ASMSquare pch1;
  REQUIRE(pch1.generateFEMTopology());

  std::stringstream str;
  REQUIRE(pch1.write(str, 1));
  REQUIRE(str.str() == ASMSquare::square);

  REQUIRE(!pch1.write(str, 2));

  str.str("");
  REQUIRE(pch1.write(str, ASM::GEOMETRY_BASIS));
  REQUIRE(str.str() == ASMSquare::square);

  str.str("");
  REQUIRE(pch1.write(str, ASM::PROJECTION_BASIS));
  REQUIRE(str.str() == ASMSquare::square);

  REQUIRE(!pch1.write(str, ASM::PROJECTION_BASIS_2));
  REQUIRE(!pch1.write(str, ASM::REFINEMENT_BASIS));

  str.str("");
  REQUIRE(pch1.write(str, ASM::INTEGRATION_BASIS));
  REQUIRE(str.str() == ASMSquare::square);
}


TEST_CASE("TestASMs2D.Connect")
{
  const int param = GENERATE(0,1);
  SECTION("Orient " + std::to_string(param)) {
    SIM2D sim(1);
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
    str << param << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.createFEMmodel());
  }
}


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
    REQUIRE(this->read(geo));
  }
  virtual ~ASMdegenerate2D() {}
};


TEST_CASE("TestASMs2D.Collapse")
{
  for (int iedge = 1; iedge <= 4; iedge++)
  {
    ASMbase::resetNumbering();
    ASMdegenerate2D pch(iedge);
    REQUIRE(pch.uniformRefine(0,2));
    REQUIRE(pch.uniformRefine(1,1));
    REQUIRE(pch.generateFEMTopology());
    std::cout <<"Degenerating E"<< iedge << std::endl;
#ifdef SP_DEBUG
    pch.write(std::cout,0);
#endif
    REQUIRE(pch.collapseEdge(iedge));
  }
}


TEST_CASE("TestASMs2D.ElmNodes")
{
  ASMbase::resetNumbering();

  ASMSquare pch1;
  pch1.createProjectionBasis(true);
  pch1.raiseOrder(1,1);
  pch1.uniformRefine(0,1);
  pch1.uniformRefine(1,1);
  pch1.createProjectionBasis(false);
  REQUIRE(pch1.uniformRefine(0,1));
  REQUIRE(pch1.uniformRefine(1,1));
  REQUIRE(pch1.generateFEMTopology());

  const IntMat mnpc = pch1.getElmNodes(1);

  const auto ref = std::array{
      std::array{0,1,3,4},
      std::array{1,2,4,5},
      std::array{3,4,6,7},
      std::array{4,5,7,8},
  };
  REQUIRE(mnpc.size() == ref.size());
  for (size_t i = 0; i < mnpc.size(); ++i) {
    REQUIRE(mnpc[i].size() == ref[i].size());
    for (size_t j = 0; j < mnpc[i].size(); ++j)
      REQUIRE(mnpc[i][j] == ref[i][j]);
  }

  const auto ref_proj = std::array{
      std::array{0,1,2,4,5,6,8,9,10},
      std::array{1,2,3,5,6,7,9,10,11},
      std::array{4,5,6,8,9,10,12,13,14},
      std::array{5,6,7,9,10,11,13,14,15},
  };
  const IntMat mnpc_proj = pch1.getElmNodes(ASM::PROJECTION_BASIS);
  REQUIRE(mnpc_proj.size() == ref_proj.size());
  for (size_t i = 0; i < mnpc_proj.size(); ++i) {
    REQUIRE(mnpc_proj[i].size() == ref_proj[i].size());
    for (size_t j = 0; j < mnpc_proj[i].size(); ++j)
      REQUIRE(mnpc_proj[i][j] == ref_proj[i][j]);
  }
}
