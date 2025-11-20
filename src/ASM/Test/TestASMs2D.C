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
#include "ASMTrap.h"
#include "ASM2DTests.h"

namespace {

class ASMdegenerate2D : public ASMs2D
{
public:
  explicit ASMdegenerate2D(int iedge)
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
};

}


TEST_CASE("TestASMs2D.BoundaryElements")
{
  ASM2DTests<ASMSquare>::BoundaryElements();
}


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


TEST_CASE("TestASMs2D.Connect")
{
  ASM2DTests<ASMSquare>::Connect(ASM::Spline);
}


TEST_CASE("TestASMs2D.ConstrainEdge")
{
  ASM2DTests<ASMSquare>::ConstrainEdge();
}


TEST_CASE("TestASMs2D.ConstrainEdgeOpen")
{
  ASM2DTests<ASMSquare>::ConstrainEdgeOpen();
}


TEST_CASE("TestASMs2D.ElementConnectivities")
{
  const auto ref = std::array{
    std::array{-1,  1, -1,  2},
    std::array{ 0, -1, -1,  3},
    std::array{-1,  3,  0, -1},
    std::array{ 2, -1,  1, -1}
  };
  ASM2DTests<ASMSquare>::GetElementConnectivities(ref);
}


TEST_CASE("TestASMs2D.ElmNodes")
{
  const auto ref = std::array{
    std::array{0,1,3,4},
    std::array{1,2,4,5},
    std::array{3,4,6,7},
    std::array{4,5,7,8},
  };
  const auto ref_proj = std::array{
    std::array{0,1,2,4,5,6,8,9,10},
    std::array{1,2,3,5,6,7,9,10,11},
    std::array{4,5,6,8,9,10,12,13,14},
    std::array{5,6,7,9,10,11,13,14,15},
  };

  ASM2DTests<ASMSquare>::ElmNodes(ref, ref_proj);
}


TEST_CASE("TestASMs2D.GetElementCorners")
{
  ASM2DTests<ASMTrap<ASMs2D>>::GetElementCorners();
}


TEST_CASE("TestASMs2D.Write")
{
  ASM2DTests<ASMSquare>::Write(false);
}
