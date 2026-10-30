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

#include "ASM3DTests.h"
#include "Functions.h"
#include "MPC.h"
#include "ASMCube.h"
#include "ASMTrap3.h"

#include <array>
#include <functional>

namespace {

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
    REQUIRE(this->read(geo));
  }
};

}


TEST_CASE("TestASMs3D.BoundaryNodes")
{
  using Func = std::function<int(int,int)>;
  const auto ref = std::array{
    Func{[](int i, int j) { return 1+4*(4*j+i); }},
    Func{[](int i, int j) { return 4+4*(4*j+i); }},
    Func{[](int i, int j) { return 16*j+i+1; }},
    Func{[](int i, int j) { return 16*j+i+13; }},
    Func{[](int i, int j) { return 4*j+i+1; }},
    Func{[](int i, int j) { return 16*3+4*j+i+1; }},
  };

  ASM3DTests<ASMCube>::BoundaryNodes(ref);
}


TEST_CASE("TestASMs3D.BoundaryElements")
{
  ASM3DTests<ASMCube>::BoundaryElements();
}


TEST_CASE("TestASMs3D.Collapse")
{
  // Face-to-edge topology, see
  // https://github.com/OPM/IFEM/blob/master/doc/sim-input.pdf Figures 2 and 3.
  const auto faceTop = std::array{
      std::array{5, 7,  9, 11},
      std::array{6, 8, 10, 12},
      std::array{1, 3,  9, 10},
      std::array{2, 4, 11, 12},
      std::array{1, 2,  5,  6},
      std::array{3, 4,  7,  8},
  };

  for (int iface = 1; iface <= 6; iface++)
    for (int iedge = 0; iedge <= 12; iedge++)
    {
      ASMbase::resetNumbering();
      ASMdegenerate3D pch(iface,iedge);
      REQUIRE(pch.uniformRefine(0,4));
      REQUIRE(pch.uniformRefine(1,1));
      REQUIRE(pch.uniformRefine(2,3));
      REQUIRE(pch.generateFEMTopology());
      std::cout <<"Degenerating F"<< iface <<" onto E"<< iedge << std::endl;
#ifdef SP_DEBUG
      pch.write(std::cout,0);
#endif
      const std::array<int,4>& face = faceTop[iface-1];
      if (iedge == 0 || std::find(face.begin(),face.end(),iedge) != face.end())
        REQUIRE(pch.collapseFace(iface,iedge));
      else
        REQUIRE(!pch.collapseFace(iface,iedge));
    }
}


TEST_CASE("TestASMs3D.Connect")
{
  ASM3DTests<ASMCube>::Connect(ASM::Spline);
}


TEST_CASE("TestASMs3D.ConnectUneven")
{
  ASM3DTests<ASMCube>::ConnectUneven(ASM::Spline);
}


TEST_CASE("TestASMs3D.ConstrainFace")
{
  ASM3DTests<ASMCube>::ConstrainFace();
}


TEST_CASE("TestASMs3D.ElementConnectivities")
{
  const auto ref = std::array{
    std::array{-1,  1, -1,  2, -1,  4},
    std::array{ 0, -1, -1,  3, -1,  5},
    std::array{-1,  3,  0, -1, -1,  6},
    std::array{ 2, -1,  1, -1, -1,  7},
    std::array{-1,  5, -1,  6,  0, -1},
    std::array{ 4, -1, -1,  7,  1, -1},
    std::array{-1,  7,  4, -1,  2, -1},
    std::array{ 6, -1,  5, -1,  3, -1},
  };

  ASM3DTests<ASMCube>::ElementConnectivities(ref);
}


TEST_CASE("TestASMs3D.ElmNodes")
{
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

  ASM3DTests<ASMCube>::ElmNodes(ref, ref_proj);
}


TEST_CASE("TestASMs3D.FaceCorners")
{
  ASM3DTests<ASMCube>::FaceCorners();
}


TEST_CASE("TestASMs3D.FaceIntegrate")
{
  ASM3DTests<ASMCube>::FaceIntegrate();
}


TEST_CASE("TestASMs3D.GetElementCorners")
{
  ASM3DTests<ASMTrap3<ASMs3D>>::GetElementCorners();
}


TEST_CASE("TestASMs3D.EvalSecondary")
{
  ASM3DTests<ASMCube>::EvalSecondary();
}


TEST_CASE("TestASMs3D.Write")
{
  ASM3DTests<ASMCube>::Write(false);
}


TEST_CASE("TestASMs3D.ConstrainFaceLocal")
{
  // A linear function of x, prescribed in local axes on a face, is reproduced
  // by the projection in each node, and not only in the first one
  MPCLess::compareSlaveDofOnly = true;
  ASMbase::resetNumbering();
  ASMCube pch;
  REQUIRE(pch.uniformRefine(0,2));
  REQUIRE(pch.uniformRefine(1,2));
  REQUIRE(pch.uniformRefine(2,2));
  REQUIRE(pch.generateFEMTopology());
  const size_t nnod = pch.getNoNodes();
  REQUIRE(pch.constrainFaceLocal(-3,false,1,2,true) == 16);

  LinearXFunc f(1.0);
  REQUIRE(pch.updateDirichlet({{2,&f}},std::map<int,VecFunc*>(),0.0,nullptr,false));
  size_t prescribed = 0;
  for (size_t n = nnod+1; n <= pch.getNoNodes(); n++)
    if (const MPC* mpc = pch.findMPC(pch.getNodeID(n),1); mpc && mpc->getNoMaster() == 0)
    {
      INFO("node " << n << " at " << pch.getCoord(n));
      CHECK_THAT(mpc->getSlave().coeff, WithinAbs(pch.getCoord(n).x, 1e-12));
      ++prescribed;
    }
  CHECK(prescribed == 16);
}


TEST_CASE("TestASMs3D.PairedDirection")
{
  // Faces 1 and 2 hold u constant, 3 and 4 hold v, 5 and 6 hold w. Patches
  // meeting at opposite ends of their boxes run that parameter the same way
  for (int face = 1; face <= 6; face++)
  {
    const int dir = (face+1)/2;
    const int facing = face%2 ? face+1 : face-1;
    REQUIRE(ASMbase::pairedDirection(face,facing,0,dir) == dir);
    REQUIRE(ASMbase::pairedDirection(face,face,0,dir) == -dir);
  }

  // A face of constant u has v and w along it, in that order. Without the
  // swap each pairs with the same parameter of the neighbour, and the two
  // digits of the orientation reverse the first and the second
  REQUIRE(ASMbase::pairedDirection(2,1,0,2) ==  2);
  REQUIRE(ASMbase::pairedDirection(2,1,0,3) ==  3);
  REQUIRE(ASMbase::pairedDirection(2,1,2,2) == -2);
  REQUIRE(ASMbase::pairedDirection(2,1,2,3) ==  3);
  REQUIRE(ASMbase::pairedDirection(2,1,1,2) ==  2);
  REQUIRE(ASMbase::pairedDirection(2,1,1,3) == -3);

  // The swap pairs each with the other one, which is a tie this patch
  // cannot make: the bases are one per parameter and a connection ties a
  // basis to the basis of the same number
  REQUIRE(ASMbase::pairedDirection(2,1,4,2) ==  3);
  REQUIRE(ASMbase::pairedDirection(2,1,4,3) ==  2);

  // A face of constant v has u and w along it, and one of constant w has u
  // and v, the normal never being among them
  REQUIRE(ASMbase::pairedDirection(4,3,0,1) == 1);
  REQUIRE(ASMbase::pairedDirection(4,3,0,3) == 3);
  REQUIRE(ASMbase::pairedDirection(6,5,0,1) == 1);
  REQUIRE(ASMbase::pairedDirection(6,5,0,2) == 2);

  // Half a turn about the y axis: the second patch meets the interface at
  // the same end of its box, with the first and third parameters reversed
  REQUIRE(ASMbase::pairedDirection(2,2,1,1) == -1);
  REQUIRE(ASMbase::pairedDirection(2,2,1,2) ==  2);
  REQUIRE(ASMbase::pairedDirection(2,2,1,3) == -3);
}
