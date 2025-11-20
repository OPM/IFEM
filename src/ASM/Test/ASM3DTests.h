//==============================================================================
//!
//! \file ASM3DTests.h
//!
//! \date Nov 20 2025
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for 3D FE models.
//!
//==============================================================================

#include "ASMbase.h"
#include "GlobalIntegral.h"
#include "IntegrandBase.h"
#include "MPCLess.h"
#include "SIM3D.h"
#include "Vec3Oper.h"

#include "Catch2Support.h"

#include <sstream>


namespace {

class NoProblem : public IntegrandBase
{
public:
  NoProblem() : IntegrandBase(3) {}
  virtual ~NoProblem() {}
protected:
  bool evalBou(LocalIntegral&, const FiniteElement&,
               const Vec3&, const Vec3&) const override
  { return true; }
};

}


template<class Patch>
class ASM3DTests
{
public:
  static void BoundaryElements()
  {
    const auto ref = std::array{
      std::array{0, 2, 4, 6},
      std::array{1, 3, 5, 7},
      std::array{0, 1, 4, 5},
      std::array{2, 3, 6, 7},
      std::array{0, 1, 2, 3},
      std::array{4, 5, 6, 7},
    };

    Patch pch1;
    ASMbase::resetNumbering();
    REQUIRE(pch1.uniformRefine(0,1));
    REQUIRE(pch1.uniformRefine(1,1));
    REQUIRE(pch1.uniformRefine(2,1));
    REQUIRE(pch1.generateFEMTopology());

    std::array<IntVec,6> n;
    for (size_t i = 1; i <= 6; ++i) {
      pch1.getBoundaryElms(i, n[i-1]);
      REQUIRE(n[i-1].size() == ref[i-1].size());
      for (size_t j = 0; j < ref[i-1].size(); ++j)
        REQUIRE(n[i-1][j] == ref[i-1][j]);
    }
  }

  template<class Ref>
  static void BoundaryNodes(const Ref& ref)
  {
    ASMbase::resetNumbering();
    Patch pch;
    REQUIRE(pch.uniformRefine(0, 2));
    REQUIRE(pch.uniformRefine(1, 2));
    REQUIRE(pch.uniformRefine(2, 2));
    REQUIRE(pch.generateFEMTopology());

    const int param = GENERATE(1,2,3,4,5,6);

    SECTION("Face " + std::to_string(param))
    {
      std::vector<int> vec;
      pch.getBoundaryNodes(param, vec, 1, 1, -1, false);
      REQUIRE(vec.size() == 16);
      size_t k = 0;
      for (int j = 0; j < 4; ++j)
        for (int i = 0; i < 4; ++i, ++k)
          REQUIRE(vec[k] == ref[param-1](i,j));
    }
  }

  static void Connect(ASM::Discretization disc)
  {
    const int param = GENERATE(0,1,2,3,4,5,6,7);

    SECTION("Orient " + std::to_string(param)) {
      SIM3D sim(3);
      sim.opt.discretization = disc;
      REQUIRE(sim.read(("src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient" +
                        std::to_string(param) + ".xinp").c_str()));
      REQUIRE(sim.createFEMmodel());
    }
  }

  static void ConnectUneven(ASM::Discretization disc)
  {
    const int param = GENERATE(0,1,2,3,4,5,6,7,8,9,10,11,12);

    SECTION("Uneven " + std::to_string(param)) {
      SIM3D sim(1);
      sim.opt.discretization = disc;
      std::stringstream str;
      str << "src/ASM/Test/refdata/3d_uneven";
      int idx = param-1;
      if (idx >= 0)
        str << "_" << idx/3 << idx%3;
      str << ".xinp";
      REQUIRE(sim.read(str.str().c_str()));
      REQUIRE(sim.createFEMmodel());
    }
  }

  static void ConstrainFace()
  {
    struct FaceTest
    {
      int face;
      int dir;
    };
    const FaceTest param = GENERATE(
      FaceTest{1, -1},
      FaceTest{2,  1},
      FaceTest{3, -2},
      FaceTest{4,  2},
      FaceTest{5, -3},
      FaceTest{6,  3}
    );

    MPCLess::compareSlaveDofOnly = true;

    SECTION("Face" + std::to_string(param.face))
    {
      ASMbase::resetNumbering();
      Patch pch;
      REQUIRE(pch.uniformRefine(0, 2));
      REQUIRE(pch.uniformRefine(1, 2));
      REQUIRE(pch.uniformRefine(2, 2));
      REQUIRE(pch.generateFEMTopology());

      pch.constrainFace(param.dir, false, 1, 1, 1);

      std::vector<int> glbNodes;
      pch.getBoundaryNodes(param.face, glbNodes, 1, 1, false, false);

      for (int node : glbNodes)
        REQUIRE(pch.findMPC(node,1) != nullptr);
    }
  }

  template<class Ref>
  static void ElementConnectivities(const Ref& ref)
  {
    Patch pch1;
    ASMbase::resetNumbering();
    REQUIRE(pch1.uniformRefine(0,1));
    REQUIRE(pch1.uniformRefine(1,1));
    REQUIRE(pch1.uniformRefine(2,1));
    REQUIRE(pch1.generateFEMTopology());
    const size_t nel = pch1.getNoElms();
    pch1.shiftElemNumbers(nel);
    IntMat neighGlb(2*nel), neighLoc(nel);
    pch1.getElmConnectivities(neighGlb);
    pch1.getElmConnectivities(neighLoc, ASM::GEOMETRY_BASIS);
    REQUIRE(neighGlb.size() == 2*nel);
    REQUIRE(neighLoc.size() == nel);
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

  template<class Ref, class RefProj>
  static void ElmNodes(const Ref& ref, const RefProj& ref_proj)
  {
    ASMbase::resetNumbering();
    Patch pch1;
    REQUIRE(pch1.createProjectionBasis(true));
    REQUIRE(pch1.raiseOrder(1,1,1,false));
    REQUIRE(pch1.uniformRefine(0,1));
    REQUIRE(pch1.uniformRefine(1,1));
    REQUIRE(pch1.uniformRefine(2,1));
    REQUIRE(pch1.createProjectionBasis(false));
    REQUIRE(pch1.uniformRefine(0,1));
    REQUIRE(pch1.uniformRefine(1,1));
    REQUIRE(pch1.uniformRefine(2,1));
    REQUIRE(pch1.generateFEMTopology());

    const IntMat mnpc = pch1.getElmNodes(1);

    REQUIRE(mnpc.size() == ref.size());
    for (size_t i = 0; i < mnpc.size(); ++i) {
      REQUIRE(mnpc[i].size() == ref[i].size());
      for (size_t j = 0; j < mnpc[i].size(); ++j)
        REQUIRE(mnpc[i][j] == ref[i][j]);
    }

    const IntMat mnpc_proj = pch1.getElmNodes(ASM::PROJECTION_BASIS);
    REQUIRE(mnpc_proj.size() == ref_proj.size());
    for (size_t i = 0; i < mnpc_proj.size(); ++i) {
      REQUIRE(mnpc_proj[i].size() == ref_proj[i].size());
      for (size_t j = 0; j < mnpc_proj[i].size(); ++j)
        REQUIRE(mnpc_proj[i][j] == ref_proj[i][j]);
    }
  }

  static void FaceCorners()
  {
    const auto refArr = std::array{
      std::array{1, 2, 3, 4}, // Bottom
      std::array{1, 2, 5, 6}, // South
      std::array{1, 3, 5, 7}, // West
      std::array{0, 0, 0, 0}, // None
      std::array{2, 4, 6, 8}, // East
      std::array{3, 4, 7, 8}, // North
      std::array{5, 6, 7, 8}, // Top
    };

    Patch pch1;
    ASMbase::resetNumbering();
    REQUIRE(pch1.generateFEMTopology());

    int dir = GENERATE(-1,1,-2,2,-3,3);
    SECTION("Direction "+ std::to_string(dir))
    {
      REQUIRE(pch1.getFaceCorners(dir) == refArr[dir+3]);
    }
  }


  static void FaceIntegrate()
  {
    GlobalIntegral dummy;
    NoProblem prb;

    SECTION("Uneven order")
    {
      ASMbase::resetNumbering();
      Patch patch;
      ASMbase::resetNumbering();
      REQUIRE(patch.raiseOrder(2,1,0,false));
      REQUIRE(patch.generateFEMTopology());
      for (int lIndex = 1; lIndex <= 6; lIndex++)
      {
        patch.generateThreadGroups(lIndex,false,false);
        REQUIRE(patch.integrate(prb,lIndex,dummy,TimeDomain()));
      }
    }

    SECTION("Uneven refine")
    {
      ASMbase::resetNumbering();
      Patch patch;

      REQUIRE(patch.uniformRefine(0,3));
      REQUIRE(patch.uniformRefine(1,2));
      REQUIRE(patch.uniformRefine(2,1));
      REQUIRE(patch.generateFEMTopology());
      for (int lIndex = 1; lIndex <= 6; lIndex++)
      {
        patch.generateThreadGroups(lIndex,false,false);
        REQUIRE(patch.integrate(prb,lIndex,dummy,TimeDomain()));
      }
    }
  }

  static void GetElementCorners()
  {
    Patch pch;
    REQUIRE(pch.raiseOrder(1,1,1,false));
    REQUIRE(pch.uniformRefine(0, 1));
    REQUIRE(pch.uniformRefine(1, 1));
    REQUIRE(pch.uniformRefine(2, 1));
    REQUIRE(pch.generateFEMTopology());

    struct Ref {
      int iel, c1, c2;
      std::array<double,2> u, v, w;
      std::array<Vec3,8> XC;
    };

    auto makePtArray = [](const std::array<double,2>& u,
                          const std::array<double,2>& v,
                          const std::array<double,2>& w)
    {
      return std::array{
        u[0], v[0], w[0],
        u[1], v[0], w[0],
        u[0], v[1], w[0],
        u[1], v[1], w[0],
        u[0], v[0], w[1],
        u[1], v[0], w[1],
        u[0], v[1], w[1],
        u[1], v[1], w[1],
      };
    };

    const Ref param = GENERATE(
        Ref{
          1, 8, 1,
          {0.0, 0.5}, {0.0, 0.5}, {0.0, 0.5},
          {Vec3{0.0, 0.0, 0.0}, Vec3{3.0, 0.0, 0.0},
           Vec3{0.5, 1.5, 0.0}, Vec3{2.75, 1.5, 0.0},
           Vec3{0.5, 0.0, 0.5}, Vec3{3.0, 0.0, 0.5},
           Vec3{1.0, 1.5, 0.5}, Vec3{2.875, 1.5, 0.5}},
        },
        Ref{
          2, 7, 2,
          {0.5, 1.0}, {0.0, 0.5}, {0.0, 0.5},
          {Vec3{3.0, 0.0, 0.0}, Vec3{6.0, 0.0, 0.0},
           Vec3{2.75, 1.5, 0.0}, Vec3{5, 1.5, 0.0},
           Vec3{3.0, 0.0, 0.5}, Vec3{5.5, 0.0, 0.5},
           Vec3{2.875, 1.5, 0.5}, Vec3{4.75, 1.5, 0.5}},
        },
        Ref{
          3, 8, 1,
          {0.0, 0.5}, {0.5, 1.0}, {0.0, 0.5},
          {Vec3{0.5, 1.5, 0.0}, Vec3{2.75, 1.5, 0.0},
           Vec3{1.0, 3.0, 0.0}, Vec3{2.5, 3.0, 0.0},
           Vec3{1.0, 1.5, 0.5}, Vec3{2.875, 1.5, 0.5},
           Vec3{1.5, 3.0, 0.5}, Vec3{2.75, 3.0, 0.5}},
        },
        Ref{
          4, 7, 2,
          {0.5, 1.0}, {0.5, 1.0}, {0.0, 0.5},
          {Vec3{2.75, 1.5, 0.0}, Vec3{5.0, 1.5, 0.0},
           Vec3{2.5, 3.0, 0.0}, Vec3{4.0, 3.0, 0.0},
           Vec3{2.875, 1.5, 0.5}, Vec3{4.75, 1.5, 0.5},
           Vec3{2.75, 3.0, 0.5}, Vec3{4.0, 3.0, 0.5}},
        },
        Ref{
          5, 8, 1,
          {0.0, 0.5}, {0.0, 0.5}, {0.5, 1.0},
          {Vec3{0.5, 0.0, 0.5}, Vec3{3.0, 0.0, 0.5},
           Vec3{1.0, 1.5, 0.5}, Vec3{2.875, 1.5, 0.5},
           Vec3{1.0, 0.0, 1.0}, Vec3{3.0, 0.0, 1.0},
           Vec3{1.5, 1.5, 1.0}, Vec3{3.0, 1.5, 1.0}},
        },
        Ref{
          6, 7, 2,
          {0.5, 1.0}, {0.0, 0.5}, {0.5, 1.0},
          {Vec3{3.0, 0.0, 0.5}, Vec3{5.5, 0.0, 0.5},
           Vec3{2.875, 1.5, 0.5}, Vec3{4.75, 1.5, 0.5},
           Vec3{3.0, 0.0, 1.0}, Vec3{5.0, 0.0, 1.0},
           Vec3{3.0, 1.5, 1.0}, Vec3{4.5, 1.5, 1.0}},
        },
        Ref{
          7, 8, 1,
          {0.0, 0.5}, {0.5, 1.0}, {0.5, 1.0},
          {Vec3{1.0, 1.5, 0.5}, Vec3{2.875, 1.5, 0.5},
           Vec3{1.5, 3.0, 0.5}, Vec3{2.75, 3.0, 0.5},
           Vec3{1.5, 1.5, 1.0}, Vec3{3.0, 1.5, 1.0},
           Vec3{2.0, 3.0, 1.0}, Vec3{3.0, 3.0, 1.0}},
        },
        Ref{
          8, 6, 3,
          {0.5, 1.0}, {0.5, 1.0}, {0.5, 1.0},
          {Vec3{2.875, 1.5, 0.5}, Vec3{4.75, 1.5, 0.5},
           Vec3{2.75, 3.0, 0.5}, Vec3{4.0, 3.0, 0.5},
           Vec3{3.0, 1.5, 1.0}, Vec3{4.5, 1.5, 1.0},
           Vec3{3.0, 3.0, 1.0}, Vec3{4.0, 3.0, 1.0}},
        }
    );

    SECTION("Elem " + std::to_string(param.iel))
    {
      const auto& [size, XC, prm] = pch.getElementCorners(param.iel);
      REQUIRE(XC.size() == 8);
      REQUIRE(prm.size() == 24);
      const std::array<double,24> ref_prm = makePtArray(param.u, param.v, param.w);
      REQUIRE(prm.size() == ref_prm.size());
      for (size_t i = 0; i < prm.size(); ++i)
        REQUIRE_THAT(prm[i], WithinRel(ref_prm[i]));
      REQUIRE(XC.size() == param.XC.size());
      for (size_t i = 0; i < XC.size(); ++i)
        REQUIRE(XC[i] == param.XC[i]);
    }
  }

  static void Write(const bool expectRef)
  {
    ASMbase::resetNumbering();
    Patch pch1;
    REQUIRE(pch1.generateFEMTopology());

    std::stringstream str;
    REQUIRE(pch1.write(str, 1));
    REQUIRE(str.str() == Patch::cube);

    REQUIRE(!pch1.write(str, 2));

    str.str("");
    REQUIRE(pch1.write(str, ASM::GEOMETRY_BASIS));
    REQUIRE(str.str() == Patch::cube);

    str.str("");
    REQUIRE(pch1.write(str, ASM::PROJECTION_BASIS));
    REQUIRE(str.str() == Patch::cube);

    REQUIRE(!pch1.write(str, ASM::PROJECTION_BASIS_2));

    str.str("");
    REQUIRE(pch1.write(str, ASM::REFINEMENT_BASIS) == expectRef);
    if (expectRef)
      REQUIRE(str.str() == Patch::cube);

    str.str("");
    REQUIRE(pch1.write(str, ASM::INTEGRATION_BASIS));
    REQUIRE(str.str() == Patch::cube);
  }
};
