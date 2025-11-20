//==============================================================================
//!
//! \file ASM2DTests.h
//!
//! \date Nov 20 2025
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for 2D FE models.
//!
//==============================================================================

#include "ASMbase.h"
#include "MPCLess.h"
#include "SIM2D.h"
#include "Vec3Oper.h"

#include "Catch2Support.h"

#include <sstream>


template<class Patch>
class ASM2DTests
{
private:
  struct EdgeTest
  {
    int edge;
    int edgeIdx;
    std::array<int,2> c1;
    std::array<int,2> c2;
  };

public:
  static void BoundaryElements()
  {
    const auto ref = std::array{
      std::array{0, 2},
      std::array{1, 3},
      std::array{0, 1},
      std::array{2, 3},
    };
    Patch pch1;
    ASMbase::resetNumbering();
    REQUIRE(pch1.uniformRefine(0,1));
    REQUIRE(pch1.uniformRefine(1,1));
    REQUIRE(pch1.generateFEMTopology());

    std::array<IntVec,4> n;
    for (size_t i = 1; i <= 4; ++i) {
      pch1.getBoundaryElms(i, n[i-1]);
      REQUIRE(n[i-1].size() == ref[i-1].size());
      for (size_t j = 0; j < ref[i-1].size(); ++j)
        REQUIRE(n[i-1][j] == ref[i-1][j]);
    }
  }

  static void Connect(ASM::Discretization disc)
  {
    const int param = GENERATE(0,1);
    SECTION("Orient " + std::to_string(param))
    {
      SIM2D sim(1);
      sim.opt.discretization = disc;
      REQUIRE(sim.read(("src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient" +
                        std::to_string(param) + ".xinp").c_str()));
      REQUIRE(sim.createFEMmodel());
    }
  }

  static void ConstrainEdge()
  {
    const EdgeTest param = GENERATE(
      EdgeTest{1, -1, {{-1, -1}}, {{-1 , 1}}},
      EdgeTest{2,  1, {{ 1, -1}}, {{ 1 , 1}}},
      EdgeTest{3, -2, {{-1, -1}}, {{ 1, -1}}},
      EdgeTest{4,  2, {{-1,  1}}, {{ 1,  1}}}
    );

    MPCLess::compareSlaveDofOnly = true;

    SECTION("Edge" + std::to_string(param.edge))
    {
      ASMbase::resetNumbering();
      Patch pch;
      REQUIRE(pch.uniformRefine(0, 2));
      REQUIRE(pch.uniformRefine(1, 2));
      REQUIRE(pch.generateFEMTopology());

      pch.constrainEdge(param.edgeIdx, false, 1, 1, 1);

      std::vector<int> glbNodes;
      pch.getBoundaryNodes(param.edge, glbNodes, 1);

      for (int node : glbNodes)
        REQUIRE(pch.findMPC(node,1) != nullptr);
    }
  }

  static void ConstrainEdgeOpen()
  {
    const EdgeTest param = GENERATE(
      EdgeTest{1, -1, {{-1, -1}}, {{-1 , 1}}},
      EdgeTest{2,  1, {{ 1, -1}}, {{ 1 , 1}}},
      EdgeTest{3, -2, {{-1, -1}}, {{ 1, -1}}},
      EdgeTest{4,  2, {{-1,  1}}, {{ 1,  1}}}
    );

    MPCLess::compareSlaveDofOnly = true;

    SECTION("Edge" + std::to_string(param.edge))
    {
      ASMbase::resetNumbering();
      Patch pch;
      REQUIRE(pch.uniformRefine(0, 2));
      REQUIRE(pch.uniformRefine(1, 2));
      REQUIRE(pch.generateFEMTopology());
      pch.constrainEdge(param.edgeIdx, true, 1, 1, 1);

      std::vector<int> glbNodes;
      pch.getBoundaryNodes(param.edge, glbNodes, 1);

      int crn = pch.getCorner(param.c1[0], param.c1[1], 1);
      REQUIRE(pch.findMPC(crn,1) == nullptr);
      glbNodes.erase(std::find(glbNodes.begin(), glbNodes.end(), crn));
      crn = pch.getCorner(param.c2[0], param.c2[1], 1);
      REQUIRE(pch.findMPC(crn,1) == nullptr);
      glbNodes.erase(std::find(glbNodes.begin(), glbNodes.end(), crn));

      for (int node : glbNodes)
        REQUIRE(pch.findMPC(node,1) != nullptr);
    }
  }

  template<class Ref, class RefProj>
  static void ElmNodes(const Ref& ref, const RefProj& ref_proj)
  {
    ASMbase::resetNumbering();

    Patch pch1;
    pch1.createProjectionBasis(true);
    pch1.raiseOrder(1,1);
    pch1.uniformRefine(0,1);
    pch1.uniformRefine(1,1);
    pch1.createProjectionBasis(false);
    REQUIRE(pch1.uniformRefine(0,1));
    REQUIRE(pch1.uniformRefine(1,1));
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

  template<class Ref>
  static void GetElementConnectivities(const Ref& ref)
  {
    Patch pch1;
    ASMbase::resetNumbering();
    REQUIRE(pch1.uniformRefine(0,1));
    REQUIRE(pch1.uniformRefine(1,1));
    REQUIRE(pch1.generateFEMTopology());
    const size_t nel = pch1.getNoElms();
    pch1.shiftElemNumbers(nel);
    IntMat neighGlb(2*nel), neighLoc(nel);
    pch1.getElmConnectivities(neighGlb);
    pch1.getElmConnectivities(neighLoc, ASM::GEOMETRY_BASIS);
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

  static void GetElementCorners()
  {
    Patch pch;
    REQUIRE(pch.raiseOrder(1,1));
    REQUIRE(pch.uniformRefine(0, 1));
    REQUIRE(pch.uniformRefine(1, 1));
    REQUIRE(pch.generateFEMTopology());

    struct Ref {
      int iel, c1, c2;
      std::array<double, 2> u, v;
      std::array<Vec3,4> XC;
    };

    auto makePtArray = [](const std::array<double,2>& u,
                          const std::array<double,2>& v)
    {
      return std::array{
        u[0], v[0],
        u[1], v[0],
        u[0], v[1],
        u[1], v[1]
      };
    };

    const Ref param = GENERATE(
      Ref{
        1, 4, 1,
        {0.0, 0.5}, {0.0, 0.5},
        {Vec3{0.0, 0.0}, Vec3{3.0, 0.0}, Vec3{0.5, 1.5}, Vec3{2.75, 1.5}},
      },
      Ref{
        2, 3, 2,
        {0.5, 1.0}, {0.0, 0.5},
        {Vec3{3.0, 0.0}, Vec3{6.0, 0.0}, Vec3{2.75, 1.5}, Vec3{5.0, 1.5}},
      },
      Ref{
        3, 4, 1,
        {0.0, 0.5}, {0.5, 1.0},
        {Vec3{0.5, 1.5}, Vec3{2.75, 1.5}, Vec3{1.0, 3.0}, Vec3{2.5, 3.0}},
      },
      Ref{
        4, 3, 2,
        {0.5, 1.0}, {0.5, 1.0},
        {Vec3{2.75, 1.5}, Vec3{5.0, 1.5}, Vec3{2.5, 3.0}, Vec3{4.0, 3.0}},
      }
    );

    SECTION("Elem " + std::to_string(param.iel))
    {
      const auto& [size, XC, prm] = pch.getElementCorners(param.iel);
      REQUIRE(XC.size() == 4);
      REQUIRE(prm.size() == 8);
      const std::array<double,8> ref_prm = makePtArray(param.u, param.v);
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
    REQUIRE(str.str() == Patch::square);

    REQUIRE(!pch1.write(str, 2));

    str.str("");
    REQUIRE(pch1.write(str, ASM::GEOMETRY_BASIS));
    REQUIRE(str.str() == Patch::square);

    str.str("");
    REQUIRE(pch1.write(str, ASM::PROJECTION_BASIS));
    REQUIRE(str.str() == Patch::square);

    str.str("");
    REQUIRE(!pch1.write(str, ASM::PROJECTION_BASIS_2));
    REQUIRE(pch1.write(str, ASM::REFINEMENT_BASIS) == expectRef);
    if (expectRef)
      REQUIRE(str.str() == Patch::square);

    str.str("");
    REQUIRE(pch1.write(str, ASM::INTEGRATION_BASIS));
    REQUIRE(str.str() == Patch::square);
  }
};
