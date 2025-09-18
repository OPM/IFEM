//==============================================================================
//!
//! \file TestSIMNodalConstraint.C
//!
//! \date Nov 26 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for base class for simulators constraining a topologyset
//!        to a given node.
//!
//==============================================================================

#include "SIMNodalConstraint.h"
#include "IFEM.h"
#include "MPC.h"
#include "SIM2D.h"

#include "LR/ASMu2D.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

namespace {

auto&& check_mpc = [](MPC* mpc, int node)
{
  REQUIRE(mpc != nullptr);
  REQUIRE(mpc->getNoMaster() == 1);
  REQUIRE(mpc->getMaster(0).node == node);
  REQUIRE(mpc->getMaster(0).dof == 1);
};

}


TEST_CASE("TestSIMNodalConstraint.Vertex2DLR")
{
  const int param = GENERATE(1,2,3,4);

  SECTION("Vertex " + std::to_string(param)) {
    IFEM::getOptions().discretization = ASM::LRSpline;
    SIMNodalConstraint<SIM2D> s({2});
    std::stringstream str;
    str << "refdata/nodal_2D_V" << param << ".xinp";
    REQUIRE(s.read(str.str().c_str()));
    s.preprocess();

    const ASMu2D& pch = static_cast<const ASMu2D&>(*s.getPatch(1));
    for (size_t i = 1; i <= pch.getNoNodes(); ++i) {
      if (i == 1 && param == 1)
        check_mpc(pch.findMPC(i, 1), pch.getCorner(1,-1,1));
      else if (i == static_cast<size_t>(pch.getCorner(1,-1,1)) && param == 2)
        check_mpc(pch.findMPC(i, 1), 1);
      else if (i == static_cast<size_t>(pch.getCorner(-1,1,1)) && param == 3)
        check_mpc(pch.findMPC(i, 1), 1);
      else if (i == static_cast<size_t>(pch.getCorner(1,1,1)) && param == 4)
        check_mpc(pch.findMPC(i, 1), 1);
      else
        REQUIRE(pch.findMPC(i, 1) == nullptr);
      REQUIRE(pch.findMPC(i, 2) == nullptr);
    }
  }
}


TEST_CASE("TestSIMNodalConstraint.Edge2DLR")
{
  const int param = GENERATE(1,2,3,4);

  SECTION("Edge " + std::to_string(param)) {
    IFEM::getOptions().discretization = ASM::LRSpline;
    SIMNodalConstraint<SIM2D> s({2});
    std::stringstream str;
    str << "refdata/nodal_2D_E" << param << ".xinp";
    REQUIRE(s.read(str.str().c_str()));
    s.preprocess();

    const ASMu2D& pch = static_cast<const ASMu2D&>(*s.getPatch(1));
    IntVec nodes;
    pch.getBoundaryNodes(param,nodes,1);
    for (size_t i = 1; i <= pch.getNoNodes(); ++i) {
      if (std::find(nodes.begin(), nodes.end(), i) != nodes.end() &&
          i != static_cast<size_t>(pch.getCorner(1,-1,1)) && param == 1)
        check_mpc(pch.findMPC(i, 1), pch.getCorner(1,-1,1));
      else if (std::find(nodes.begin(), nodes.end(), i) != nodes.end() && param == 3)
        check_mpc(pch.findMPC(i, 1), pch.getCorner(1,1,1));
      else if (std::find(nodes.begin(), nodes.end(), i) != nodes.end())
        check_mpc(pch.findMPC(i, 1), 1);
      else
        REQUIRE(pch.findMPC(i, 1) == nullptr);
      REQUIRE(pch.findMPC(i, 2) == nullptr);
    }
  }
}


TEST_CASE("TestSIMNodalConstraint.Edge2DLRmx")
{
  const int param = GENERATE(1,2,3,4);

  SECTION("Edge " + std::to_string(param)) {
    IFEM::getOptions().discretization = ASM::LRSpline;
    SIMNodalConstraint<SIM2D> s({2,2});
    std::stringstream str;
    str << "refdata/nodal_2D_E" << param << "_mixed.xinp";
    REQUIRE(s.read(str.str().c_str()));
    s.preprocess();

    const ASMu2D& pch = static_cast<const ASMu2D&>(*s.getPatch(1));
    IntVec nodes;
    pch.getBoundaryNodes(param,nodes,2);
    size_t ofs = pch.getNoNodes(1);
    for (size_t i = 1; i <= pch.getNoNodes(1); ++i) {
      REQUIRE(pch.findMPC(i, 1) == nullptr);
      REQUIRE(pch.findMPC(i, 2) == nullptr);
    }
    for (size_t i = 1; i <= pch.getNoNodes(2); ++i) {
      if (std::find(nodes.begin(), nodes.end(), i+ofs) != nodes.end() &&
          i+ofs != static_cast<size_t>(pch.getCorner(1,-1,2)) && param == 1)
        check_mpc(pch.findMPC(i+ofs, 1), pch.getCorner(1,-1,2));
      else if (std::find(nodes.begin(), nodes.end(), i+ofs) != nodes.end() && param == 3)
        check_mpc(pch.findMPC(i+ofs, 1), pch.getCorner(1,1,2));
      else if (std::find(nodes.begin(), nodes.end(), i+ofs) != nodes.end())
        check_mpc(pch.findMPC(i+ofs, 1), ofs+1);
      else
        REQUIRE(pch.findMPC(i+ofs, 1) == nullptr);
      REQUIRE(pch.findMPC(i+ofs, 2) == nullptr);
    }
  }
}
