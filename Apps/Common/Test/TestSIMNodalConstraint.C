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
#include "ASMs2D.h"
#include "ASMs3D.h"
#include "ASMs3Dmx.h"
#include "MPC.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"

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


TEST_CASE("TestSIMNodalConstraint.Vertex1D")
{
  const int param = GENERATE(1,2);

  SECTION("Vertex " + std::to_string(param)) {
    SIMNodalConstraint<SIM1D> s({2});
    std::stringstream str;
    str << "refdata/nodal_1D_V" << param << ".xinp";
    REQUIRE(s.read(str.str().c_str()));
    s.preprocess();

    const ASMbase& pch = *s.getPatch(1);
    for (size_t i=1; i <= pch.getNoNodes(); ++i) {
      if (i == 1 && param == 1)
        check_mpc(pch.findMPC(i,1), pch.getNoNodes());
      else if (i == pch.getNoNodes() && param == 2)
        check_mpc(pch.findMPC(i,1), 1);
      else
        REQUIRE(pch.findMPC(i, 1) == nullptr);
      REQUIRE(pch.findMPC(i, 2) == nullptr);
    }
  }
}


TEST_CASE("TestSIMNodalConstraint.Vertex2D")
{
  const int param = GENERATE(1,2,3,4);

  SECTION("Vertex " + std::to_string(param)) {
    SIMNodalConstraint<SIM2D> s({2});
    std::stringstream str;
    str << "refdata/nodal_2D_V" << param << ".xinp";
    REQUIRE(s.read(str.str().c_str()));
    s.preprocess();

    const ASMs2D& pch = static_cast<const ASMs2D&>(*s.getPatch(1));
    int n1, n2;
    pch.getSize(n1,n2);
    for (size_t i = 1; i <= pch.getNoNodes(); ++i) {
      if (i == 1 && param == 1)
        check_mpc(pch.findMPC(i, 1), n1);
      else if (i == static_cast<size_t>(n1) && param == 2)
        check_mpc(pch.findMPC(i, 1), 1);
      else if (i == static_cast<size_t>(n1*(n2-1)+1) && param == 3)
        check_mpc(pch.findMPC(i, 1), 1);
        else if (i == pch.getNoNodes() && param == 4)
        check_mpc(pch.findMPC(i, 1), 1);
      else
        REQUIRE(pch.findMPC(i, 1) == nullptr);
      REQUIRE(pch.findMPC(i, 2) == nullptr);
    }
  }
}


TEST_CASE("TestSIMNodalConstraint.Edge2D")
{
  const int param = GENERATE(1,2,3,4);

  SECTION("Edge " + std::to_string(param)) {
    SIMNodalConstraint<SIM2D> s({2});
    std::stringstream str;
    str << "refdata/nodal_2D_E" << param << ".xinp";
    REQUIRE(s.read(str.str().c_str()));
    s.preprocess();

    const ASMs2D& pch = static_cast<const ASMs2D&>(*s.getPatch(1));
    int n1, n2;
    pch.getSize(n1,n2);
    for (size_t i = 1; i <= pch.getNoNodes(); ++i) {
      if (i % n1 == 1 && param == 1)
        check_mpc(pch.findMPC(i, 1), n1);
      else if (i % n1 == 0 && param == 2)
        check_mpc(pch.findMPC(i, 1), 1);
      else if (i <= static_cast<size_t>(n1) && param == 3)
        check_mpc(pch.findMPC(i, 1), n1*n2);
      else if (i > static_cast<size_t>(n1*(n2-1)) && param == 4)
        check_mpc(pch.findMPC(i, 1), 1);
      else
        REQUIRE(pch.findMPC(i, 1) == nullptr);
      REQUIRE(pch.findMPC(i, 2) == nullptr);
    }
  }
}


TEST_CASE("TestSIMNodalConstraint.Edge2Dmx")
{
  const int param = GENERATE(1,2,3,4);

  SECTION("Edge " + std::to_string(param)) {
    SIMNodalConstraint<SIM2D> s({2,2});
    std::stringstream str;
    str << "refdata/nodal_2D_E" << param << "_mixed.xinp";
    REQUIRE(s.read(str.str().c_str()));
    s.preprocess();

    const ASMs2D& pch = static_cast<const ASMs2D&>(*s.getPatch(1));
    int n1, n2;
    pch.getSize(n1,n2,2);
    size_t ofs = pch.getNoNodes(1);
    for (size_t i = 1; i <= pch.getNoNodes(1); ++i) {
      REQUIRE(pch.findMPC(i, 1) == nullptr);
      REQUIRE(pch.findMPC(i, 2) == nullptr);
    }
    for (size_t i = 1; i <= pch.getNoNodes(2); ++i) {
      if (i % n1 == 1 && param == 1)
        check_mpc(pch.findMPC(i+ofs, 1), n1+ofs);
      else if (i % n1 == 0 && param == 2)
        check_mpc(pch.findMPC(i+ofs, 1), 1+ofs);
      else if (i <= static_cast<size_t>(n1) && param == 3)
        check_mpc(pch.findMPC(i+ofs, 1), n1*n2+ofs);
      else if (i > static_cast<size_t>(n1*(n2-1)) && param == 4)
        check_mpc(pch.findMPC(i+ofs, 1), 1+ofs);
      else
        REQUIRE(pch.findMPC(i+ofs, 1) == nullptr);
      REQUIRE(pch.findMPC(i+ofs, 2) == nullptr);
    }
  }
}


TEST_CASE("TestSIMNodalConstraint.Vertex3D")
{
  const int param = GENERATE(1,2,3,4,5,6,7,8);

  SECTION("Vertex " + std::to_string(param)) {
    SIMNodalConstraint<SIM3D> s({2});
    std::stringstream str;
    str << "refdata/nodal_3D_V" << param << ".xinp";
    REQUIRE(s.read(str.str().c_str()));
    s.preprocess();

    const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
    int n1, n2, n3;
    pch.getSize(n1,n2,n3);
    for (size_t i = 1; i <= pch.getNoNodes(); ++i) {
      if (i == 1 && param == 1)
        check_mpc(pch.findMPC(i, 1), n1);
      else if ((i == static_cast<size_t>(n1) && param == 2)
            || (i == static_cast<size_t>(n1*(n2-1)+1) && param == 3)
            || (i == static_cast<size_t>(n1*n2) && param == 4)
            || (i == static_cast<size_t>(1+n1*n2*(n3-1)) && param == 5)
            || (i == static_cast<size_t>(n1+n1*n2*(n3-1)) && param == 6)
            || (i == static_cast<size_t>(1+n1*(n2-1)+n1*n2*(n3-1)) && param == 7)
            || (i == pch.getNoNodes() && param == 8))
        check_mpc(pch.findMPC(i, 1), 1);
      else
        REQUIRE(pch.findMPC(i, 1) == nullptr);
      REQUIRE(pch.findMPC(i, 2) == nullptr);
    }
  }
}


TEST_CASE("TestSIMNodalConstraint.Edge3D")
{
  const int param = GENERATE(1,2,3,4,5,6,7,8,9,10,11,12);

  SECTION("Edge " + std::to_string(param)) {
    SIMNodalConstraint<SIM3D> s({2});
    std::stringstream str;
    str << "refdata/nodal_3D_E" << param << ".xinp";
    REQUIRE(s.read(str.str().c_str()));
    s.preprocess();

    const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
    int n1, n2, n3;
    pch.getSize(n1,n2,n3);
    for (size_t i = 1; i <= pch.getNoNodes(); ++i) {
      if (i <= static_cast<size_t>(n1) && param == 1)
        check_mpc(pch.findMPC(i, 1), n1*n2*n3);
      else if (i >= static_cast<size_t>(n1*(n2-1)+1) &&
               i <= static_cast<size_t>(n1*n2) && param == 2)
        check_mpc(pch.findMPC(i, 1), 1);
      else if (i >= static_cast<size_t>(n1*n2*(n3-1)+1) &&
               i <= static_cast<size_t>(n1*n2*(n3-1)+n1) && param == 3)
        check_mpc(pch.findMPC(i, 1), 1);
      else if (i >= static_cast<size_t>(n1*n2*(n3-1)+n1*(n2-1)+1) && param == 4)
        check_mpc(pch.findMPC(i, 1), 1);
      else if (i <= static_cast<size_t>(n1*n2) && i % n2 == 1 && param == 5)
        check_mpc(pch.findMPC(i, 1), n1*n2*n3);
      else if (i <= static_cast<size_t>(n1*n2) && i % n2 == 0 && param == 6)
        check_mpc(pch.findMPC(i, 1), 1);
      else if (i >= static_cast<size_t>(n1*n2*(n3-1)+1) && i % n2 == 1 && param == 7)
        check_mpc(pch.findMPC(i, 1), 1);
      else if (i >= static_cast<size_t>(n1*n2*(n3-1)+1) && i % n2 == 0 && param == 8)
        check_mpc(pch.findMPC(i, 1), 1);
      else if (i % (n1*n2) == 1 && param == 9)
        check_mpc(pch.findMPC(i, 1), n1*n2*n3);
      else if (i % (n1*n2) == static_cast<size_t>(n1) && param == 10)
        check_mpc(pch.findMPC(i, 1), 1);
      else if (i % (n1*n2) == static_cast<size_t>(n1*(n2-1)+1) && param == 11)
        check_mpc(pch.findMPC(i, 1), 1);
      else if (i % (n1*n2) == 0 && param == 12)
        check_mpc(pch.findMPC(i, 1), 1);
      else
        REQUIRE(pch.findMPC(i, 1) == nullptr);
      REQUIRE(pch.findMPC(i, 2) == nullptr);
    }
  }
}


TEST_CASE("TestSIMNodalConstraint.Face3D")
{
  const int param = GENERATE(1,2,3,4,5,6);

  SECTION("Face " + std::to_string(param)) {
    SIMNodalConstraint<SIM3D> s({2});
    std::stringstream str;
    str << "refdata/nodal_3D_F" << param << ".xinp";
    REQUIRE(s.read(str.str().c_str()));
    s.preprocess();

    const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
    int n1, n2, n3;
    pch.getSize(n1,n2,n3);
    for (size_t i = 1; i <= pch.getNoNodes(); ++i) {
      if (i % n1 == 1 && param == 1)
        check_mpc(pch.findMPC(i, 1), n1*n2*n3);
      else if (i % n1 == 0 && param == 2)
        check_mpc(pch.findMPC(i, 1), 1);
      else if (i % (n1*n2) > 0 && i % (n1*n2) <= static_cast<size_t>(n1) && param == 3)
        check_mpc(pch.findMPC(i, 1), n1*n2*n3);
      else if ((i % (n1*n2) >= static_cast<size_t>(n1*(n2-1)+1) || i % (n1*n2) == 0) && param == 4)
        check_mpc(pch.findMPC(i, 1), 1);
      else if (i <= static_cast<size_t>(n1*n2) && param == 5)
        check_mpc(pch.findMPC(i, 1), n1*n2*n3);
      else if (i >= static_cast<size_t>(n1*n2*(n3-1)+1) && param == 6)
        check_mpc(pch.findMPC(i, 1), 1);
      else
        REQUIRE(pch.findMPC(i, 1) == nullptr);
      REQUIRE(pch.findMPC(i, 2) == nullptr);
    }
  }
}


TEST_CASE("TestSIMNodalConstraint.ASMs3DmxE3")
{
  SIMNodalConstraint<SIM3D> s({2,2});
  REQUIRE(s.read("refdata/nodal_3D_E3_mixed.xinp"));
  s.preprocess();

  const ASMs3Dmx& pch = static_cast<const ASMs3Dmx&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3,2);
  size_t ofs = pch.getNoNodes(1);
  for (size_t i = 1; i <= pch.getNoNodes(1); ++i) {
    REQUIRE(pch.findMPC(i, 1) == nullptr);
    REQUIRE(pch.findMPC(i, 2) == nullptr);
  }
  for (size_t i = 1; i <= pch.getNoNodes(2); ++i) {
    if (i >= static_cast<size_t>(n1*n2*(n3-1)+1) &&
        i <= static_cast<size_t>(n1*n2*(n3-1)+n1))
      check_mpc(pch.findMPC(i+ofs, 1), 1+ofs);
    else
      REQUIRE(pch.findMPC(i+ofs, 1) == nullptr);
    REQUIRE(pch.findMPC(i+ofs, 2) == nullptr);
  }
}
