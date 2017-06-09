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
#ifdef HAS_LRSPLINE
#include "ASMu2D.h"
#include "LRSpline/LRSpline.h"
#endif
#include "ASMs2Dmx.h"
#include "ASMs3Dmx.h"
#include "MPC.h"

#include "gtest/gtest.h"

auto&& check_mpc = [](MPC* mpc, int node)
{
  ASSERT_TRUE(mpc != nullptr);
  ASSERT_TRUE(mpc->getNoMaster() == 1);
  ASSERT_TRUE(mpc->getMaster(0).node == node);
  ASSERT_TRUE(mpc->getMaster(0).dof == 1);
};


class TestSIMNodalConstraint : public testing::Test,
                               public testing::WithParamInterface<int>
{
};


TEST_P(TestSIMNodalConstraint, Vertex1D)
{
  if (GetParam() > 2)
    return;

  SIMNodalConstraint<SIM1D> s({2});
  std::stringstream str;
  str << "refdata/nodal_1D_V" << GetParam() << ".xinp";
  ASSERT_TRUE(s.read(str.str().c_str()));
  s.preprocess();

  const ASMbase& pch = *s.getPatch(1);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == 1 && GetParam() == 1)
      check_mpc(pch.findMPC(i,1), pch.getNoNodes());
    else if (i == pch.getNoNodes() && GetParam() == 2)
      check_mpc(pch.findMPC(i,1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST_P(TestSIMNodalConstraint, Vertex2D)
{
  if (GetParam() > 4)
    return;

  SIMNodalConstraint<SIM2D> s({2});
  std::stringstream str;
  str << "refdata/nodal_2D_V" << GetParam() << ".xinp";
  ASSERT_TRUE(s.read(str.str().c_str()));
  s.preprocess();

  const ASMs2D& pch = static_cast<const ASMs2D&>(*s.getPatch(1));
  int n1, n2;
  pch.getSize(n1,n2);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == 1 && GetParam() == 1)
      check_mpc(pch.findMPC(i, 1), n1);
    else if (i == (size_t)n1 && GetParam() == 2)
      check_mpc(pch.findMPC(i, 1), 1);
    else if (i == (size_t)n1*(n2-1)+1 && GetParam() == 3)
      check_mpc(pch.findMPC(i, 1), 1);
      else if (i == pch.getNoNodes() && GetParam() == 4)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST_P(TestSIMNodalConstraint, Edge2D)
{
  if (GetParam() > 4)
    return;

  SIMNodalConstraint<SIM2D> s({2});
  std::stringstream str;
  str << "refdata/nodal_2D_E" << GetParam() << ".xinp";
  ASSERT_TRUE(s.read(str.str().c_str()));
  s.preprocess();

  const ASMs2D& pch = static_cast<const ASMs2D&>(*s.getPatch(1));
  int n1, n2;
  pch.getSize(n1,n2);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i % n1 == 1 && GetParam() == 1)
      check_mpc(pch.findMPC(i, 1), n1);
    else if (i % n1 == 0 && GetParam() == 2)
      check_mpc(pch.findMPC(i, 1), 1);
    else if (i <= (size_t)n1 && GetParam() == 3)
      check_mpc(pch.findMPC(i, 1), n1*n2);
    else if (i > (size_t)n1*(n2-1) && GetParam() == 4)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST_P(TestSIMNodalConstraint, Edge2Dmx)
{
  if (GetParam() > 4)
    return;

  SIMNodalConstraint<SIM2D> s({2,2});
  std::stringstream str;
  str << "refdata/nodal_2D_E" << GetParam() << "_mixed.xinp";
  ASSERT_TRUE(s.read(str.str().c_str()));
  s.preprocess();

  const ASMs2D& pch = static_cast<const ASMs2D&>(*s.getPatch(1));
  int n1, n2;
  pch.getSize(n1,n2,2);
  size_t ofs = pch.getNoNodes(1);
  for (size_t i=1; i <= pch.getNoNodes(1); ++i) {
    ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
  for (size_t i=1; i <= pch.getNoNodes(2); ++i) {
    if (i % n1 == 1 && GetParam() == 1)
      check_mpc(pch.findMPC(i+ofs, 1), n1+ofs);
    else if (i % n1 == 0 && GetParam() == 2)
      check_mpc(pch.findMPC(i+ofs, 1), 1+ofs);
    else if (i <= (size_t)n1 && GetParam() == 3)
      check_mpc(pch.findMPC(i+ofs, 1), n1*n2+ofs);
    else if (i > (size_t)n1*(n2-1) && GetParam() == 4)
      check_mpc(pch.findMPC(i+ofs, 1), 1+ofs);
    else
      ASSERT_TRUE(pch.findMPC(i+ofs, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i+ofs, 2) == nullptr);
  }
}


#ifdef HAS_LRSPLINE
TEST_P(TestSIMNodalConstraint, Vertex2DLR)
{
  if (GetParam() > 4)
    return;

  auto old = IFEM::getOptions().discretization;
  IFEM::getOptions().discretization = ASM::LRSpline;
  SIMNodalConstraint<SIM2D> s({2});
  std::stringstream str;
  str << "refdata/nodal_2D_V" << GetParam() << ".xinp";
  ASSERT_TRUE(s.read(str.str().c_str()));
  s.preprocess();

  const ASMu2D& pch = static_cast<const ASMu2D&>(*s.getPatch(1));
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == 1 && GetParam() == 1)
      check_mpc(pch.findMPC(i, 1), pch.getCorner(1,-1,1));
    else if (i == (size_t)pch.getCorner(1,-1,1) && GetParam() == 2)
      check_mpc(pch.findMPC(i, 1), 1);
    else if (i == (size_t)pch.getCorner(-1,1,1) && GetParam() == 3)
      check_mpc(pch.findMPC(i, 1), 1);
    else if (i == (size_t)pch.getCorner(1,1,1) && GetParam() == 4)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
  IFEM::getOptions().discretization = old;
}


TEST_P(TestSIMNodalConstraint, Edge2DLR)
{
  if (GetParam() > 4)
    return;

  auto old = IFEM::getOptions().discretization;
  IFEM::getOptions().discretization = ASM::LRSpline;
  SIMNodalConstraint<SIM2D> s({2});
  std::stringstream str;
  str << "refdata/nodal_2D_E" << GetParam() << ".xinp";
  ASSERT_TRUE(s.read(str.str().c_str()));
  s.preprocess();

  const ASMu2D& pch = static_cast<const ASMu2D&>(*s.getPatch(1));
  static const std::vector<LR::parameterEdge> E
        {LR::WEST, LR::EAST, LR::SOUTH, LR::NORTH};
  auto nodes = pch.getEdgeNodes(E[GetParam()-1], 1, 0);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (std::find(nodes.begin(), nodes.end(), i) != nodes.end() &&
        i != (size_t)pch.getCorner(1,-1,1) && GetParam() == 1)
      check_mpc(pch.findMPC(i, 1), pch.getCorner(1,-1,1));
    else if (std::find(nodes.begin(), nodes.end(), i) != nodes.end() && GetParam() == 3)
      check_mpc(pch.findMPC(i, 1), pch.getCorner(1,1,1));
    else if (std::find(nodes.begin(), nodes.end(), i) != nodes.end())
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
  IFEM::getOptions().discretization = old;
}


TEST_P(TestSIMNodalConstraint, Edge2DLRmx)
{
  if (GetParam() > 4)
    return;

  auto old = IFEM::getOptions().discretization;
  IFEM::getOptions().discretization = ASM::LRSpline;
  SIMNodalConstraint<SIM2D> s({2,2});
  std::stringstream str;
  str << "refdata/nodal_2D_E" << GetParam() << "_mixed.xinp";
  ASSERT_TRUE(s.read(str.str().c_str()));
  s.preprocess();

  const ASMu2D& pch = static_cast<const ASMu2D&>(*s.getPatch(1));
  static const std::vector<LR::parameterEdge> E
        {LR::WEST, LR::EAST, LR::SOUTH, LR::NORTH};
  auto nodes = pch.getEdgeNodes(E[GetParam()-1], 2, 0);
  size_t ofs = pch.getNoNodes(1);
  for (size_t i=1; i <= pch.getNoNodes(1); ++i) {
    ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
  for (size_t i=1; i <= pch.getNoNodes(2); ++i) {
    if (std::find(nodes.begin(), nodes.end(), i+ofs) != nodes.end() &&
        i+ofs != (size_t)pch.getCorner(1,-1,2) && GetParam() == 1)
      check_mpc(pch.findMPC(i+ofs, 1), pch.getCorner(1,-1,2));
    else if (std::find(nodes.begin(), nodes.end(), i+ofs) != nodes.end() && GetParam() == 3)
      check_mpc(pch.findMPC(i+ofs, 1), pch.getCorner(1,1,2));
    else if (std::find(nodes.begin(), nodes.end(), i+ofs) != nodes.end())
      check_mpc(pch.findMPC(i+ofs, 1), ofs+1);
    else
      ASSERT_TRUE(pch.findMPC(i+ofs, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i+ofs, 2) == nullptr);
  }
  IFEM::getOptions().discretization = old;
}
#endif


TEST_P(TestSIMNodalConstraint, Vertex3D)
{
  if (GetParam() > 8)
    return;

  SIMNodalConstraint<SIM3D> s({2});
  std::stringstream str;
  str << "refdata/nodal_3D_V" << GetParam() << ".xinp";
  ASSERT_TRUE(s.read(str.str().c_str()));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == 1 && GetParam() == 1)
      check_mpc(pch.findMPC(i, 1), n1);
    else if ((i == (size_t)n1 && GetParam() == 2)
          || (i == (size_t)n1*(n2-1)+1 && GetParam() == 3)
          || (i == (size_t)n1*n2 && GetParam() == 4)
          || (i == (size_t)1+n1*n2*(n3-1) && GetParam() == 5)
          || (i == (size_t)n1+n1*n2*(n3-1) && GetParam() == 6)
          || (i == (size_t)1+n1*(n2-1)+n1*n2*(n3-1) && GetParam() == 7)
          || (i == pch.getNoNodes() && GetParam() == 8))
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST_P(TestSIMNodalConstraint, Edge3D)
{
  SIMNodalConstraint<SIM3D> s({2});
  std::stringstream str;
  str << "refdata/nodal_3D_E" << GetParam() << ".xinp";
  ASSERT_TRUE(s.read(str.str().c_str()));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i <= (size_t)n1 && GetParam() == 1)
      check_mpc(pch.findMPC(i, 1), n1*n2*n3);
    else if (i >= (size_t)n1*(n2-1)+1 && i <= (size_t)n1*n2 && GetParam() == 2)
      check_mpc(pch.findMPC(i, 1), 1);
    else if (i >= (size_t)n1*n2*(n3-1)+1 &&
             i <= (size_t)n1*n2*(n3-1)+n1 && GetParam() == 3)
      check_mpc(pch.findMPC(i, 1), 1);
    else if (i >= (size_t)n1*n2*(n3-1)+n1*(n2-1)+1 && GetParam() == 4)
      check_mpc(pch.findMPC(i, 1), 1);
    else if (i <= (size_t)n1*n2 && i % n2 == 1 && GetParam() == 5)
      check_mpc(pch.findMPC(i, 1), n1*n2*n3);
    else if (i <= (size_t)n1*n2 && i % n2 == 0 && GetParam() == 6)
      check_mpc(pch.findMPC(i, 1), 1);
    else if (i >= (size_t)n1*n2*(n3-1)+1 && i % n2 == 1 && GetParam() == 7)
      check_mpc(pch.findMPC(i, 1), 1);
    else if (i >= (size_t)n1*n2*(n3-1)+1 && i % n2 == 0 && GetParam() == 8)
      check_mpc(pch.findMPC(i, 1), 1);
    else if (i % (n1*n2) == 1 && GetParam() == 9)
      check_mpc(pch.findMPC(i, 1), n1*n2*n3);
    else if (i % (n1*n2) == (size_t)n1 && GetParam() == 10)
      check_mpc(pch.findMPC(i, 1), 1);
    else if (i % (n1*n2) == (size_t)n1*(n2-1)+1 && GetParam() == 11)
      check_mpc(pch.findMPC(i, 1), 1);
    else if (i % (n1*n2) == 0 && GetParam() == 12)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST_P(TestSIMNodalConstraint, Face3D)
{
  if (GetParam() > 6)
    return;

  SIMNodalConstraint<SIM3D> s({2});
  std::stringstream str;
  str << "refdata/nodal_3D_F" << GetParam() << ".xinp";
  ASSERT_TRUE(s.read(str.str().c_str()));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i % n1 == 1 && GetParam() == 1)
      check_mpc(pch.findMPC(i, 1), n1*n2*n3);
    else if (i % n1 == 0 && GetParam() == 2)
      check_mpc(pch.findMPC(i, 1), 1);
    else if (i % (n1*n2) > 0 && i % (n1*n2) <= (size_t)n1 && GetParam() == 3)
      check_mpc(pch.findMPC(i, 1), n1*n2*n3);
    else if ((i % (n1*n2) >= (size_t)n1*(n2-1)+1 || i % (n1*n2) == 0) && GetParam() == 4)
      check_mpc(pch.findMPC(i, 1), 1);
    else if (i <= (size_t)n1*n2 && GetParam() == 5)
      check_mpc(pch.findMPC(i, 1), n1*n2*n3);
    else if (i >= (size_t)n1*n2*(n3-1)+1 && GetParam() == 6)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DmxE3)
{
  SIMNodalConstraint<SIM3D> s({2,2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_E3_mixed.xinp"));
  s.preprocess();

  const ASMs3Dmx& pch = static_cast<const ASMs3Dmx&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3,2);
  size_t ofs = pch.getNoNodes(1);
  for (size_t i=1; i <= pch.getNoNodes(1); ++i) {
    ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
  for (size_t i=1; i <= pch.getNoNodes(2); ++i) {
    if (i >= (size_t)n1*n2*(n3-1)+1 && i <= (size_t)n1*n2*(n3-1)+n1)
      check_mpc(pch.findMPC(i+ofs, 1), 1+ofs);
    else
      ASSERT_TRUE(pch.findMPC(i+ofs, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i+ofs, 2) == nullptr);
  }
}

const std::vector<int> tests = {1,2,3,4,5,6,7,8,9,10,11,12};
INSTANTIATE_TEST_CASE_P(TestSIMNodalConstraint, TestSIMNodalConstraint,
                        testing::ValuesIn(tests));
