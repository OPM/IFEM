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
#include "ASMbase.h"
#include "ASMs2D.h"
#ifdef HAS_LRSPLINE
#include "ASMu2D.h"
#endif
#include "ASMs2Dmx.h"
#include "ASMs3D.h"
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

TEST(TestSIMNodalConstraint, ASMs1DV1)
{
  SIMNodalConstraint<SIM1D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_1D_V1.xinp"));
  s.preprocess();

  const ASMbase& pch = *s.getPatch(1);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == 1)
      check_mpc(pch.findMPC(i,1), pch.getNoNodes());
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs1DV2)
{
  SIMNodalConstraint<SIM1D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_1D_V2.xinp"));
  s.preprocess();

  const ASMbase& pch = *s.getPatch(1);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == pch.getNoNodes())
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs2DV1)
{
  SIMNodalConstraint<SIM2D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_2D_V1.xinp"));
  s.preprocess();

  const ASMs2D& pch = static_cast<const ASMs2D&>(*s.getPatch(1));
  int n1, n2;
  pch.getSize(n1,n2);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == 1)
      check_mpc(pch.findMPC(i, 1), n1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs2DV2)
{
  SIMNodalConstraint<SIM2D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_2D_V2.xinp"));
  s.preprocess();

  const ASMs2D& pch = static_cast<const ASMs2D&>(*s.getPatch(1));
  int n1, n2;
  pch.getSize(n1,n2);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == (size_t)n1)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs2DV3)
{
  SIMNodalConstraint<SIM2D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_2D_V3.xinp"));
  s.preprocess();

  const ASMs2D& pch = static_cast<const ASMs2D&>(*s.getPatch(1));
  int n1, n2;
  pch.getSize(n1,n2);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == (size_t)n1*(n2-1)+1)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs2DV4)
{
  SIMNodalConstraint<SIM2D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_2D_V4.xinp"));
  s.preprocess();

  const ASMs2D& pch = static_cast<const ASMs2D&>(*s.getPatch(1));
  int n1, n2;
  pch.getSize(n1,n2);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == pch.getNoNodes())
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs2DE1)
{
  SIMNodalConstraint<SIM2D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_2D_E1.xinp"));
  s.preprocess();

  const ASMs2D& pch = static_cast<const ASMs2D&>(*s.getPatch(1));
  int n1, n2;
  pch.getSize(n1,n2);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i % n1 == 1)
      check_mpc(pch.findMPC(i, 1), n1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}

TEST(TestSIMNodalConstraint, ASMs2DE2)
{
  SIMNodalConstraint<SIM2D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_2D_E2.xinp"));
  s.preprocess();

  const ASMs2D& pch = static_cast<const ASMs2D&>(*s.getPatch(1));
  int n1, n2;
  pch.getSize(n1,n2);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i % n1 == 0)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs2DE3)
{
  SIMNodalConstraint<SIM2D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_2D_E3.xinp"));
  s.preprocess();

  const ASMs2D& pch = static_cast<const ASMs2D&>(*s.getPatch(1));
  int n1, n2;
  pch.getSize(n1,n2);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i <= (size_t)n1)
      check_mpc(pch.findMPC(i, 1), n1*n2);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs2DE4)
{
  SIMNodalConstraint<SIM2D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_2D_E4.xinp"));
  s.preprocess();

  const ASMs2D& pch = static_cast<const ASMs2D&>(*s.getPatch(1));
  int n1, n2;
  pch.getSize(n1,n2);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i > (size_t)n1*(n2-1))
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}

#ifdef HAS_LRSPLINE
TEST(TestSIMNodalConstraint, ASMu2DV1)
{
  auto old = IFEM::getOptions().discretization;
  IFEM::getOptions().discretization = ASM::LRSpline;
  SIMNodalConstraint<SIM2D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_2D_V1.xinp"));
  s.preprocess();

  const ASMu2D& pch = static_cast<const ASMu2D&>(*s.getPatch(1));
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == 1)
      check_mpc(pch.findMPC(i, 1), pch.getCorner(1,-1,1));
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
  IFEM::getOptions().discretization = old;
}


TEST(TestSIMNodalConstraint, ASMu2DV2)
{
  auto old = IFEM::getOptions().discretization;
  IFEM::getOptions().discretization = ASM::LRSpline;
  SIMNodalConstraint<SIM2D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_2D_V2.xinp"));
  s.preprocess();

  const ASMu2D& pch = static_cast<const ASMu2D&>(*s.getPatch(1));
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == (size_t)pch.getCorner(1,-1,1))
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
  IFEM::getOptions().discretization = old;
}


TEST(TestSIMNodalConstraint, ASMu2DV3)
{
  auto old = IFEM::getOptions().discretization;
  IFEM::getOptions().discretization = ASM::LRSpline;
  SIMNodalConstraint<SIM2D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_2D_V3.xinp"));
  s.preprocess();

  const ASMu2D& pch = static_cast<const ASMu2D&>(*s.getPatch(1));
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == (size_t)pch.getCorner(-1,1,1))
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
  IFEM::getOptions().discretization = old;
}


TEST(TestSIMNodalConstraint, ASMu2DV4)
{
  auto old = IFEM::getOptions().discretization;
  IFEM::getOptions().discretization = ASM::LRSpline;
  SIMNodalConstraint<SIM2D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_2D_V4.xinp"));
  s.preprocess();

  const ASMu2D& pch = static_cast<const ASMu2D&>(*s.getPatch(1));
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == (size_t)pch.getCorner(1,1,1))
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
  IFEM::getOptions().discretization = old;
}


TEST(TestSIMNodalConstraint, ASMu2DE1)
{
  auto old = IFEM::getOptions().discretization;
  IFEM::getOptions().discretization = ASM::LRSpline;
  SIMNodalConstraint<SIM2D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_2D_E1.xinp"));
  s.preprocess();

  const ASMu2D& pch = static_cast<const ASMu2D&>(*s.getPatch(1));
  auto nodes = pch.getEdgeNodes(-1, 1);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (std::find(nodes.begin(), nodes.end(), i) != nodes.end() && i != (size_t)pch.getCorner(1,-1,1))
      check_mpc(pch.findMPC(i, 1), pch.getCorner(1,-1,1));
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
  IFEM::getOptions().discretization = old;
}

TEST(TestSIMNodalConstraint, ASMu2DE2)
{
  auto old = IFEM::getOptions().discretization;
  IFEM::getOptions().discretization = ASM::LRSpline;
  SIMNodalConstraint<SIM2D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_2D_E2.xinp"));
  s.preprocess();

  const ASMu2D& pch = static_cast<const ASMu2D&>(*s.getPatch(1));
  auto nodes = pch.getEdgeNodes(1, 1);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (std::find(nodes.begin(), nodes.end(), i) != nodes.end())
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
  IFEM::getOptions().discretization = old;
}


TEST(TestSIMNodalConstraint, ASMu2DE3)
{
  auto old = IFEM::getOptions().discretization;
  IFEM::getOptions().discretization = ASM::LRSpline;
  SIMNodalConstraint<SIM2D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_2D_E3.xinp"));
  s.preprocess();

  const ASMu2D& pch = static_cast<const ASMu2D&>(*s.getPatch(1));
  auto nodes = pch.getEdgeNodes(-2, 1);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (std::find(nodes.begin(), nodes.end(), i) != nodes.end())
      check_mpc(pch.findMPC(i, 1), pch.getCorner(1,1,1));
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
  IFEM::getOptions().discretization = old;
}


TEST(TestSIMNodalConstraint, ASMu2DE4)
{
  auto old = IFEM::getOptions().discretization;
  IFEM::getOptions().discretization = ASM::LRSpline;
  SIMNodalConstraint<SIM2D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_2D_E4.xinp"));
  s.preprocess();

  const ASMu2D& pch = static_cast<const ASMu2D&>(*s.getPatch(1));
  auto nodes = pch.getEdgeNodes(2, 1);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (std::find(nodes.begin(), nodes.end(), i) != nodes.end())
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
  IFEM::getOptions().discretization = old;
}


TEST(TestSIMNodalConstraint, ASMu2DmxE3)
{
  auto old = IFEM::getOptions().discretization;
  IFEM::getOptions().discretization = ASM::LRSpline;
  SIMNodalConstraint<SIM2D> s({2,2});
  ASSERT_TRUE(s.read("refdata/nodal_2D_E3_mixed.xinp"));
  s.preprocess();

  const ASMu2D& pch = static_cast<const ASMu2D&>(*s.getPatch(1));
  auto nodes = pch.getEdgeNodes(-2, 2);
  size_t ofs = pch.getNoNodes(1);
  for (size_t i=1; i <= pch.getNoNodes(1); ++i) {
    ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
  for (size_t i=1; i <= pch.getNoNodes(2); ++i) {
    if (std::find(nodes.begin(), nodes.end(), i+ofs) != nodes.end())
      check_mpc(pch.findMPC(i+ofs, 1), ofs+pch.getCorner(1,1,2));
    else
      ASSERT_TRUE(pch.findMPC(i+ofs, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i+ofs, 2) == nullptr);
  }
  IFEM::getOptions().discretization = old;
}
#endif


TEST(TestSIMNodalConstraint, ASMs3DV1)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_V1.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == 1)
      check_mpc(pch.findMPC(i, 1), n1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DV2)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_V2.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == (size_t)n1)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DV3)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_V3.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == (size_t)n1*(n2-1)+1)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DV4)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_V4.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == (size_t)n1*n2)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}

TEST(TestSIMNodalConstraint, ASMs3DV5)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_V5.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == (size_t)1+n1*n2*(n3-1))
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DV6)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_V6.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == (size_t)n1+n1*n2*(n3-1))
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DV7)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_V7.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == (size_t)1+n1*(n2-1)+n1*n2*(n3-1))
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DV8)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_V8.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i == pch.getNoNodes())
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DE1)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_E1.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i <= (size_t)n1)
      check_mpc(pch.findMPC(i, 1), n1*n2*n3);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DE2)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_E2.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i >= (size_t)n1*(n2-1)+1 && i <= (size_t)n1*n2)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DE3)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_E3.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i >= (size_t)n1*n2*(n3-1)+1 && i <= (size_t)n1*n2*(n3-1)+n1)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DE4)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_E4.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i >= (size_t)n1*n2*(n3-1)+n1*(n2-1)+1)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DE5)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_E5.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i <= (size_t)n1*n2 && i % n2 == 1)
      check_mpc(pch.findMPC(i, 1), n1*n2*n3);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}

TEST(TestSIMNodalConstraint, ASMs3DE6)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_E6.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i <= (size_t)n1*n2 && i % n2 == 0)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DE7)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_E7.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i >= (size_t)n1*n2*(n3-1)+1 && i % n2 == 1)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DE8)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_E8.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i >= (size_t)n1*n2*(n3-1)+1 && i % n2 == 0)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}

TEST(TestSIMNodalConstraint, ASMs3DE9)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_E9.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i % (n1*n2) == 1)
      check_mpc(pch.findMPC(i, 1), n1*n2*n3);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}

TEST(TestSIMNodalConstraint, ASMs3DE10)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_E10.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i % (n1*n2) == (size_t)n1)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}

TEST(TestSIMNodalConstraint, ASMs3DE11)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_E11.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i % (n1*n2) == (size_t)n1*(n2-1)+1)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DE12)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_E12.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i % (n1*n2) == 0)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DF1)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_F1.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i % n1 == 1)
      check_mpc(pch.findMPC(i, 1), n1*n2*n3);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DF2)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_F2.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i % n1 == 0)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DF3)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_F3.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i % (n1*n2) > 0 && i % (n1*n2) <= (size_t)n1)
      check_mpc(pch.findMPC(i, 1), n1*n2*n3);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DF4)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_F4.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i % (n1*n2) >= (size_t)n1*(n2-1)+1 || i % (n1*n2) == 0)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DF5)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_F5.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i <= (size_t)n1*n2)
      check_mpc(pch.findMPC(i, 1), n1*n2*n3);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs3DF6)
{
  SIMNodalConstraint<SIM3D> s({2});
  ASSERT_TRUE(s.read("refdata/nodal_3D_F6.xinp"));
  s.preprocess();

  const ASMs3D& pch = static_cast<const ASMs3D&>(*s.getPatch(1));
  int n1, n2, n3;
  pch.getSize(n1,n2,n3);
  for (size_t i=1; i <= pch.getNoNodes(); ++i) {
    if (i >= (size_t)n1*n2*(n3-1)+1)
      check_mpc(pch.findMPC(i, 1), 1);
    else
      ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
}


TEST(TestSIMNodalConstraint, ASMs2DmxE3)
{
  SIMNodalConstraint<SIM2D> s({2,2});
  ASSERT_TRUE(s.read("refdata/nodal_2D_E3_mixed.xinp"));
  s.preprocess();

  const ASMs2Dmx& pch = static_cast<const ASMs2Dmx&>(*s.getPatch(1));
  int n1, n2;
  pch.getSize(n1,n2,2);
  size_t ofs = pch.getNoNodes(1);
  for (size_t i=1; i <= pch.getNoNodes(1); ++i) {
    ASSERT_TRUE(pch.findMPC(i, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i, 2) == nullptr);
  }
  for (size_t i=1; i <= pch.getNoNodes(2); ++i) {
    if (i <= (size_t)n1)
      check_mpc(pch.findMPC(i+ofs, 1), ofs+n1*n2);
    else
      ASSERT_TRUE(pch.findMPC(i+ofs, 1) == nullptr);
    ASSERT_TRUE(pch.findMPC(i+ofs, 2) == nullptr);
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
