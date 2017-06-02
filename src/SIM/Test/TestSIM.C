//==============================================================================
//!
//! \file TestSIM.C
//!
//! \date Apr 27 2016
//!
//! \author Arne Morten Kvarving / NTNU
//!
//! \brief Tests for various expected SIM behavior.
//!
//==============================================================================

#include "SIM2D.h"
#include "SIM3D.h"
#include "ASMs3D.h"
#include "ASMmxBase.h"
#include "IntegrandBase.h"

#include "gtest/gtest.h"


template<class Dim> class TestProjectSIM : public Dim
{
public:
  TestProjectSIM(const SIMinput::CharVec& nf) : Dim(nf)
  {
    Dim::myProblem = new TestProjectIntegrand(Dim::dimension);
    EXPECT_TRUE(this->createDefaultModel());
    EXPECT_TRUE(this->preprocess());
  }
  virtual ~TestProjectSIM() {}

  bool addMixedMADOF(unsigned char basis, unsigned char nndof)
  {
    return this->addMADOF(basis, nndof);
  }

private:
  class TestProjectIntegrand : public IntegrandBase
  {
  public:
    TestProjectIntegrand(int dim) : IntegrandBase(dim) {}

    using IntegrandBase::evalSol;
    virtual bool evalSol(Vector& s, const FiniteElement&, const Vec3& X,
                         const std::vector<int>&) const
    {
      s.resize(1);
      s(1) = X.sum();
      return true;
    }

    virtual size_t getNoFields(int) const { return 1; }
  };
};


TEST(TestSIM2D, UniqueBoundaryNodes)
{
  const char* boundary_nodes = "<geometry>"
    "<patchfile>src/ASM/Test/refdata/square-4-orient0.g2</patchfile>"
    "<topology>"
    "  <connection master='1' medge='4' slave='2' sedge='3'/>"
    "  <connection master='1' medge='2' slave='3' sedge='1'/>"
    "  <connection master='2' medge='2' slave='4' sedge='1'/>"
    "  <connection master='3' medge='4' slave='4' sedge='3'/>"
    "</topology>"
    "<topologysets>"
    "  <set name='dir' type='edge'>"
    "    <item patch='1'>3</item>"
    "    <item patch='3'>3</item>"
    "  </set>"
    "</topologysets>"
    "</geometry>";

  SIM2D sim(1);
  ASSERT_TRUE(sim.loadXML(boundary_nodes));
  ASSERT_TRUE(sim.createFEMmodel());

  std::vector<int> vec;
  sim.getBoundaryNodes(sim.getUniquePropertyCode("dir",0),vec);

  std::sort(vec.begin(), vec.end());
  ASSERT_TRUE(std::unique(vec.begin(), vec.end()) == vec.end());
}


TEST(TestSIM2D, Periodic)
{
  const char* geometry = "<geometry>"
    "<patchfile>src/SIM/Test/cylinder-shell.g2</patchfile>"
    "<refine patch='1' u='3' v='3'/>"
    "<periodic patch='1' dir='2'/>"
    "</geometry>";

  SIM2D sim(1);
  ASSERT_TRUE(sim.loadXML(geometry));
  ASSERT_TRUE(sim.createFEMmodel());
}


TEST(TestSIM2D, ProjectSolution)
{
  TestProjectSIM<SIM2D> sim({1});

  Vector ssol, psol(sim.getNoDOFs());
  ASSERT_TRUE(sim.project(ssol,psol));

  size_t n = 1;
  for (size_t j = 0; j < 2; ++j)
    for (size_t i = 0; i < 2; ++i)
      EXPECT_FLOAT_EQ(ssol(n++), i + j);
}


TEST(TestSIM2D, ProjectSolutionMixed)
{
  TestProjectSIM<SIM2D> sim({1,1});

  Vector ssol, psol(sim.getNoDOFs());
  ASSERT_TRUE(sim.project(ssol,psol));

  size_t n = 1;
  for (size_t j = 0; j < 3; ++j)
    for (size_t i = 0; i < 3; ++i)
      EXPECT_FLOAT_EQ(ssol(n++), i/2.0 + j/2.0);
}


TEST(TestSIM3D, ProjectSolution)
{
  TestProjectSIM<SIM3D> sim({1});

  Vector ssol, psol(sim.getNoDOFs());
  ASSERT_TRUE(sim.project(ssol,psol));

  size_t n = 1;
  for (size_t k = 0; k < 2; ++k)
    for (size_t j = 0; j < 2; ++j)
      for (size_t i = 0; i < 2; ++i)
        EXPECT_FLOAT_EQ(ssol(n++), i + j + k);
}


TEST(TestSIM3D, ProjectSolutionMixed)
{
  TestProjectSIM<SIM3D> sim({1,1});

  Vector ssol, psol(sim.getNoDOFs());
  ASSERT_TRUE(sim.project(ssol,psol));

  size_t n = 1;
  for (size_t k = 0; k < 3; ++k)
    for (size_t j = 0; j < 3; ++j)
      for (size_t i = 0; i < 3; ++i)
        EXPECT_FLOAT_EQ(ssol(n++), i/2.0 + j/2.0 + k/2.0);
}


TEST(TestSIM2D, InjectPatchSolution)
{
  ASMmxBase::Type = ASMmxBase::REDUCED_CONT_RAISE_BASIS1;
  ASMmxBase::geoBasis = 2;
  TestProjectSIM<SIM2D> sim({1,1});
  ASMbase* pch = sim.getPatch(1);
  size_t n1 = sim.getNoNodes(1);
  size_t n2 = sim.getNoNodes(2);

  Vector sol(2*n1 + n2);
  Vector lsol(2*n1);
  size_t i, ofs;
  for (i = 0; i < n1; i++)
    lsol[2*i] = lsol[2*i+1] = i+1;

  ASSERT_TRUE(sim.addMixedMADOF(1, 2));
  sim.injectPatchSolution(sol, lsol, pch, 2, 1);
  for (i = ofs = 0; i < n1; i++) {
    EXPECT_FLOAT_EQ(sol[ofs++], i+1);
    EXPECT_FLOAT_EQ(sol[ofs++], i+1);
  }
  for (i = 0; i < n2; i++)
    EXPECT_FLOAT_EQ(sol[ofs++], 0);

  Vector sol2(n1 + 2*n2);
  Vector lsol2(2*n2);
  for (i = 0; i < n2; i++)
    lsol2[2*i] = lsol2[2*i+1] = i+1;

  ASSERT_TRUE(sim.addMixedMADOF(2, 2));
  sim.injectPatchSolution(sol2, lsol2, pch, 2, 2);
  for (i = ofs = 0; i < n1; i++)
    EXPECT_FLOAT_EQ(sol2[ofs++], 0);
  for (i = 0; i < n2; i++) {
    EXPECT_FLOAT_EQ(sol2[ofs++], i+1);
    EXPECT_FLOAT_EQ(sol2[ofs++], i+1);
  }
}


TEST(TestSIM3D, Periodic)
{
  ASMmxBase::Type = ASMmxBase::REDUCED_CONT_RAISE_BASIS1;

  SIM3D sim({3,1});
  ASSERT_TRUE(sim.createDefaultModel());

  ASMs3D* pch = dynamic_cast<ASMs3D*>(sim.getPatch(1));
  ASSERT_TRUE(pch != nullptr);
  ASSERT_TRUE(pch->uniformRefine(0,1));
  ASSERT_TRUE(pch->uniformRefine(1,1));
  ASSERT_TRUE(pch->uniformRefine(2,1));
  ASSERT_TRUE(sim.createFEMmodel());

  pch->closeBoundaries(1,0,1);
  pch->closeBoundaries(3,0,1);

  ASSERT_TRUE(sim.preprocess());
}
