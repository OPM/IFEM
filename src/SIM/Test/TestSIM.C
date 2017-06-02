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
#include "ASMmxBase.h"

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


class DummyIntegrand : public IntegrandBase {};


TEST(TestSIM2D, UniqueBoundaryNodes)
{
  SIM2D sim(new DummyIntegrand(),1);
  ASSERT_TRUE(sim.read("src/SIM/Test/refdata/boundary_nodes.xinp"));
  ASSERT_TRUE(sim.preprocess());

  std::vector<int> vec;
  sim.getBoundaryNodes(sim.getUniquePropertyCode("dir",0),vec);

  std::sort(vec.begin(), vec.end());
  ASSERT_TRUE(std::unique(vec.begin(), vec.end()) == vec.end());
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
