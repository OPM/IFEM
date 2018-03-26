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
#include "IntegrandBase.h"
#include "ASMmxBase.h"

#include "gtest/gtest.h"
#include "tinyxml.h"


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
      s(1) = X[0] + X[1] + X[2];
      return true;
    }

    virtual size_t getNoFields(int) const { return 1; }
  };
};


class DummyIntegrand : public IntegrandBase {};


TEST(TestSIM, UniqueBoundaryNodes)
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

  Matrix ssol;
  ASSERT_TRUE(sim.project(ssol, Vector(sim.getNoDOFs())));

  size_t n = 1;
  for (size_t j = 0; j < 2; ++j)
    for (size_t i = 0; i < 2; ++i)
      EXPECT_FLOAT_EQ(ssol(1, n++), i + j);
}


TEST(TestSIM2D, ProjectSolutionMixed)
{
  TestProjectSIM<SIM2D> sim({1,1});

  Matrix ssol;
  ASSERT_TRUE(sim.project(ssol, Vector(sim.getNoDOFs())));

  size_t n = 1;
  for (size_t j = 0; j < 3; ++j)
    for (size_t i = 0; i < 3; ++i)
      EXPECT_FLOAT_EQ(ssol(1, n++), i/2.0 + j/2.0);
}


TEST(TestSIM3D, ProjectSolution)
{
  TestProjectSIM<SIM3D> sim({1});

  Matrix ssol;
  ASSERT_TRUE(sim.project(ssol, Vector(sim.getNoDOFs())));

  size_t n = 1;
  for (size_t k = 0; k < 2; ++k)
    for (size_t j = 0; j < 2; ++j)
      for (size_t i = 0; i < 2; ++i)
        EXPECT_FLOAT_EQ(ssol(1, n++), i + j + k);
}


TEST(TestSIM3D, ProjectSolutionMixed)
{
  TestProjectSIM<SIM3D> sim({1,1});

  Matrix ssol;
  ASSERT_TRUE(sim.project(ssol, Vector(sim.getNoDOFs())));

  size_t n = 1;
  for (size_t k = 0; k < 3; ++k)
    for (size_t j = 0; j < 3; ++j)
      for (size_t i = 0; i < 3; ++i)
        EXPECT_FLOAT_EQ(ssol(1, n++), i/2.0 + j/2.0 + k/2.0);
}


TEST(TestSIM, InjectPatchSolution)
{
  ASMmxBase::Type = ASMmxBase::REDUCED_CONT_RAISE_BASIS1;
  ASMmxBase::elmBasis = 2;
  TestProjectSIM<SIM2D> sim({1,1});
  ASMbase* pch = sim.getPatch(1);

  Vector sol(2*sim.getNoNodes(1) + sim.getNoNodes(2));
  Vector lsol(2*sim.getNoNodes(1));
  size_t i, ofs;
  for (i = 0; i < sim.getNoNodes(1); i++)
    lsol[2*i] = lsol[2*i+1] = i+1;

  ASSERT_TRUE(sim.addMixedMADOF(1, 2));
  sim.injectPatchSolution(sol, lsol, pch, 2, 1);
  for (i = ofs = 0; i < sim.getNoNodes(1); i++, ofs += 2) {
    EXPECT_FLOAT_EQ(sol[ofs], i+1);
    EXPECT_FLOAT_EQ(sol[ofs+1], i+1);
  }
  for (i = 0; i < sim.getNoNodes(2); i++, ofs++)
    EXPECT_FLOAT_EQ(sol[ofs], 0);

  ASSERT_TRUE(sim.addMixedMADOF(2, 2));
  Vector sol2(sim.getNoNodes(1) + 2*sim.getNoNodes(2));
  Vector lsol2(2*sim.getNoNodes(2));
  for (i = 0; i < sim.getNoNodes(2); i++)
    lsol2[2*i] = lsol2[2*i+1] = i+1;

  sim.injectPatchSolution(sol2, lsol2, pch, 2, 2);
  for (i = ofs = 0; i < sim.getNoNodes(1); i++, ofs++)
    EXPECT_FLOAT_EQ(sol2[ofs], 0);

  for (i = 0; i < sim.getNoNodes(2); i++, ofs += 2) {
    EXPECT_FLOAT_EQ(sol2[ofs], i+1);
    EXPECT_FLOAT_EQ(sol2[ofs+1], i+1);
  }
}
