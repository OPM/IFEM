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

#include "IntegrandBase.h"
#include "IFEM.h"
#include "SIM2D.h"
#include "SIM3D.h"

#include "gtest/gtest.h"
#include "tinyxml.h"


template<class Dim>
class TestProjectSIM : public Dim {
public:
  TestProjectSIM(const SIMbase::CharVec& nf) : Dim(nf)
  {
    Dim::myProblem = new TestProjectIntegrand(Dim::dimension);
  }

  bool addMixedMADOF(unsigned char basis, unsigned char nndof)
  {
    return this->addMADOF(basis, nndof);
  }
private:
  class TestProjectIntegrand : public IntegrandBase {
  public:
    TestProjectIntegrand(int dim) : IntegrandBase(dim) {}

    //! \brief Evaluates the secondary solution at a result point.
    //! \param[out] s The solution field values at current point
    //! \param[in] fe Finite element data at current point
    //! \param[in] X Cartesian coordinates of current point
    //! \param[in] MNPC Nodal point correspondance for the basis function values
    bool evalSol(Vector& s, const FiniteElement& fe, const Vec3& X,
                 const std::vector<int>& MNPC) const override
    {
      s.resize(1);
      s(1) = X[0] + X[1] + X[2];

      return true;
    }

    size_t getNoFields(int = 2) const override { return 1; }
  };
};


TEST(TestSIM, UniqueBoundaryNodes)
{
  SIM2D sim(1);
  ASSERT_TRUE(sim.read("src/SIM/Test/refdata/boundary_nodes.xinp"));
  sim.preprocess();

  int bcode = sim.getUniquePropertyCode("dir",0);
  std::vector<int> vec;
  sim.getBoundaryNodes(bcode, vec);

  std::sort(vec.begin(), vec.end());
  ASSERT_TRUE(std::unique(vec.begin(), vec.end()) == vec.end());
}


TEST(TestSIM2D, ProjectSolution)
{
  TestProjectSIM<SIM2D> sim({1});
  sim.createDefaultModel();
  ASSERT_TRUE(sim.preprocess());

  Matrix ssol;
  ASSERT_TRUE(sim.project(ssol, Vector(sim.getNoDOFs())));

  size_t n = 1;
  for (size_t j = 0; j < 2; ++j)
    for (size_t i = 0; i < 2; ++i)
      ASSERT_FLOAT_EQ(ssol(1, n++), i + j);
}


TEST(TestSIM2D, ProjectSolutionMixed)
{
  TestProjectSIM<SIM2D> sim({1,1});
  sim.createDefaultModel();
  ASSERT_TRUE(sim.preprocess());

  Matrix ssol;
  ASSERT_TRUE(sim.project(ssol, Vector(sim.getNoDOFs())));

  size_t n = 1;
  for (size_t j = 0; j < 3; ++j)
    for (size_t i = 0; i < 3; ++i)
      ASSERT_FLOAT_EQ(ssol(1, n++), i/2.0 + j/2.0);
}


TEST(TestSIM3D, ProjectSolution)
{
  TestProjectSIM<SIM3D> sim({1});
  sim.createDefaultModel();
  ASSERT_TRUE(sim.preprocess());

  Matrix ssol;
  ASSERT_TRUE(sim.project(ssol, Vector(sim.getNoDOFs())));

  size_t n = 1;
  for (size_t k = 0; k < 2; ++k)
    for (size_t j = 0; j < 2; ++j)
      for (size_t i = 0; i < 2; ++i)
        ASSERT_FLOAT_EQ(ssol(1, n++), i + j + k);
}


TEST(TestSIM3D, ProjectSolutionMixed)
{
  TestProjectSIM<SIM3D> sim({1,1});
  sim.createDefaultModel();
  ASSERT_TRUE(sim.preprocess());

  Matrix ssol;
  ASSERT_TRUE(sim.project(ssol, Vector(sim.getNoDOFs())));

  size_t n = 1;
  for (size_t k = 0; k < 3; ++k)
    for (size_t j = 0; j < 3; ++j)
      for (size_t i = 0; i < 3; ++i)
        ASSERT_FLOAT_EQ(ssol(1, n++), i/2.0 + j/2.0 + k/2.0);
}


TEST(TestSIM, InjectPatchSolution)
{
  TestProjectSIM<SIM2D> sim({1,1});
  sim.createDefaultModel();
  ASSERT_TRUE(sim.preprocess());

  Vector sol(2*sim.getNoNodes(true, 1) + sim.getNoNodes(true, 2));
  Vector lsol(2*sim.getNoNodes(true, 1));
  for (size_t i = 0; i < sim.getNoNodes(true,1); ++i)
    lsol[2*i] = lsol[2*i+1] = i+1;

  ASSERT_TRUE(sim.addMixedMADOF(1, 2));
  sim.injectPatchSolution(sol, lsol, 0, 2, 1);
  size_t ofs = 0;
  for (size_t i = 0; i < sim.getNoNodes(true,1); ++i, ofs += 2) {
    ASSERT_FLOAT_EQ(sol[ofs], i+1);
    ASSERT_FLOAT_EQ(sol[ofs+1], i+1);
  }
  for (size_t i = 0; i < sim.getNoNodes(true,2); ++i, ++ofs)
    ASSERT_FLOAT_EQ(sol[ofs], 0);

  ASSERT_TRUE(sim.addMixedMADOF(2, 2));
  Vector sol2(sim.getNoNodes(true, 1) + 2*sim.getNoNodes(true, 2));
  Vector lsol2(2*sim.getNoNodes(true, 2));
  for (size_t i = 0; i < sim.getNoNodes(true,2); ++i)
    lsol2[2*i] = lsol2[2*i+1] = i+1;

  sim.injectPatchSolution(sol2, lsol2, 0, 2, 2);
  ofs = 0;
  for (size_t i = 0; i < sim.getNoNodes(true,1); ++i, ++ofs)
    ASSERT_FLOAT_EQ(sol2[ofs], 0);

  for (size_t i = 0; i < sim.getNoNodes(true,2); ++i, ofs += 2) {
    ASSERT_FLOAT_EQ(sol2[ofs], i+1);
    ASSERT_FLOAT_EQ(sol2[ofs+1], i+1);
  }
}
