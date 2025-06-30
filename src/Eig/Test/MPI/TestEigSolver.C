//==============================================================================
//!
//! \file TestEigSolver.C
//!
//! \date Jun 25 2025
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for parallel eigenvalue solvers.
//!
//==============================================================================

#include "EigSolver.h"

#include "ProcessAdm.h"
#include "PETScMatrix.h"
#include <tinyxml2.h>

#include "IFEM.h"

#include "gtest/gtest.h"


class TestEigSolverMPI : public testing::Test,
                         public testing::WithParamInterface<int>
{
};



TEST_P(TestEigSolverMPI, ArpackPETSc)
{
  LinSolParams spar;
  tinyxml2::XMLDocument doc;
  EXPECT_EQ(doc.Parse(R"(<linearsolver>
                           <type>preonly</type>
                           <pc>lu</pc>
                           <package>mumps</package>
                        </linearsolver>)"), tinyxml2::XML_SUCCESS);
  spar.read(doc.RootElement());
  ProcessAdm adm(PETSC_COMM_WORLD);
  adm.dd.setElms({1}, ""); // just to flag partitioning
  PETScMatrix A(adm, spar);
  const IntMat elms {
      {0}, {1}, {2}, {3},
  };
  const IntMat neigh {
    {1},
    {0,2},
    {1,3},
    {2},
  };
  A.init(4, &elms, &neigh, nullptr);
  std::unique_ptr<PETScMatrix> B(static_cast<PETScMatrix*>(A.copy()));

  Matrix eA(1,1);
  Matrix eB(1,1);
  eB(1,1) = 1.0;

  for (int iel : A.getDD().getElms()) {
    eA(1,1) = iel + 1;
    A.assemble(eA, elms[iel]);
    B->assemble(eB, elms[iel]);
  }

  A.endAssembly();
  B->endAssembly();

  Vector eigs;
  Matrix eigVec;
  eig::solve(&A, B.get(), eigs, eigVec, 3, 4, GetParam());

  for (size_t i = 1; i <= 3; ++i)
    EXPECT_FLOAT_EQ(eigs(i), i + (GetParam() == 3 ? 1 : 0));
}


#ifdef HAS_SLEPC
TEST_P(TestEigSolverMPI, SLEPc)
{
  LinSolParams spar;
  tinyxml2::XMLDocument doc;
  EXPECT_EQ(doc.Parse(R"(<linearsolver>
                           <type>preonly</type>
                           <pc>lu</pc>
                           <package>mumps</package>
                        </linearsolver>)"), tinyxml2::XML_SUCCESS);
  spar.read(doc.RootElement());
  ProcessAdm adm(PETSC_COMM_WORLD);
  adm.dd.setElms({1}, ""); // just to flag partitioning
  PETScMatrix A(adm, spar);
  const IntMat elms {
      {0}, {1}, {2}, {3},
  };
  const IntMat neigh {
    {1},
    {0,2},
    {1,3},
    {2},
  };
  A.init(4, &elms, &neigh, nullptr);
  std::unique_ptr<PETScMatrix> B(static_cast<PETScMatrix*>(A.copy()));

  Matrix eA(1,1);
  Matrix eB(1,1);
  eB(1,1) = 1.0;

  for (int iel : A.getDD().getElms()) {
    eA(1,1) = iel + 1;
    A.assemble(eA, elms[iel]);
    B->assemble(eB, elms[iel]);
  }

  A.endAssembly();
  B->endAssembly();

  Vector eigs;
  Matrix eigVec;
  eig::solve(&A, B.get(), eigs, eigVec, 3, 4, GetParam() + 10);
  EXPECT_EQ(eigs.size(), 3);

  for (size_t i = 1; i <= 3; ++i)
    EXPECT_FLOAT_EQ(eigs(i), i);
}
#endif


const std::vector<int> modes = {1,2,3,4};
INSTANTIATE_TEST_SUITE_P(TestEigSolver, TestEigSolverMPI, testing::ValuesIn(modes));
