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
  PETScMatrix A(adm, spar), B(adm, spar);
  A.redim(4,4);
  B.redim(4,4);

  for (size_t i = 1; i <= 4; ++i) {
    A(i,i) = i;
    B(i,i) = 1.0;
  }

  A.endAssembly();
  B.endAssembly();

  Vector eigs;
  Matrix eigVec;
  eig::solve(&A, &B, eigs, eigVec, 3, 4, GetParam());

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
  PETScMatrix A(adm, spar), B(adm, spar);
  A.redim(4,4);
  B.redim(4,4);

  for (size_t i = 1; i <= 4; ++i) {
    A(i,i) = i;
    B(i,i) = 1.0;
  }

  A.endAssembly();
  B.endAssembly();

  Vector eigs;
  Matrix eigVec;
  eig::solve(&A, &B, eigs, eigVec, 3, 4, GetParam() + 10);
  EXPECT_EQ(eigs.size(), 3);

  for (size_t i = 1; i <= 3; ++i)
    EXPECT_FLOAT_EQ(eigs(i), i);
}
#endif


const std::vector<int> modes = {1,2,3,4};
INSTANTIATE_TEST_SUITE_P(TestEigSolver, TestEigSolverMPI, testing::ValuesIn(modes));
