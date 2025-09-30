//==============================================================================
//!
//! \file TestEigSolver.C
//!
//! \date Jul 3 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for eigenvalue solvers.
//!
//==============================================================================

#include "EigSolver.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"

#ifdef HAS_PETSC
#include "ProcessAdm.h"
#include "PETScMatrix.h"
#include <tinyxml2.h>
#endif

#include "Catch2Support.h"

#include <string>


TEST_CASE("TestEigSolver.Lapack")
{
  DenseMatrix A, B;
  A.redim(3,3);
  B.redim(3,3);

  for (size_t i = 1; i <= 3; ++i) {
    A(i,i) = i;
    B(i,i) = 1.0;
  }

  Vector eigs(3);
  Matrix eigVec(3,3);
  eig::solve(&A, &B, eigs, eigVec, 3);

  for (size_t i = 1; i <= 3; ++i)
    REQUIRE(eigs(i) == i);

  for (size_t i = 1; i <= 3; ++i)
    for (size_t j = 1; j <= 3; ++j)
      REQUIRE(eigVec(i,j) == (i==j ? 1.0 : 0.0));
}


TEST_CASE("TestEigSolver.Arpack")
{
  const int param = GENERATE(1,2,3,4);
  SECTION("Mode " + std::to_string(param)) {
    SparseMatrix A(SparseMatrix::SUPERLU), B(SparseMatrix::SUPERLU);
    A.redim(4,4);
    B.redim(4,4);

    for (size_t i = 1; i <= 4; ++i) {
      A(i,i) = i;
      B(i,i) = 1.0;
    }

    Vector eigs;
    Matrix eigVec;
    eig::solve(&A, &B, eigs, eigVec, 3, 4, param);

    for (size_t i = 1; i <= 3; ++i)
      REQUIRE_THAT(eigs(i), WithinRel(i + (param == 3 ? 1.0 : 0.0)));
  }
}


#ifdef HAS_PETSC
TEST_CASE("TestEigSolver.ArpackPETSc")
{
  const int param = GENERATE(1,2,3,4);
  SECTION("Mode " + std::to_string(param)) {
    LinSolParams spar;
    tinyxml2::XMLDocument doc;
    REQUIRE(doc.Parse(R"(<linearsolver>
                             <type>preonly</type>
                             <pc>lu</pc>
                          </linearsolver>)") == tinyxml2::XML_SUCCESS);
    spar.read(doc.RootElement());
    ProcessAdm adm;
    PETScMatrix A(adm, spar), B(adm, spar);
    const IntMat elms {
        {0}, {1}, {2}, {3},
    };
    A.init(4);
    B.init(4);

    Matrix eA(1,1);
    Matrix eB(1,1);
    eB(1,1) = 1.0;

    for (int iel = 0; iel < 4; ++iel) {
      eA(1,1) = iel + 1;
      A.assemble(eA, elms[iel]);
      B.assemble(eB, elms[iel]);
    }

    A.endAssembly();
    B.endAssembly();

    Vector eigs;
    Matrix eigVec;
    eig::solve(&A, &B, eigs, eigVec, 3, 4, param);

    for (size_t i = 1; i <= 3; ++i)
      REQUIRE_THAT(eigs(i), WithinRel(i + (param == 3 ? 1.0 : 0.0)));
  }
}
#endif


#ifdef HAS_SLEPC
TEST_CASE("TestEigSolver.SLEPc")
{
  const int param = GENERATE(1,2,3,4);
  SECTION("Mode " + std::to_string(param)) {
    LinSolParams spar;
    tinyxml2::XMLDocument doc;
    REQUIRE(doc.Parse(R"(<linearsolver>
                             <type>preonly</type>
                             <pc>lu</pc>
                          </linearsolver>)") == tinyxml2::XML_SUCCESS);
    spar.read(doc.RootElement());
    ProcessAdm adm;
    PETScMatrix A(adm, spar), B(adm, spar);
    const IntMat elms {
        {0}, {1}, {2}, {3},
    };
    A.init(4);
    B.init(4);

    Matrix eA(1,1);
    Matrix eB(1,1);
    eB(1,1) = 1.0;

    for (int iel = 0; iel < 4; ++iel) {
      eA(1,1) = iel + 1;
      A.assemble(eA, elms[iel]);
      B.assemble(eB, elms[iel]);
    }

    A.endAssembly();
    B.endAssembly();

    Vector eigs;
    Matrix eigVec;
    eig::solve(&A, &B, eigs, eigVec, 3, 4, param + 10);
    REQUIRE(eigs.size() == 3);

    for (size_t i = 1; i <= 3; ++i)
      REQUIRE_THAT(eigs(i), WithinRel(static_cast<double>(i)));
  }
}
#endif
