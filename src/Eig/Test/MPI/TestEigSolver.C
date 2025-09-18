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

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinRel;


TEST_CASE("TestEigSolver.ArpackPETScMPI")
{
  const int param = GENERATE(1,2,3,4);

  SECTION("Mode " + std::to_string(param)) {
    LinSolParams spar;
    tinyxml2::XMLDocument doc;
    REQUIRE(doc.Parse(R"(<linearsolver>
                             <type>preonly</type>
                             <pc>lu</pc>
                             <package>mumps</package>
                          </linearsolver>)") == tinyxml2::XML_SUCCESS);
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
    eig::solve(&A, B.get(), eigs, eigVec, 3, 4, param);

    for (size_t i = 1; i <= 3; ++i)
      REQUIRE_THAT(eigs(i), WithinRel(i + (param == 3 ? 1.0 : 0.0)));
  }
}


#ifdef HAS_SLEPC
TEST_CASE("TestEigSolver.SLEPcMPI")
{
  const int param = GENERATE(1,2,3,4);

  SECTION("Mode " + std::to_string(param)) {
    LinSolParams spar;
    tinyxml2::XMLDocument doc;
    REQUIRE(doc.Parse(R"(<linearsolver>
                             <type>preonly</type>
                             <pc>lu</pc>
                             <package>mumps</package>
                          </linearsolver>)") == tinyxml2::XML_SUCCESS);
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
    eig::solve(&A, B.get(), eigs, eigVec, 3, 4, param + 10);
    REQUIRE(eigs.size() == 3);

    for (size_t i = 1; i <= 3; ++i)
      REQUIRE_THAT(eigs(i), WithinRel(static_cast<Real>(i)));
  }
}
#endif
