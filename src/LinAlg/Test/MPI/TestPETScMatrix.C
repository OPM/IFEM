//==============================================================================
//!
//! \file TestPETScMatrix.C
//!
//! \date Apr 27 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for parallel PETSc matrices
//!
//==============================================================================

#include "PETScMatrix.h"
#include "SIM2D.h"
#include "SAM.h"
#include "readIntVec.h"

#include "Catch2Support.h"


TEST_CASE("TestPETScMatrix.AssembleMPI")
{
  SIM2D sim(1);
  sim.read("src/LinAlg/Test/refdata/petsc_test.xinp");
  sim.opt.solver = LinAlg::PETSC;
  REQUIRE(sim.preprocess());
  REQUIRE(sim.initSystem(sim.opt.solver));

  PETScMatrix* myMat = dynamic_cast<PETScMatrix*>(sim.getLHSmatrix());
  REQUIRE(myMat != nullptr);

  Matrix stencil(4,4);
  stencil.diag(1.0);

  for (int iel = 1; iel <= sim.getSAM()->getNoElms(); ++iel)
    myMat->assemble(stencil, *sim.getSAM(), iel);

  myMat->endAssembly();

  // now inspect the matrix
  Mat& mat = myMat->getMatrix();

  PetscInt r, c;
  MatGetLocalSize(mat, &r, &c);

  Vec vec;
  VecCreate(PETSC_COMM_WORLD, &vec);
  VecSetFromOptions(vec);
  VecSetType(vec, VECMPI);
  VecSetSizes(vec, r, PETSC_DECIDE);
  MatGetDiagonal(mat, vec);
  PetscScalar* a;
  VecGetArray(vec, &a);
  IntVec v = readIntVector("src/LinAlg/Test/refdata/petsc_matrix_diagonal.ref");
  VecGetOwnershipRange(vec, &r, &c);

  // check that we have the correct diagonal values
  for (int i = r; i < c; ++i)
    REQUIRE_THAT(v[i], WithinRel(a[i-r]));

  VecRestoreArray(vec, &a);
  VecDestroy(&vec);

  // check that no values outside the diagonal are != 0
  for (; r < c; ++r) {
    const PetscInt* cols;
    PetscInt ncols;
    const PetscScalar* vals;
    MatGetRow(mat, r, &ncols, &cols, &vals);
    for (PetscInt i = 0; i < ncols; ++i)
      if (cols[i] != r)
        REQUIRE_THAT(vals[i], WithinAbs(0.0, 1e-14));
    MatRestoreRow(mat, r, &ncols, &cols, &vals);
  }
}


TEST_CASE("TestPETScMatrix.AssembleBasisBlocksMPI")
{
  SIM2D sim({1,1});
  sim.read("src/LinAlg/Test/refdata/petsc_test_blocks_basis.xinp");
  sim.opt.solver = LinAlg::PETSC;
  REQUIRE(sim.preprocess());
  REQUIRE(sim.initSystem(sim.opt.solver));

  PETScMatrix* myMat = dynamic_cast<PETScMatrix*>(sim.getLHSmatrix());
  REQUIRE(myMat != nullptr);

  Matrix stencil(13,13);
  stencil.diag(1.0);
  for (size_t i = 10; i <= 13; i++)
    stencil(i,i) = 2.0;

  for (int iel = 1; iel <= sim.getSAM()->getNoElms(); ++iel)
    myMat->assemble(stencil, *sim.getSAM(), iel);

  myMat->endAssembly();

  // now inspect the matrix blocks
  const std::vector<Mat>& mat = myMat->getBlockMatrices();
  const ProcessAdm& adm = sim.getProcessAdm();

  for (size_t b = 0; b < mat.size(); ++b) {
    PetscInt mr = adm.dd.getMaxEq(b/2+1)-adm.dd.getMinEq(b/2+1)+1;

    if (b == 0 || b == 3) { // diagonal blocks
      Vec vec;
      VecCreate(PETSC_COMM_WORLD, &vec);
      VecSetFromOptions(vec);
      VecSetType(vec, VECMPI);
      VecSetSizes(vec, mr, PETSC_DECIDE);
      MatGetDiagonal(mat[b], vec);
      PetscScalar* a;
      VecGetArray(vec, &a);
      // check that we have the correct diagonal values
      std::stringstream str;
      str << "src/LinAlg/Test/refdata/petsc_matrix_diagonal_basis_block" << b/2 + 1 << ".ref";
      IntVec v = readIntVector(str.str());
      PetscInt r,c;
      VecGetOwnershipRange(vec, &r, &c);
      for (int i = r; i < c; ++i)
        REQUIRE_THAT(v[i], WithinRel(a[i-r]));

      VecRestoreArray(vec, &a);
      VecDestroy(&vec);
    }

    // check that no values outside the diagonal are != 0
    for (PetscInt r = 0; r < mr; ++r) {
      const PetscInt* cols;
      PetscInt ncols;
      const PetscScalar* vals;
      MatGetRow(mat[b], r+adm.dd.getMinEq(b/2+1)-1, &ncols, &cols, &vals);
      for (PetscInt i = 0; i < ncols; ++i)
        if (cols[i] != r+adm.dd.getMinEq(b/2+1)-1)
          REQUIRE_THAT(vals[i], WithinAbs(0.0, 1e-14));
      MatRestoreRow(mat[b], r+adm.dd.getMinEq(b/2+1)-1, &ncols, &cols, &vals);
    }
  }
}


TEST_CASE("TestPETScMatrix.AssembleComponentBlocksMPI")
{
  SIM2D sim(2);
  sim.read("src/LinAlg/Test/refdata/petsc_test_blocks_components.xinp");
  sim.opt.solver = LinAlg::PETSC;
  REQUIRE(sim.preprocess());
  REQUIRE(sim.initSystem(sim.opt.solver));

  PETScMatrix* myMat = dynamic_cast<PETScMatrix*>(sim.getLHSmatrix());
  REQUIRE(myMat != nullptr);

  Matrix stencil(4*2,4*2);
  for (size_t i = 1; i<= 4; ++i) {
    stencil(2*i-1,2*i-1) = 1.0;
    stencil(2*i,2*i) = 2.0;
  }

  for (int iel = 1; iel <= sim.getSAM()->getNoElms(); ++iel)
    myMat->assemble(stencil, *sim.getSAM(), iel);

  myMat->endAssembly();

  // now inspect the matrix blocks
  const std::vector<Mat>& mat = myMat->getBlockMatrices();
  const ProcessAdm& adm = sim.getProcessAdm();

  for (size_t b = 0; b < mat.size(); ++b) {
    PetscInt mr = adm.dd.getMaxEq(b/2+1)-adm.dd.getMinEq(b/2+1)+1;

    if (b == 0 || b == 3) { // diagonal blocks
      Vec vec;
      VecCreate(PETSC_COMM_WORLD, &vec);
      VecSetFromOptions(vec);
      VecSetType(vec, VECMPI);
      VecSetSizes(vec, mr, PETSC_DECIDE);
      MatGetDiagonal(mat[b], vec);
      PetscScalar* a;
      VecGetArray(vec, &a);
      // check that we have the correct diagonal values
      std::stringstream str;
      str << "src/LinAlg/Test/refdata/petsc_matrix_diagonal_components_block" << b/2 + 1 << ".ref";
      IntVec v = readIntVector(str.str());
      PetscInt r,c;
      VecGetOwnershipRange(vec, &r, &c);
      for (int i = r; i < c; ++i)
        REQUIRE_THAT(v[i], WithinRel(a[i-r]));

      VecRestoreArray(vec, &a);
      VecDestroy(&vec);
    }

    // check that no values outside the diagonal are != 0
    for (PetscInt r = 0; r < mr; ++r) {
      const PetscInt* cols;
      PetscInt ncols;
      const PetscScalar* vals;
      MatGetRow(mat[b], r+adm.dd.getMinEq(b/2+1)-1, &ncols, &cols, &vals);
      for (PetscInt i = 0; i < ncols; ++i)
        if (cols[i] != r+adm.dd.getMinEq(b/2+1)-1)
          REQUIRE_THAT(vals[i], WithinAbs(0.0, 1e-14));
      MatRestoreRow(mat[b], r+adm.dd.getMinEq(b/2+1)-1, &ncols, &cols, &vals);
    }
  }
}
