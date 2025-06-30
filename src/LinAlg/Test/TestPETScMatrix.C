//==============================================================================
//!
//! \file TestPETScMatrix.C
//!
//! \date Apr 27 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for PETSc matrices
//!
//==============================================================================

#include "ASMbase.h"
#include "PETScMatrix.h"
#include "SIM2D.h"
#include "ASMmxBase.h"
#include "SAM.h"
#include "readIntVec.h"

#include "gtest/gtest.h"

#include <memory>
#include <numeric>


TEST(TestPETScMatrix, Assemble)
{
  SIM2D sim(1);
  sim.read("src/LinAlg/Test/refdata/petsc_test.xinp");
  sim.opt.solver = LinAlg::PETSC;
  ASSERT_TRUE(sim.preprocess());
  ASSERT_TRUE(sim.initSystem(sim.opt.solver));

  PETScMatrix* myMat = dynamic_cast<PETScMatrix*>(sim.getLHSmatrix());
  ASSERT_TRUE(myMat != nullptr);

  myMat->init();

  Matrix stencil(4,4);
  stencil.diag(1.0);

  for (int iel = 1; iel <= sim.getSAM()->getNoElms(); ++iel)
    myMat->assemble(stencil, *sim.getSAM(), iel);

  myMat->endAssembly();

  // now inspect the matrix
  Mat& mat = myMat->getMatrix();

  PetscInt mr, mc;
  MatGetSize(mat, &mr, &mc);

  Vec vec;
  VecCreate(PETSC_COMM_SELF, &vec);
  VecSetFromOptions(vec);
  VecSetType(vec, VECSEQ);
  VecSetSizes(vec, mr, PETSC_DECIDE);
  MatGetDiagonal(mat, vec);
  PetscScalar* a;
  VecGetArray(vec, &a);
  // check that we have the correct diagonal values
  IntVec v = readIntVector("src/LinAlg/Test/refdata/petsc_matrix_diagonal.ref");
  ASSERT_EQ(mr, (int)v.size());
  for (int i = 0; i < mr; ++i)
    EXPECT_FLOAT_EQ(v[i], a[i]);

  VecRestoreArray(vec, &a);
  VecDestroy(&vec);

  // check that no values outside the diagonal are != 0
  for (PetscInt r = 0; r < mr; ++r) {
    const PetscInt* cols;
    PetscInt ncols;
    const PetscScalar* vals;
    MatGetRow(mat, r, &ncols, &cols, &vals);
    for (PetscInt i = 0; i < ncols; ++i)
      if (cols[i] != r)
        EXPECT_FLOAT_EQ(vals[i], 0.0);
    MatRestoreRow(mat, r, &ncols, &cols, &vals);
  }
}


TEST(TestPETScMatrix, AssembleBasisBlocks)
{
  ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS1;
  ASMmxBase::itgBasis = 2;
  SIM2D sim({1,1});
  sim.read("src/LinAlg/Test/refdata/petsc_test_blocks_basis.xinp");
  sim.opt.solver = LinAlg::PETSC;
  ASSERT_TRUE(sim.preprocess());
  ASSERT_TRUE(sim.initSystem(sim.opt.solver));

  PETScMatrix* myMat = dynamic_cast<PETScMatrix*>(sim.getLHSmatrix());
  ASSERT_TRUE(myMat != nullptr);

  Matrix stencil(13,13);
  stencil.diag(1.0);
  for (size_t i = 10; i <= 13; i++)
    stencil(i,i) = 2.0;

  for (int iel = 1; iel <= sim.getSAM()->getNoElms(); ++iel)
    myMat->assemble(stencil, *sim.getSAM(), iel);

  myMat->endAssembly();

  // now inspect the matrix blocks
  const std::vector<Mat>& mat = myMat->getBlockMatrices();

  for (size_t b = 0; b < mat.size(); ++b) {
    PetscInt mr, mc;
    MatGetSize(mat[b], &mr, &mc);

    if (b == 0 || b == 3) { // diagonal blocks
      Vec vec;
      VecCreate(PETSC_COMM_SELF, &vec);
      VecSetFromOptions(vec);
      VecSetType(vec, VECSEQ);
      VecSetSizes(vec, mr, PETSC_DECIDE);
      MatGetDiagonal(mat[b], vec);
      PetscScalar* a;
      VecGetArray(vec, &a);
      // check that we have the correct diagonal values
      std::stringstream str;
      str << "src/LinAlg/Test/refdata/petsc_matrix_diagonal_basis_block" << b/2 + 1 << ".ref";
      IntVec v = readIntVector(str.str());
      ASSERT_EQ(mr, (int)v.size());
      for (int i = 0; i < mr; ++i)
        EXPECT_FLOAT_EQ(v[i], a[i]);

      VecRestoreArray(vec, &a);
      VecDestroy(&vec);
    }

    // check that no values outside the diagonal are != 0
    for (PetscInt r = 0; r < mr; ++r) {
      const PetscInt* cols;
      PetscInt ncols;
      const PetscScalar* vals;
      MatGetRow(mat[b], r, &ncols, &cols, &vals);
      for (PetscInt i = 0; i < ncols; ++i)
        if (cols[i] != r)
          EXPECT_FLOAT_EQ(vals[i], 0.0);
      MatRestoreRow(mat[b], r, &ncols, &cols, &vals);
    }
  }
}


TEST(TestPETScMatrix, AssembleComponentBlocks)
{
  SIM2D sim(2);
  sim.read("src/LinAlg/Test/refdata/petsc_test_blocks_components.xinp");
  sim.opt.solver = LinAlg::PETSC;
  ASSERT_TRUE(sim.preprocess());
  ASSERT_TRUE(sim.initSystem(sim.opt.solver));

  PETScMatrix* myMat = dynamic_cast<PETScMatrix*>(sim.getLHSmatrix());
  ASSERT_TRUE(myMat != nullptr);

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

  for (size_t b = 0; b < mat.size(); ++b) {
    PetscInt mr, mc;
    MatGetSize(mat[b], &mr, &mc);

    if (b == 0 || b == 3) { // diagonal blocks
      Vec vec;
      VecCreate(PETSC_COMM_SELF, &vec);
      VecSetFromOptions(vec);
      VecSetType(vec, VECSEQ);
      VecSetSizes(vec, mr, PETSC_DECIDE);
      MatGetDiagonal(mat[b], vec);
      PetscScalar* a;
      VecGetArray(vec, &a);
      // check that we have the correct diagonal values
      std::stringstream str;
      str << "src/LinAlg/Test/refdata/petsc_matrix_diagonal_components_block" << b/2 + 1 << ".ref";
      IntVec v = readIntVector(str.str());
      ASSERT_EQ(mr, (int)v.size());
      for (int i = 0; i < mr; ++i)
        EXPECT_FLOAT_EQ(v[i], a[i]);

      VecRestoreArray(vec, &a);
      VecDestroy(&vec);
    }

    // check that no values outside the diagonal are != 0
    for (PetscInt r = 0; r < mr; ++r) {
      const PetscInt* cols;
      PetscInt ncols;
      const PetscScalar* vals;
      MatGetRow(mat[b], r, &ncols, &cols, &vals);
      for (PetscInt i = 0; i < ncols; ++i)
        if (cols[i] != r)
          EXPECT_FLOAT_EQ(vals[i], 0.0);
      MatRestoreRow(mat[b], r, &ncols, &cols, &vals);
    }
  }
}


TEST(TestPETScMatrix, Copy)
{
  ProcessAdm adm;
  LinSolParams spar;
  PETScMatrix orig(adm, spar);

  orig.init(8);
  for (size_t r = 1; r <= orig.dim(1); ++r)
    for (size_t c = std::max(r - 1, size_t{1}); c <= std::min(r + 1, orig.dim(2)); ++c)
      MatSetValue(orig.getMatrix(), r-c, c-1, r == c ? 1.0 : -2.0, INSERT_VALUES);

  orig.endAssembly();

  std::unique_ptr<PETScMatrix> copy(static_cast<PETScMatrix*>(orig.copy()));
  EXPECT_EQ(orig.dim(1), copy->dim(1));
  EXPECT_EQ(orig.dim(2), copy->dim(2));
  for (size_t r = 1; r <= orig.dim(1); ++r)
    for (size_t c = 1; c <= orig.dim(2); ++c) {
      PetscScalar a, o;
      MatGetValue(orig.getMatrix(), r-1, c-1, &o);
      MatGetValue(copy->getMatrix(), r-1, c-1, &a);
      EXPECT_FLOAT_EQ(o, a);
    }
}


TEST(TestPETScVector, Copy)
{
  ProcessAdm adm;
  adm.dd.setup(5);
  PETScVector orig(adm, 5);

  for (size_t i = 1; i <= orig.dim(); ++i)
    orig(i) = i;

  std::unique_ptr<PETScVector> copy1(static_cast<PETScVector*>(orig.copy()));

  EXPECT_EQ(orig.dim(), copy1->dim());
  PetscScalar* cv;
  VecGetArray(copy1->getVector(), &cv);
  for (size_t i = 1; i <= orig.dim(); ++i) {
      EXPECT_EQ((*copy1)(i), orig(i));
      EXPECT_EQ(cv[i-1], 0.0);
  }

  orig.endAssembly();
  std::unique_ptr<PETScVector> copy2(static_cast<PETScVector*>(orig.copy()));

  EXPECT_EQ(orig.dim(), copy2->dim());
  VecGetArray(copy2->getVector(), &cv);
  for (size_t i = 1; i <= orig.dim(); ++i) {
      EXPECT_EQ((*copy2)(i), orig(i));
      EXPECT_EQ(cv[i-1], orig(i));
  }
}


TEST(TestPETScVectors, Assemble)
{
  SIM2D sim(1);
  sim.read("src/LinAlg/Test/refdata/petsc_test.xinp");
  sim.opt.solver = LinAlg::PETSC;
  ASSERT_TRUE(sim.preprocess());
  ASSERT_TRUE(sim.initSystem(sim.opt.solver));

  PETScMatrix* myMat = dynamic_cast<PETScMatrix*>(sim.getLHSmatrix());
  ASSERT_TRUE(myMat != nullptr);

  myMat->init();
  PETScVectors v(*myMat, 2);

  IntMat mnpc(sim.getPatch(1)->begin_elm(), sim.getPatch(1)->end_elm());
  Vectors el_data(2, Vector(mnpc[0].size()));
  std::iota(el_data[0].begin(), el_data[0].end(), 1.0);
  std::iota(el_data[1].begin(), el_data[1].end(), 5.0);

  v.assemble(el_data, mnpc[0]);
  for (int c = 0; c < 2; ++c) {
    VecAssemblyBegin(v.get(c));
    VecAssemblyEnd(v.get(c));
    for (PetscInt i = 0; i < static_cast<PetscInt>(v.dim()); ++i) {
      PetscScalar a;
      VecGetValues(v.get(c), 1, &i, &a);
      const auto it = std::find(mnpc[0].begin(), mnpc[0].end(), i);
      if (it == mnpc[0].end())
        EXPECT_DOUBLE_EQ(a, 0.0);
      else
        EXPECT_DOUBLE_EQ(a, 1 + 4*c + std::distance(mnpc[0].begin(), it));
    }
  }
}
