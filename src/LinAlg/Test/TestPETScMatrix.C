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

#include "AlgEqSystem.h"
#include "DenseMatrix.h"
#include "PETScMatrix.h"
#include "SIM2D.h"
#include "ASMs2D.h"
#include "IntegrandBase.h"
#include "SAM.h"

#include "gtest/gtest.h"

#include <fstream>


typedef std::vector<int> IntVec;

static IntVec readIntVector(const std::string& file)
{
  std::vector<int> result;
  std::ifstream f(file);
  size_t size;
  f >> size;
  result.resize(size);
  for (size_t j=0;j<size;++j)
    f >> result[j];

  return result;
}


class DummyIntegrandForDummies : public IntegrandBase {
public:
  DummyIntegrandForDummies(unsigned short int n = 0) : IntegrandBase(n) {}
};


class InspectMatrixSIM : public SIM2D {
public:
  InspectMatrixSIM(unsigned char n1 = 2, bool check = false) :
    SIM2D(n1, check) { myProblem = new DummyIntegrandForDummies; }

  InspectMatrixSIM(const std::vector<unsigned char>& n, bool check = false) :
    SIM2D(n, check) { myProblem = new DummyIntegrandForDummies; }

  SystemMatrix* getMatrix() { return myEqSys->getMatrix(0); }
};


TEST(TestPETScMatrix, Assemble)
{
  InspectMatrixSIM sim(1);
  sim.read("src/LinAlg/Test/refdata/petsc_test.xinp");
  sim.opt.solver = SystemMatrix::PETSC;
  sim.preprocess();
  sim.initSystem(SystemMatrix::PETSC);

  Matrix stencil(4,4);
  stencil(1,1) = stencil(2,2) = stencil(3,3) = stencil(4,4) = 1.0;

  for (int iel = 1; iel <= sim.getSAM()->getNoElms(); ++iel)
    sim.getMatrix()->assemble(stencil, *sim.getSAM(), iel);

  sim.getMatrix()->beginAssembly();
  sim.getMatrix()->endAssembly();

  // now inspect the matrix
  Mat& mat = static_cast<PETScMatrix*>(sim.getMatrix())->getMatrix();

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
  std::vector<int> ref;
  // check that we have the correct diagonal values
  IntVec v = readIntVector("src/LinAlg/Test/refdata/petsc_matrix_diagonal.ref");
  ASSERT_EQ(mr, (int)v.size());
  for (int i = 0; i < mr; ++i)
    ASSERT_FLOAT_EQ(double(v[i]), a[i]);

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
        ASSERT_FLOAT_EQ(vals[i], 0.0);
    MatRestoreRow(mat, r, &ncols, &cols, &vals);
  }
}


TEST(TestPETScMatrix, AssembleBasisBlocks)
{
  InspectMatrixSIM sim({1,1});
  sim.read("src/LinAlg/Test/refdata/petsc_test_blocks_basis.xinp");
  sim.opt.solver = SystemMatrix::PETSC;
  sim.preprocess();
  sim.initSystem(SystemMatrix::PETSC);

  Matrix stencil(13,13);
  for (size_t i = 1; i<= 9; ++i)
    stencil(i,i) = 1.0;
  for (size_t i = 10; i<= 13; ++i)
    stencil(i,i) = 2.0;

  for (int iel = 1; iel <= sim.getSAM()->getNoElms(); ++iel)
    sim.getMatrix()->assemble(stencil, *sim.getSAM(), iel);

  sim.getMatrix()->beginAssembly();
  sim.getMatrix()->endAssembly();

  // now inspect the matrix blocks
  const std::vector<Mat>& mat = static_cast<PETScMatrix*>(sim.getMatrix())->getBlockMatrices();

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
        ASSERT_FLOAT_EQ(double(v[i]), a[i]);

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
          ASSERT_FLOAT_EQ(vals[i], 0.0);
      MatRestoreRow(mat[b], r, &ncols, &cols, &vals);
    }
  }
}


TEST(TestPETScMatrix, AssembleComponentBlocks)
{
  InspectMatrixSIM sim(2);
  sim.read("src/LinAlg/Test/refdata/petsc_test_blocks_components.xinp");
  sim.opt.solver = SystemMatrix::PETSC;
  sim.preprocess();
  sim.initSystem(SystemMatrix::PETSC);

  Matrix stencil(4*2,4*2);
  for (size_t i = 1; i<= 4; ++i) {
    stencil(2*i-1,2*i-1) = 1.0;
    stencil(2*i,2*i) = 2.0;
  }

  for (int iel = 1; iel <= sim.getSAM()->getNoElms(); ++iel)
    sim.getMatrix()->assemble(stencil, *sim.getSAM(), iel);

  sim.getMatrix()->beginAssembly();
  sim.getMatrix()->endAssembly();

  // now inspect the matrix blocks
  const std::vector<Mat>& mat = static_cast<PETScMatrix*>(sim.getMatrix())->getBlockMatrices();

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
      std::vector<int> ref;
      // check that we have the correct diagonal values
      std::stringstream str;
      str << "src/LinAlg/Test/refdata/petsc_matrix_diagonal_components_block" << b/2 + 1 << ".ref";
      IntVec v = readIntVector(str.str());
      ASSERT_EQ(mr, (int)v.size());
      for (int i = 0; i < mr; ++i)
        ASSERT_FLOAT_EQ(double(v[i]), a[i]);

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
          ASSERT_FLOAT_EQ(vals[i], 0.0);
      MatRestoreRow(mat[b], r, &ncols, &cols, &vals);
    }
  }
}
