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

#include "AlgEqSystem.h"
#include "PETScMatrix.h"
#include "SIM2D.h"
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


TEST(TestPETScMatrix, AssembleMPI)
{
  InspectMatrixSIM sim(1);
  sim.read("src/LinAlg/Test/refdata/petsc_test.xinp");
  sim.opt.solver = LinAlg::PETSC;
  ASSERT_TRUE(sim.preprocess());
  ASSERT_TRUE(sim.initSystem(sim.opt.solver));

  Matrix stencil(4,4);
  stencil(1,1) = stencil(2,2) = stencil(3,3) = stencil(4,4) = 1.0;

  for (int iel = 1; iel <= sim.getSAM()->getNoElms(); ++iel)
    sim.getMatrix()->assemble(stencil, *sim.getSAM(), iel);

  sim.getMatrix()->beginAssembly();
  sim.getMatrix()->endAssembly();

  // now inspect the matrix
  Mat& mat = static_cast<PETScMatrix*>(sim.getMatrix())->getMatrix();

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
    ASSERT_FLOAT_EQ(v[i], a[i-r]);

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
        ASSERT_FLOAT_EQ(vals[i], 0.0);
    MatRestoreRow(mat, r, &ncols, &cols, &vals);
  }
}


TEST(TestPETScMatrix, AssembleBasisBlocksMPI)
{
  InspectMatrixSIM sim({1,1});
  sim.read("src/LinAlg/Test/refdata/petsc_test_blocks_basis.xinp");
  sim.opt.solver = LinAlg::PETSC;
  ASSERT_TRUE(sim.preprocess());
  ASSERT_TRUE(sim.initSystem(sim.opt.solver));

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
        ASSERT_FLOAT_EQ(double(v[i]), a[i-r]);

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
          ASSERT_FLOAT_EQ(vals[i], 0.0);
      MatRestoreRow(mat[b], r+adm.dd.getMinEq(b/2+1)-1, &ncols, &cols, &vals);
    }
  }
}


TEST(TestPETScMatrix, AssembleComponentBlocksMPI)
{
  InspectMatrixSIM sim(2);
  sim.read("src/LinAlg/Test/refdata/petsc_test_blocks_components.xinp");
  sim.opt.solver = LinAlg::PETSC;
  ASSERT_TRUE(sim.preprocess());
  ASSERT_TRUE(sim.initSystem(sim.opt.solver));

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
        ASSERT_FLOAT_EQ(double(v[i]), a[i-r]);

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
          ASSERT_FLOAT_EQ(vals[i], 0.0);
      MatRestoreRow(mat[b], r+adm.dd.getMinEq(b/2+1)-1, &ncols, &cols, &vals);
    }
  }
}
