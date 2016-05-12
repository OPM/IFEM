//==============================================================================
//!
//! \file TestISTLPETScMatrix.C
//!
//! \date Apr 27 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests comparing ISTL and PETSc matrices
//!
//==============================================================================

#include "AlgEqSystem.h"
#include "ISTLMatrix.h"
#include "PETScMatrix.h"
#include "SIM2D.h"
#include "ASMs2D.h"
#include "IntegrandBase.h"
#include "SAM.h"
#include <dune/istl/io.hh>

#include "gtest/gtest.h"

#include <fstream>


class DummyIntegrandForDummies : public IntegrandBase {
public:
  DummyIntegrandForDummies(unsigned short int n = 0) : IntegrandBase(n) {}
};


class InspectMatrixSIM : public SIM2D {
public:
  InspectMatrixSIM() :
    SIM2D({1,1}, false) { myProblem = new DummyIntegrandForDummies; }

  SystemMatrix* getMatrix() { return myEqSys->getMatrix(0); }
};



TEST(TestISTLPETScMatrix, SchurComplement)
{
  Matrix stencil(13,13);
  for (size_t i = 1; i<= 13; ++i)
    for (size_t j = 1; j <= 13; ++j)
      stencil(i,j) = 1.0;

  std::array<InspectMatrixSIM,2> sim;
  for (size_t i = 0; i < 2; ++i) {
    sim[i].read("src/LinAlg/Test/refdata/petsc_test_blocks_basis.xinp");
    sim[i].opt.solver = i == 0 ? SystemMatrix::PETSC : SystemMatrix::ISTL;
    sim[i].preprocess();
    sim[i].initSystem(i == 0 ? SystemMatrix::PETSC : SystemMatrix::ISTL);

    for (int iel = 1; iel <= sim[i].getSAM()->getNoElms(); ++iel)
      sim[i].getMatrix()->assemble(stencil, *sim[i].getSAM(), iel);

    sim[i].getMatrix()->beginAssembly();
    sim[i].getMatrix()->endAssembly();
  }

  const ProcessAdm& adm = sim[1].getProcessAdm();
  ISTL::Mat& A = static_cast<ISTLMatrix*>(sim[1].getMatrix())->getMatrix();
  ISTL::BlockPreconditioner block(A, adm.dd, "upper");

  ISTL::Mat& S = block.getBlock(1);
  PETScSolParams params(LinSolParams(), adm);
  params.setupSchurComplement(static_cast<PETScMatrix*>(sim[0].getMatrix())->getBlockMatrices());

  // check that matrices are the same
  for (size_t r = 0; r < S.N(); ++r) {
    const PetscInt* cols;
    PetscInt ncols;
    const PetscScalar* vals;
    MatGetRow(params.getSchurComplement(), r, &ncols, &cols, &vals);
    for (PetscInt i = 0; i < ncols; ++i)
      ASSERT_FLOAT_EQ(vals[i], S[r][cols[i]]);
    MatRestoreRow(params.getSchurComplement(), r, &ncols, &cols, &vals);
  }
}
