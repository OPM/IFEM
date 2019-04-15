//==============================================================================
//!
//! \file TestISTLMatrix.C
//!
//! \date Apr 27 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for ISTL matrices
//!
//==============================================================================

#include "AlgEqSystem.h"
#include "ISTLMatrix.h"
#include "SIM2D.h"
#include "ASMs2D.h"
#include "IntegrandBase.h"
#include "SAM.h"
#include <dune/istl/io.hh>

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


class InspectBlockPreconditioner : public ISTL::BlockPreconditioner {
public:
  InspectBlockPreconditioner(const ISTL::Mat& A,
                             const DomainDecomposition& dd) :
    ISTL::BlockPreconditioner(A, dd, "upper") {}

  void extractBlock(ISTL::Mat& B, const ISTL::Mat& A,
                    const std::set<int>& eqs_row,
                    const std::set<int>& eqs_col,
                    const std::map<int,int>& eqs_col_g2l)
  {
    ISTL::BlockPreconditioner::extractBlock(B, A, eqs_row, eqs_col, eqs_col_g2l);
  }
};


TEST(TestISTLMatrix, Assemble)
{
  InspectMatrixSIM sim(1);
  sim.read("src/LinAlg/Test/refdata/petsc_test.xinp");
  sim.opt.solver = SystemMatrix::ISTL;
  sim.preprocess();
  sim.initSystem(SystemMatrix::ISTL);

  Matrix stencil(4,4);
  stencil(1,1) = stencil(2,2) = stencil(3,3) = stencil(4,4) = 1.0;

  for (int iel = 1; iel <= sim.getSAM()->getNoElms(); ++iel)
    sim.getMatrix()->assemble(stencil, *sim.getSAM(), iel);

  sim.getMatrix()->beginAssembly();
  sim.getMatrix()->endAssembly();

  // now inspect the matrix
  ISTL::Mat& mat = static_cast<ISTLMatrix*>(sim.getMatrix())->getMatrix();
  IntVec v = readIntVector("src/LinAlg/Test/refdata/petsc_matrix_diagonal.ref");

  ASSERT_EQ(mat.N(), v.size());

  for (size_t i = 0; i < v.size(); ++i)
    ASSERT_FLOAT_EQ(double(v[i]), mat[i][i]);
}


TEST(TestISTLMatrix, AssembleBasisBlocks)
{
  InspectMatrixSIM sim({1,1});
  sim.read("src/LinAlg/Test/refdata/petsc_test_blocks_basis.xinp");
  sim.opt.solver = SystemMatrix::ISTL;
  sim.preprocess();
  sim.initSystem(SystemMatrix::ISTL);

  Matrix stencil(13,13);
  for (size_t i = 1; i<= 9; ++i)
    stencil(i,i) = 1.0;
  for (size_t i = 10; i<= 13; ++i)
    stencil(i,i) = 2.0;

  for (int iel = 1; iel <= sim.getSAM()->getNoElms(); ++iel)
    sim.getMatrix()->assemble(stencil, *sim.getSAM(), iel);

  sim.getMatrix()->beginAssembly();
  sim.getMatrix()->endAssembly();

  const ProcessAdm& adm = sim.getProcessAdm();
  ISTL::Mat& A = static_cast<ISTLMatrix*>(sim.getMatrix())->getMatrix();
  InspectBlockPreconditioner block(A, adm.dd);

  // now inspect the matrix blocks
  std::vector<ISTL::Mat> mat(5);
  block.extractBlock(mat[0], A, adm.dd.getBlockEqs(0), adm.dd.getBlockEqs(0), adm.dd.getG2LEQ(1));
  block.extractBlock(mat[1], A, adm.dd.getBlockEqs(0), adm.dd.getBlockEqs(1), adm.dd.getG2LEQ(2));
  block.extractBlock(mat[2], A, adm.dd.getBlockEqs(1), adm.dd.getBlockEqs(0), adm.dd.getG2LEQ(1));
  block.extractBlock(mat[3], A, adm.dd.getBlockEqs(1), adm.dd.getBlockEqs(1), adm.dd.getG2LEQ(2));
  mat[4] = block.getBlock(1); // S = D

  for (size_t b = 0; b < mat.size(); ++b) {
    if (b == 0 || b >= 3) { // diagonal blocks
      std::stringstream str;
      str << "src/LinAlg/Test/refdata/petsc_matrix_diagonal_basis_block" << std::min(b/2 + 1, 2ul) << ".ref";
      IntVec v = readIntVector(str.str());
      ASSERT_EQ(v.size(), mat[b].N());
      for (size_t i = 0; i < v.size(); ++i)
        ASSERT_FLOAT_EQ(double(v[i]), mat[b][i][i]);
    }

    // check that no values outside the diagonal are != 0
    for (auto r = mat[b].begin(); r != mat[b].end(); ++r) {
      for (auto it = r->begin(); it != r->end(); ++it)
        if (r.index() != it.index())
          ASSERT_FLOAT_EQ(*it, 0.0);
    }
  }
}


TEST(TestISTLMatrix, AssembleComponentBlocks)
{
  InspectMatrixSIM sim(2);
  sim.read("src/LinAlg/Test/refdata/petsc_test_blocks_components.xinp");
  sim.opt.solver = SystemMatrix::ISTL;
  sim.preprocess();
  sim.initSystem(SystemMatrix::ISTL);

  Matrix stencil(4*2,4*2);
  for (size_t i = 1; i<= 4; ++i) {
    stencil(2*i-1,2*i-1) = 1.0;
    stencil(2*i,2*i) = 2.0;
  }

  for (int iel = 1; iel <= sim.getSAM()->getNoElms(); ++iel)
    sim.getMatrix()->assemble(stencil, *sim.getSAM(), iel);

  sim.getMatrix()->beginAssembly();
  sim.getMatrix()->endAssembly();

  const ProcessAdm& adm = sim.getProcessAdm();
  ISTL::Mat& A = static_cast<ISTLMatrix*>(sim.getMatrix())->getMatrix();
  InspectBlockPreconditioner block(A, adm.dd);

  // now inspect the matrix blocks
  std::vector<ISTL::Mat> mat(5);
  block.extractBlock(mat[0], A, adm.dd.getBlockEqs(0), adm.dd.getBlockEqs(0), adm.dd.getG2LEQ(1));
  block.extractBlock(mat[1], A, adm.dd.getBlockEqs(0), adm.dd.getBlockEqs(1), adm.dd.getG2LEQ(2));
  block.extractBlock(mat[2], A, adm.dd.getBlockEqs(1), adm.dd.getBlockEqs(0), adm.dd.getG2LEQ(1));
  block.extractBlock(mat[3], A, adm.dd.getBlockEqs(1), adm.dd.getBlockEqs(1), adm.dd.getG2LEQ(2));
  mat[4] = block.getBlock(1); // S = D

  for (size_t b = 0; b < mat.size(); ++b) {
    if (b == 0 || b >= 3) { // diagonal blocks
      std::stringstream str;
      str << "src/LinAlg/Test/refdata/petsc_matrix_diagonal_components_block" << std::min(b/2 + 1, 2ul) << ".ref";
      IntVec v = readIntVector(str.str());
      ASSERT_EQ(v.size(), mat[b].N());
      for (size_t i = 0; i < v.size(); ++i)
        ASSERT_FLOAT_EQ(double(v[i]), mat[b][i][i]);
    }

    // check that no values outside the diagonal are != 0
    for (auto r = mat[b].begin(); r != mat[b].end(); ++r) {
      for (auto it = r->begin(); it != r->end(); ++it)
        if (r.index() != it.index())
          ASSERT_FLOAT_EQ(*it, 0.0);
    }
  }
}
