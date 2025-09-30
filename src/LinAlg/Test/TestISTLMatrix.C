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

#include "ISTLMatrix.h"
#include "SIM2D.h"
#include "ASMmxBase.h"
#include "SAM.h"
#include "readIntVec.h"
#include <dune/istl/io.hh>

#include "Catch2Support.h"


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


TEST_CASE("TestISTLMatrix.Assemble")
{
  SIM2D sim(1);
  sim.read("src/LinAlg/Test/refdata/petsc_test.xinp");
  sim.opt.solver = LinAlg::ISTL;
  REQUIRE(sim.preprocess());
  REQUIRE(sim.initSystem(sim.opt.solver));

  ISTLMatrix* myMat = dynamic_cast<ISTLMatrix*>(sim.getLHSmatrix());
  REQUIRE(myMat != nullptr);

  Matrix stencil(4,4);
  stencil.diag(1.0);

  for (int iel = 1; iel <= sim.getSAM()->getNoElms(); ++iel)
    myMat->assemble(stencil, *sim.getSAM(), iel);

  myMat->endAssembly();

  // now inspect the matrix
  ISTL::Mat& A = myMat->getMatrix();
  IntVec v = readIntVector("src/LinAlg/Test/refdata/petsc_matrix_diagonal.ref");

  REQUIRE(A.N() == v.size());
  for (size_t i = 0; i < v.size(); ++i)
    REQUIRE_THAT(double(v[i]), WithinRel(A[i][i]));
}


TEST_CASE("TestISTLMatrix.AssembleBasisBlocks")
{
  ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS1;
  ASMmxBase::itgBasis = 2;
  SIM2D sim({1,1});
  sim.read("src/LinAlg/Test/refdata/petsc_test_blocks_basis.xinp");
  sim.opt.solver = LinAlg::ISTL;
  REQUIRE(sim.preprocess());
  REQUIRE(sim.initSystem(sim.opt.solver));

  ISTLMatrix* myMat = dynamic_cast<ISTLMatrix*>(sim.getLHSmatrix());
  REQUIRE(myMat != nullptr);

  Matrix stencil(13,13);
  stencil.diag(1.0);
  for (size_t i = 10; i <= 13; i++)
    stencil(i,i) = 2.0;

  for (int iel = 1; iel <= sim.getSAM()->getNoElms(); ++iel)
    myMat->assemble(stencil, *sim.getSAM(), iel);

  myMat->endAssembly();

  ISTL::Mat& A = myMat->getMatrix();
  const ProcessAdm& adm = sim.getProcessAdm();
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
      REQUIRE(v.size() == mat[b].N());
      for (size_t i = 0; i < v.size(); ++i)
        REQUIRE_THAT(double(v[i]), WithinRel(mat[b][i][i]));
    }

    // check that no values outside the diagonal are != 0
    for (auto r = mat[b].begin(); r != mat[b].end(); ++r)
      for (auto it = r->begin(); it != r->end(); ++it)
        if (r.index() != it.index())
          REQUIRE_THAT(*it, WithinAbs(0.0, 1e-14));
  }
}


TEST_CASE("TestISTLMatrix.AssembleComponentBlocks")
{
  SIM2D sim(2);
  sim.read("src/LinAlg/Test/refdata/petsc_test_blocks_components.xinp");
  sim.opt.solver = LinAlg::ISTL;
  REQUIRE(sim.preprocess());
  REQUIRE(sim.initSystem(sim.opt.solver));

  ISTLMatrix* myMat = dynamic_cast<ISTLMatrix*>(sim.getLHSmatrix());
  REQUIRE(myMat != nullptr);

  Matrix stencil(4*2,4*2);
  for (size_t i = 1; i<= 4; ++i) {
    stencil(2*i-1,2*i-1) = 1.0;
    stencil(2*i,2*i) = 2.0;
  }

  for (int iel = 1; iel <= sim.getSAM()->getNoElms(); ++iel)
    myMat->assemble(stencil, *sim.getSAM(), iel);

  myMat->endAssembly();

  ISTL::Mat& A = myMat->getMatrix();
  const ProcessAdm& adm = sim.getProcessAdm();
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
      REQUIRE(v.size() == mat[b].N());
      for (size_t i = 0; i < v.size(); ++i)
        REQUIRE_THAT(double(v[i]), WithinRel(mat[b][i][i]));
    }

    // check that no values outside the diagonal are != 0
    for (auto r = mat[b].begin(); r != mat[b].end(); ++r)
      for (auto it = r->begin(); it != r->end(); ++it)
        if (r.index() != it.index())
          REQUIRE_THAT(*it, WithinAbs(0.0, 1e-14));
  }
}
