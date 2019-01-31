//==============================================================================
//!
//! \file TestISTLMatrix.C
//!
//! \date Apr 27 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for parallel ISTL matrices
//!
//==============================================================================

#include "ISTLMatrix.h"
#include "SIM2D.h"
#include "SAM.h"
#include "readIntVec.h"

#include "gtest/gtest.h"


TEST(TestISTLMatrix, AssembleMPI)
{
  SIM2D sim(1);
  sim.read("src/LinAlg/Test/refdata/petsc_test.xinp");
  sim.opt.solver = LinAlg::ISTL;
  ASSERT_TRUE(sim.preprocess());
  ASSERT_TRUE(sim.initSystem(sim.opt.solver));

  ISTLMatrix* myMat = dynamic_cast<ISTLMatrix*>(sim.getLHSmatrix());
  ASSERT_TRUE(myMat != nullptr);

  Matrix stencil(4,4);
  stencil.diag(1.0);

  for (int iel = 1; iel <= sim.getSAM()->getNoElms(); ++iel)
    myMat->assemble(stencil, *sim.getSAM(), iel);

  myMat->beginAssembly();
  myMat->endAssembly();

  // now inspect the matrix
  const ProcessAdm& adm = sim.getProcessAdm();
  ISTL::Mat& mat = myMat->getMatrix();
  ISTL::Vec b(mat.N()), b2(mat.N());

  try {
    Dune::OwnerOverlapCopyCommunication<int,int> comm(*adm.getCommunicator());
    comm.indexSet().beginResize();
    typedef Dune::ParallelLocalIndex<Dune::OwnerOverlapCopyAttributeSet::AttributeSet> LI;
    for (size_t i = 0; i < adm.dd.getMLGEQ().size(); ++i) {
      int gid = adm.dd.getGlobalEq(i+1);
      comm.indexSet().add(gid-1, LI(i, gid >= adm.dd.getMinEq() ?
                                       Dune::OwnerOverlapCopyAttributeSet::owner :
                                       Dune::OwnerOverlapCopyAttributeSet::overlap));
    }
    comm.indexSet().endResize();
    comm.remoteIndices().setIncludeSelf(true);
    comm.remoteIndices().template rebuild<false>();

    ISTL::ParMatrixAdapter op(mat, comm);

    b = 1.0;
    op.apply(b, b2);
  } catch (Dune::ISTLError& e) {
    std::cerr << e << std::endl;
    ASSERT_TRUE(false);
  }

  IntVec v = readIntVector("src/LinAlg/Test/refdata/petsc_matrix_diagonal.ref");
  for (size_t i = 1; i <= adm.dd.getMLGEQ().size(); ++i)
    EXPECT_FLOAT_EQ(v[adm.dd.getGlobalEq(i)-1], b2[i-1]);
}
