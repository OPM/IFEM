//==============================================================================
//!
//! \file TestDomainDecomposition.C
//!
//! \date Mar 29 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for domain decomposition related partitioning for structured models.
//!
//==============================================================================

#include "DomainDecomposition.h"
#include "SAM.h"
#include "ASMbase.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "../../LinAlg/readIntVec.h"

#include "gtest/gtest.h"


template<class Dim>
class TestGlobalLMSIM : public Dim {
public:
  TestGlobalLMSIM(unsigned char n1) : Dim(n1) {}
  TestGlobalLMSIM(const std::vector<unsigned char>& nf) : Dim(nf) {}

protected:
  bool preprocessBeforeAsmInit(int& nnod)
  {
    ++nnod;
    for (int p = 1; p <= this->getNoPatches(); ++p) {
      int idx = this->getLocalPatchIndex(p);
      if (idx > 0)
        this->getPatch(idx)->addGlobalLagrangeMultipliers({nnod}, 1);
    }

    return true;
  }
};


static void check_intvectors_equal (const IntVec& A, const IntVec& B)
{
  ASSERT_EQ(A.size(), B.size());
  for (size_t i = 0; i < A.size(); i++)
    ASSERT_EQ(A[i], B[i]);
};


class TestDomainDecomposition2D : public testing::Test,
                                  public testing::WithParamInterface<int>
{
};


class TestDomainDecomposition3D : public testing::Test,
                                  public testing::WithParamInterface<int>
{
};


TEST_P(TestDomainDecomposition2D, Corner)
{
  std::string file = "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_corner";
  if (GetParam() == 1)
    file += "_fail.xinp";
  else
    file += ".xinp";

  SIM2D sim;
  ASSERT_TRUE(sim.read(file.c_str()));

  if (GetParam() == 0) {
    ASSERT_TRUE(sim.preprocess());
    const ProcessAdm& adm = sim.getProcessAdm();
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_corner_nodes";
    str << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGN(), B);
  } else
    ASSERT_FALSE(sim.preprocess());
}


TEST_P(TestDomainDecomposition2D, CornerLR)
{
  std::stringstream str;
  std::string file = "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_corner";
  if (GetParam() == 1)
    file += "_fail.xinp";
  else
    file += ".xinp";

  SIM2D sim;
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read(file.c_str()));

  if (GetParam() == 0) {
    ASSERT_TRUE(sim.preprocess());
    const ProcessAdm& adm = sim.getProcessAdm();
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_cornerLR_nodes";
    str << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGN(), B);
  } else
    ASSERT_FALSE(sim.preprocess());
}


TEST_P(TestDomainDecomposition2D, SetupSingleBasis)
{
  SIM2D sim(2);
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.preprocess());

  const ProcessAdm& adm = sim.getProcessAdm();
  const SAM* sam = sim.getSAM();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
  str << GetParam() << "_nodes" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGN(), B);
  ASSERT_EQ(sam->getNoNodes(), 9);
  const std::vector<int> maxeqs {16, 26, 36, 44};
  ASSERT_EQ(adm.dd.getMaxEq(), maxeqs[adm.getProcId()]);
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
  str << GetParam() << "_eqs" << adm.getProcId() << ".ref";
  B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGEQ(), B);
  ASSERT_EQ(adm.dd.getMaxDOF(), 2*adm.dd.getMaxNode());
}


TEST_P(TestDomainDecomposition2D, SetupSingleBasisLR)
{
  SIM2D sim(2);
  sim.opt.discretization = ASM::LRSpline;
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.preprocess());

  const ProcessAdm& adm = sim.getProcessAdm();
  const SAM* sam = sim.getSAM();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_LR_orient";
  str << GetParam() << "_nodes" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGN(), B);
  ASSERT_EQ(sam->getNoNodes(), 9);
  const std::vector<int> maxeqs {16, 26, 36, 44};
  ASSERT_EQ(adm.dd.getMaxEq(), maxeqs[adm.getProcId()]);
}


TEST_P(TestDomainDecomposition2D, SetupSingleBasisPeriodic)
{
  SIM2D sim(2);
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_periodic";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.preprocess());

  const ProcessAdm& adm = sim.getProcessAdm();
  const SAM* sam = sim.getSAM();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_periodic";
  str << GetParam() << "_nodes" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGN(), B);
  ASSERT_EQ(sam->getNoNodes(), 9);
}


TEST_P(TestDomainDecomposition2D, SetupSingleBasisGlobalLM)
{
  TestGlobalLMSIM<SIM2D> sim(2);
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.preprocess());

  const ProcessAdm& adm = sim.getProcessAdm();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
  str << GetParam() << "_eqs" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  if (adm.getProcId() == 0)
    B.push_back(17);
  else {
    for (auto& it : B)
      if (it > 16)
        ++it;
    B.push_back(17);
  }
  check_intvectors_equal(adm.dd.getMLGEQ(), B);
}


TEST_P(TestDomainDecomposition2D, SetupSingleBasisBlockEqsComponent)
{
  SIM2D sim(3);
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_blocks_components_orient";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.preprocess());

  const ProcessAdm& adm = sim.getProcessAdm();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
  str << GetParam() << "_eqs" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGEQ(1), B);
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_blocks_components_orient";
  str << GetParam() << "_block2_eqs" << adm.getProcId() << ".ref";
  B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGEQ(2), B);
}


TEST_P(TestDomainDecomposition2D, SetupMixedBasis)
{
  SIM2D sim({2,2});
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.preprocess());

  const ProcessAdm& adm = sim.getProcessAdm();
  const SAM* sam = sim.getSAM();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_mixed_orient";
  str << GetParam() << "_nodes" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGN(), B);
  ASSERT_EQ(sam->getNoNodes(), 25);
  ASSERT_EQ(sam->getNoDOFs(), 50);
  const std::vector<int> maxeqs {48, 81, 114, 140};
  ASSERT_EQ(adm.dd.getMaxEq(), maxeqs[adm.getProcId()]);
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_mixed_orient";
  str << GetParam() << "_eqs" << adm.getProcId() << ".ref";
  B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGEQ(), B);
  ASSERT_EQ(adm.dd.getMinDOF(), 2*(adm.dd.getMinNode()-1)+1);
  ASSERT_EQ(adm.dd.getMaxDOF(), 2*adm.dd.getMaxNode());
}


TEST_P(TestDomainDecomposition2D, SetupMixedBasisLR)
{
  SIM2D sim({2,2});
  sim.opt.discretization = ASM::LRSpline;
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.preprocess());

  const ProcessAdm& adm = sim.getProcessAdm();
  const SAM* sam = sim.getSAM();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_mixed_LR_orient";
  str << GetParam() << "_nodes" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGN(), B);
  ASSERT_EQ(sam->getNoNodes(), 25);
  ASSERT_EQ(sam->getNoDOFs(), 50);
  const std::vector<int> maxeqs {48, 81, 114, 140};
  ASSERT_EQ(adm.dd.getMaxEq(), maxeqs[adm.getProcId()]);
}


TEST_P(TestDomainDecomposition2D, SetupMixedBasisPeriodic)
{
  SIM2D sim({2,2});
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_periodic";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.preprocess());

  const ProcessAdm& adm = sim.getProcessAdm();
  const SAM* sam = sim.getSAM();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_mixed_periodic";
  str << GetParam() << "_nodes" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGN(), B);
  ASSERT_EQ(sam->getNoNodes(), 25);
  ASSERT_EQ(sam->getNoDOFs(), 50);
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_mixed_periodic";
  str << GetParam() << "_eqs" << adm.getProcId() << ".ref";
  B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGEQ(), B);
}


TEST_P(TestDomainDecomposition2D, SetupMixedBasisBlockEqsBasis)
{
  SIM2D sim({2,2});
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_blocks_basis_orient";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.preprocess());

  const ProcessAdm& adm = sim.getProcessAdm();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_mixed_orient";
  str << GetParam() << "_eqs" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGEQ(0), B);
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_blocks_basis_orient";
  str << GetParam() << "_block1_eqs" << adm.getProcId() << ".ref";
  B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGEQ(1), B);
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_blocks_basis_orient";
  str << GetParam() << "_block2_eqs" << adm.getProcId() << ".ref";
  B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGEQ(2), B);
}


TEST_P(TestDomainDecomposition2D, SetupMixedBasisBlockEqsBasisGlobalLM)
{
  TestGlobalLMSIM<SIM2D> sim({2,2});
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_blocks_basis_orient";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.preprocess());

  const ProcessAdm& adm = sim.getProcessAdm();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_mixed_orient";
  str << GetParam() << "_eqs" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  if (adm.getProcId() == 0)
    B.push_back(49);
  else {
    for (auto& it : B)
      if (it > 48)
        ++it;
    B.push_back(49);
  }
  check_intvectors_equal(adm.dd.getMLGEQ(0), B);
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_blocks_basis_orient";
  str << GetParam() << "_block1_eqs" << adm.getProcId() << ".ref";
  B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGEQ(1), B);
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_blocks_basis_orient";
  str << GetParam() << "_block2_eqs" << adm.getProcId() << ".ref";
  B = readIntVector(str.str());
  if (adm.getProcId() == 0)
    B.push_back(19);
  else {
    for (auto& it : B)
      if (it > 18)
        ++it;
    B.push_back(19);
  }
  check_intvectors_equal(adm.dd.getMLGEQ(2), B);
}


TEST_P(TestDomainDecomposition3D, SetupSingleBasis)
{
  SIM3D sim(3);
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.preprocess());

  const ProcessAdm& adm = sim.getProcessAdm();
  const SAM* sam = sim.getSAM();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
  str << GetParam() << "_nodes" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGN(), B);
  ASSERT_EQ(sam->getNoNodes(), 45);
  const std::vector<int> maxeqs {111, 195, 254, 300};
  ASSERT_EQ(adm.dd.getMaxEq(), maxeqs[adm.getProcId()]);
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
  str << GetParam() << "_eqs" << adm.getProcId() << ".ref";
  B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGEQ(), B);
  ASSERT_EQ(adm.dd.getMaxDOF(), 3*adm.dd.getMaxNode());
}


TEST_P(TestDomainDecomposition3D, SetupSingleBasisLR)
{
  SIM3D sim(3);
  sim.opt.discretization = ASM::LRSpline;
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.preprocess());

  const ProcessAdm& adm = sim.getProcessAdm();
  const SAM* sam = sim.getSAM();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_LR_orient";
  str << GetParam() << "_nodes" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGN(), B);
  ASSERT_EQ(sam->getNoNodes(), 45);
  const std::vector<int> maxeqs {111, 195, 254, 300};
  ASSERT_EQ(adm.dd.getMaxEq(), maxeqs[adm.getProcId()]);
}


TEST_P(TestDomainDecomposition3D, SetupSingleBasisPeriodic)
{
  if (GetParam() > 2)
    return;

  SIM3D sim(3);
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_periodic";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.preprocess());

  const ProcessAdm& adm = sim.getProcessAdm();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_periodic";
  str << GetParam() << "_nodes" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGN(), B);
}


TEST_P(TestDomainDecomposition3D, SetupSingleBasisBlockEqsComponent)
{
  if (GetParam() > 0)
    return;

  SIM3D sim(4);
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_blocks_components_orient";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.preprocess());

  const ProcessAdm& adm = sim.getProcessAdm();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
  str << GetParam() << "_eqs" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGEQ(1), B);
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
  str << GetParam() << "_nodes" << adm.getProcId() << ".ref";
  B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGEQ(2), B);
}


TEST_P(TestDomainDecomposition3D, SetupMixedBasis)
{
  if (GetParam() > 0)
    return;

  SIM3D sim({3,1});
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.preprocess());

  const ProcessAdm& adm = sim.getProcessAdm();
  const SAM* sam = sim.getSAM();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_mixed_orient";
  str << GetParam() << "_nodes" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGN(), B);
  ASSERT_EQ(sam->getNoNodes(), 157);
  const std::vector<int> maxeqs {337, 607, 828, 1007};
  ASSERT_EQ(adm.dd.getMaxEq(), maxeqs[adm.getProcId()]);
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_mixed_orient";
  str << GetParam() << "_eqs" << adm.getProcId() << ".ref";
  B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGEQ(), B);
}


TEST_P(TestDomainDecomposition3D, SetupMixedBasisPeriodic)
{
  if (GetParam() > 2)
    return;

  SIM3D sim({1,1});
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_periodic";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.preprocess());

  const ProcessAdm& adm = sim.getProcessAdm();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_mixed_periodic";
  str << GetParam() << "_nodes" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());

  check_intvectors_equal(adm.dd.getMLGN(), B);
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_mixed_periodic";
  str << GetParam() << "_eqs" << adm.getProcId() << ".ref";
  B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGEQ(), B);
}


TEST_P(TestDomainDecomposition3D, SetupMixedBasisPeriodicLM)
{
  if (GetParam() > 0)
    return;

  TestGlobalLMSIM<SIM3D> sim({1,1});
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_periodiclm";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.preprocess());

  const ProcessAdm& adm = sim.getProcessAdm();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_mixed_periodiclm";
  str << GetParam() << "_nodes" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGN(), B);

  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_mixed_periodiclm";
  str << GetParam() << "_eqs" << adm.getProcId() << ".ref";
  B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGEQ(), B);
}


TEST_P(TestDomainDecomposition3D, Corner)
{
  if (GetParam() > 1)
    return;

  SIM3D sim;
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_corner";
  if (GetParam() > 0)
    str << "_fail";
  str << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));

  if (GetParam() < 1) {
    ASSERT_TRUE(sim.preprocess());
    const ProcessAdm& adm = sim.getProcessAdm();
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_corner_nodes";
    str << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGN(), B);
  } else
    ASSERT_FALSE(sim.preprocess());
}


TEST_P(TestDomainDecomposition3D, CornerLR)
{
  if (GetParam() > 1)
    return;

  SIM3D sim;
  sim.opt.discretization = ASM::LRSpline;

  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_corner";
  if (GetParam() > 0)
    str << "_fail";
  str << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));

  if (GetParam() < 1) {
    ASSERT_TRUE(sim.preprocess());
    const ProcessAdm& adm = sim.getProcessAdm();
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_cornerLR_nodes";
    str << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGN(), B);
  } else
    ASSERT_FALSE(sim.preprocess());
}

const std::vector<int> orientations2D = {0,1};
INSTANTIATE_TEST_SUITE_P(TestDomainDecomposition2D, TestDomainDecomposition2D, testing::ValuesIn(orientations2D));
const std::vector<int> orientations3D = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
INSTANTIATE_TEST_SUITE_P(TestDomainDecomposition3D, TestDomainDecomposition3D, testing::ValuesIn(orientations3D));
