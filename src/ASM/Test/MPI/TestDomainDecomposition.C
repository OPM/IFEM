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
#include "SIM2D.h"
#include "SIM3D.h"
#include "IFEM.h"
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


auto&& check_intvectors_equal = [](const IntVec& A,
                                   const IntVec& B)
{
  ASSERT_EQ(A.size(), B.size());
  auto it = B.begin();
  for (auto& it2 : A) {
    ASSERT_EQ(it2, *it);
    ++it;
   }
};


class TestDomainDecomposition2D : public testing::Test,
                                  public testing::WithParamInterface<int>
{
};


class TestDomainDecomposition3D : public testing::Test,
                                  public testing::WithParamInterface<int>
{
};


TEST_P(TestDomainDecomposition2D, SetupSingleBasis)
{
  SIM2D sim(2);
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
  str << GetParam() << ".xinp";
  sim.read(str.str().c_str());
  sim.preprocess();

  const ProcessAdm& adm = sim.getProcessAdm();
  // TODO: unnecessary cast and setup call after integration.
  DomainDecomposition& dd = const_cast<DomainDecomposition&>(adm.dd);
  dd.setup(adm, sim);
  const SAM* sam = sim.getSAM();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
  str << GetParam() << "_nodes" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  check_intvectors_equal(dd.getMLGN(), B);
  ASSERT_EQ(sam->getNoNodes(), 9);
  const std::vector<int> maxeqs {16, 26, 36, 44};
  ASSERT_EQ(dd.getMaxEq(), maxeqs[adm.getProcId()]);
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
  str << GetParam() << "_eqs" << adm.getProcId() << ".ref";
  B = readIntVector(str.str());
  check_intvectors_equal(dd.getMLGEQ(), B);
  ASSERT_EQ(dd.getMaxDOF(), 2*dd.getMaxNode());
}


TEST_P(TestDomainDecomposition2D, SetupMixedBasis)
{
  SIM2D sim({2,2});
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
  str << GetParam() << ".xinp";
  sim.read(str.str().c_str());
  sim.preprocess();

  const ProcessAdm& adm = sim.getProcessAdm();
  // TODO: unnecessary cast and setup call after integration.
  DomainDecomposition& dd = const_cast<DomainDecomposition&>(adm.dd);
  dd.setup(adm, sim);
  const SAM* sam = sim.getSAM();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_mixed_orient";
  str << GetParam() << "_nodes" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  check_intvectors_equal(dd.getMLGN(), B);
  ASSERT_EQ(sam->getNoNodes(), 25);
  ASSERT_EQ(sam->getNoDOFs(), 50);
  const std::vector<int> maxeqs {48, 81, 114, 140};
  ASSERT_EQ(dd.getMaxEq(), maxeqs[adm.getProcId()]);
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_mixed_orient";
  str << GetParam() << "_eqs" << adm.getProcId() << ".ref";
  B = readIntVector(str.str());
  check_intvectors_equal(dd.getMLGEQ(), B);
  ASSERT_EQ(dd.getMinDOF(), 2*(dd.getMinNode()-1)+1);
  ASSERT_EQ(dd.getMaxDOF(), 2*dd.getMaxNode());
}


TEST_P(TestDomainDecomposition3D, SetupSingleBasis)
{
  SIM3D sim(3);
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
  str << GetParam() << ".xinp";
  sim.read(str.str().c_str());
  sim.preprocess();

  const ProcessAdm& adm = sim.getProcessAdm();
  // TODO: unnecessary cast and setup call after integration.
  DomainDecomposition& dd = const_cast<DomainDecomposition&>(adm.dd);
  dd.setup(adm, sim);
  const SAM* sam = sim.getSAM();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
  str << GetParam() << "_nodes" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  check_intvectors_equal(dd.getMLGN(), B);
  ASSERT_EQ(sam->getNoNodes(), 45);
  const std::vector<int> maxeqs {111, 195, 254, 300};
  ASSERT_EQ(dd.getMaxEq(), maxeqs[adm.getProcId()]);
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
  str << GetParam() << "_eqs" << adm.getProcId() << ".ref";
  B = readIntVector(str.str());
  check_intvectors_equal(dd.getMLGEQ(), B);
  ASSERT_EQ(dd.getMaxDOF(), 3*dd.getMaxNode());
}


TEST_P(TestDomainDecomposition3D, SetupMixedBasis)
{
  if (GetParam() > 0)
    return;

  SIM3D sim({3,1});
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
  str << GetParam() << ".xinp";
  sim.read(str.str().c_str());
  sim.preprocess();

  const ProcessAdm& adm = sim.getProcessAdm();
  // TODO: unnecessary cast and setup call after integration.
  DomainDecomposition& dd = const_cast<DomainDecomposition&>(adm.dd);
  dd.setup(adm, sim);
  const SAM* sam = sim.getSAM();
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_mixed_orient";
  str << GetParam() << "_nodes" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  check_intvectors_equal(dd.getMLGN(), B);
  ASSERT_EQ(sam->getNoNodes(), 157);
  const std::vector<int> maxeqs {337, 607, 828, 1007};
  ASSERT_EQ(dd.getMaxEq(), maxeqs[adm.getProcId()]);
  str.str("");
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_mixed_orient";
  str << GetParam() << "_eqs" << adm.getProcId() << ".ref";
  B = readIntVector(str.str());
  check_intvectors_equal(dd.getMLGEQ(), B);
}


const std::vector<int> orientations2D = {0,1};
INSTANTIATE_TEST_CASE_P(TestDomainDecomposition, TestDomainDecomposition2D, testing::ValuesIn(orientations2D));
const std::vector<int> orientations3D = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
INSTANTIATE_TEST_CASE_P(TestDomainDecomposition, TestDomainDecomposition3D, testing::ValuesIn(orientations3D));
