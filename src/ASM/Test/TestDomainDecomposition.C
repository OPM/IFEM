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
#include "gtest/gtest.h"
#include <fstream>


typedef std::vector< std::vector<int> > IntMat;

static IntMat readIntMatrix(size_t r, const std::string& file)
{
  std::vector< std::vector<int> > result;
  result.resize(r);
  std::ifstream f(file);
  for (size_t i=0;i<r;++i) {
    size_t size;
    f >> size;
    result[i].resize(size);
    for (size_t j=0;j<size;++j)
      f >> result[i][j];
  }

  return result;
}


auto&& check_intmatrices_equal = [](const std::vector<std::vector<int>>& subdomains,
                                    const std::string& path)
{
  IntMat B = readIntMatrix(subdomains.size(), path);
  for (size_t i = 0; i < subdomains.size(); ++i) {
    size_t j = 0;
    for (const auto& it2 : subdomains[i])
      ASSERT_EQ(it2, B[i][j++]);
  }
};


TEST(TestDomainDecomposition, LocalGroups1DO1)
{
  auto domains = DomainDecomposition::calcSubdomains(7, 0, 0, 3, 0, 0, 1);
  check_intmatrices_equal(domains, "src/ASM/Test/refdata/DomainDecomposition_1D_O1.ref");
}


TEST(TestDomainDecomposition, LocalGroups1DO2)
{
  auto domains = DomainDecomposition::calcSubdomains(7, 0, 0, 3, 0, 0, 2);
  check_intmatrices_equal(domains, "src/ASM/Test/refdata/DomainDecomposition_1D_O2.ref");
}


TEST(TestDomainDecomposition, LocalGroups2DO1)
{
  auto domains = DomainDecomposition::calcSubdomains(7, 7, 0, 3, 3, 0, 1);
  check_intmatrices_equal(domains, "src/ASM/Test/refdata/DomainDecomposition_2D_O1.ref");
}


TEST(TestDomainDecomposition, LocalGroups2DO2)
{
  auto domains = DomainDecomposition::calcSubdomains(7, 7, 0, 3, 3, 0, 2);
  check_intmatrices_equal(domains, "src/ASM/Test/refdata/DomainDecomposition_2D_O2.ref");
}


TEST(TestDomainDecomposition, LocalGroups3DO1)
{
  auto domains = DomainDecomposition::calcSubdomains(7, 7, 7, 3, 3, 3, 1);
  check_intmatrices_equal(domains, "src/ASM/Test/refdata/DomainDecomposition_3D_O1.ref");
}


TEST(TestDomainDecomposition, LocalGroups3DO2)
{
  auto domains = DomainDecomposition::calcSubdomains(7, 7, 7, 3, 3, 3, 2);
  check_intmatrices_equal(domains, "src/ASM/Test/refdata/DomainDecomposition_3D_O2.ref");
}


TEST(TestDomainDecomposition, Setup)
{
  SIM2D sim(1);
  sim.read("src/ASM/Test/refdata/DomainDecomposition_2D_1P.xinp");
  sim.preprocess();

  DomainDecomposition dd;
  // TODO: Remove after integration
  dd.setup(sim.getProcessAdm(), sim);
  const SAM* sam = sim.getSAM();

  ASSERT_EQ(dd.getMinEq(), 1);
  ASSERT_EQ(dd.getMaxEq(), sam->getNoEquations());
  ASSERT_EQ(dd.getMinNode(), 1);
  ASSERT_EQ(dd.getMaxNode(), sam->getNoNodes());
  ASSERT_EQ(dd.getMinDOF(), 1);
  ASSERT_EQ(dd.getMaxDOF(), sam->getNoDOFs());
  for (int i = 1; i <= sam->getNoEquations(); ++i)
    ASSERT_EQ(dd.getGlobalEq(i), i);
  ASSERT_EQ(dd.getGlobalEq(sam->getNoEquations()+1), 0);

  // disabled until integration in the SIM class is added
//  ASSERT_EQ(dd.getPatchOwner(1), 0);
//  ASSERT_EQ(dd.getPatchOwner(2), -1);

  ASSERT_TRUE(dd.getMLGEQ().empty());
  ASSERT_TRUE(dd.getMLGN().empty());
}
