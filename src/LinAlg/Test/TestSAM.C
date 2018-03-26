//==============================================================================
//!
//! \file TestSAM.C
//!
//! \date Mar 29 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for SAM
//!
//==============================================================================

#include "SAM.h"
#include "SIM2D.h"
#include "ASMbase.h"
#include "ASMmxBase.h"

#include "gtest/gtest.h"

#include <fstream>

typedef std::vector<IntVec> IntMat;


static std::istream& operator>> (std::istream& is, IntVec& v)
{
  size_t nv;
  is >> nv;
  v.resize(nv);
  for (size_t i = 0; i < nv; i++)
    is >> v[i];

  return is;
}


static IntMat readIntMatrix (size_t r, const std::string& file)
{
  IntMat result(r);
  std::ifstream f(file);
  for (size_t i = 0; i < r; i++)
    f >> result[i];

  return result;
}


static void check_intmatrices_equal (const SAM* sam, const std::string& path)
{
  std::vector<IntSet> A;
  ASSERT_TRUE(sam->getDofCouplings(A));

  size_t i, j;
  IntSet::const_iterator it;
  IntMat B = readIntMatrix(A.size(), path);
  for (i = 0; i < A.size(); i++)
    for (it = A[i].begin(), j = 0; it != A[i].end(); ++it, j++)
      ASSERT_EQ(*it, B[i][j]);
};


TEST(TestSAM, SingleBasis1P)
{
  SIM2D sim(1);
  sim.read("src/LinAlg/Test/refdata/sam_2D_1P.xinp");
  sim.preprocess();

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_singlebasis_1P.ref");

  for (int i = 1; i <= sam->getNoNodes(); ++i)
    ASSERT_EQ(sam->getEquation(i, 1), i);
}


TEST(TestSAM, SingleBasisDirichlet1P)
{
  SIM2D sim(1);
  sim.read("src/LinAlg/Test/refdata/sam_2D_dir_1P.xinp");
  sim.preprocess();

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_singlebasis_dir_1P.ref");

  int i, eq = 1;
  for (i = 1; i <= 3; ++i) {
    ASSERT_EQ(sam->getEquation(3*(i-1)+1, 1), eq++);
    ASSERT_EQ(sam->getEquation(3*(i-1)+2, 1), eq++);
    ASSERT_EQ(sam->getEquation(3*(i-1)+3, 1), 0);
  }
}


TEST(TestSAM, MixedBasis1P)
{
  ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS1;
  ASMmxBase::elmBasis = 2;

  SIM2D sim({1,1});
  sim.read("src/LinAlg/Test/refdata/sam_2D_1P.xinp");
  sim.preprocess();

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_mixedbasis_1P.ref");

  for (int i = 1; i <= sam->getNoNodes(); ++i)
    ASSERT_EQ(sam->getEquation(i, 1), i);
}


TEST(TestSAM, MixedBasisDirichlet1P)
{
  SIM2D sim({1,1});
  sim.read("src/LinAlg/Test/refdata/sam_2D_dir_1P.xinp");
  sim.preprocess();

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_mixedbasis_dir_1P.ref");

  int i, eq = 1;
  for (i = 1; i <= 4; ++i) {
    ASSERT_EQ(sam->getEquation(4*(i-1)+1, 1), eq++);
    ASSERT_EQ(sam->getEquation(4*(i-1)+2, 1), eq++);
    ASSERT_EQ(sam->getEquation(4*(i-1)+3, 1), eq++);
    ASSERT_EQ(sam->getEquation(4*(i-1)+4, 1), 0);
  }

  for (i = 17; i <= 25; ++i)
    ASSERT_EQ(sam->getEquation(i, 1), eq++);
}


TEST(TestSAM, SingleBasis2P)
{
  SIM2D sim(1);
  sim.read("src/LinAlg/Test/refdata/sam_2D_2P.xinp");
  sim.preprocess();

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_singlebasis_2P.ref");

  for (int i = 1; i <= sam->getNoNodes(); ++i)
    ASSERT_EQ(sam->getEquation(i, 1), i);
}


TEST(TestSAM, MixedBasis2P)
{
  SIM2D sim({1,1});
  sim.read("src/LinAlg/Test/refdata/sam_2D_2P.xinp");
  sim.preprocess();

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_mixedbasis_2P.ref");

  for (int i = 1; i <= sam->getNoNodes(); ++i)
    ASSERT_EQ(sam->getEquation(i, 1), i);
}


TEST(TestSAM, SingleBasisDirichlet2P)
{
  SIM2D sim(1);
  sim.read("src/LinAlg/Test/refdata/sam_2D_dir_2P.xinp");
  sim.preprocess();

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_singlebasis_dir_2P.ref");

  ASSERT_EQ(sam->getEquation(1, 1), 1);
  ASSERT_EQ(sam->getEquation(2, 1), 2);
  ASSERT_EQ(sam->getEquation(3, 1), 3);
  ASSERT_EQ(sam->getEquation(4, 1), 4);
  ASSERT_EQ(sam->getEquation(5, 1), 0);
  ASSERT_EQ(sam->getEquation(6, 1), 0);
}


TEST(TestSAM, MixedBasisDirichlet2P)
{
  SIM2D sim({1,1});
  sim.read("src/LinAlg/Test/refdata/sam_2D_dir_2P.xinp");
  sim.preprocess();

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_mixedbasis_dir_2P.ref");

  int i, eq = 1;
  for (i = 1; i <= 13; ++i)
    ASSERT_EQ(sam->getEquation(i, 1), eq++);

  const ASMbase* pch = sim.getPatch(2);
  for (i = 1; i <= 3; ++i) {
    ASSERT_EQ(sam->getEquation(pch->getNodeID((i-1)*3+2), 1), eq++);
    ASSERT_EQ(sam->getEquation(pch->getNodeID((i-1)*3+3), 1), 0);
  }
  ASSERT_EQ(sam->getEquation(20, 1), eq++);
  ASSERT_EQ(sam->getEquation(21, 1), eq++);
}
