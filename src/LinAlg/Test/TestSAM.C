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

#include <catch2/catch_test_macros.hpp>

#include <fstream>

using IntMat = std::vector<IntVec>;


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
  sam->getDofCouplings(A);

  size_t i, j;
  IntSet::const_iterator it;
  IntMat B = readIntMatrix(A.size(), path);
  for (i = 0; i < A.size(); i++)
    for (it = A[i].begin(), j = 0; it != A[i].end(); ++it, j++)
      REQUIRE(*it == B[i][j]);
};


TEST_CASE("TestSAM.SingleBasis1P")
{
  SIM2D sim(1);
  sim.read("src/LinAlg/Test/refdata/sam_2D_1P.xinp");
  sim.preprocess();

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_singlebasis_1P.ref");

  for (int i = 1; i <= sam->getNoNodes(); ++i)
    REQUIRE(sam->getEquation(i, 1) == i);
}


TEST_CASE("TestSAM.SingleBasisDirichlet1P")
{
  SIM2D sim(1);
  sim.read("src/LinAlg/Test/refdata/sam_2D_dir_1P.xinp");
  sim.preprocess();

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_singlebasis_dir_1P.ref");

  int i, eq = 1;
  for (i = 1; i <= 3; ++i) {
    REQUIRE(sam->getEquation(3*(i-1)+1, 1) == eq++);
    REQUIRE(sam->getEquation(3*(i-1)+2, 1) == eq++);
    REQUIRE(sam->getEquation(3*(i-1)+3, 1) == 0);
  }
}


TEST_CASE("TestSAM.MixedBasis1P")
{
  ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS1;
  ASMmxBase::itgBasis = 2;

  SIM2D sim({1,1});
  sim.read("src/LinAlg/Test/refdata/sam_2D_1P.xinp");
  sim.preprocess();

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_mixedbasis_1P.ref");

  for (int i = 1; i <= sam->getNoNodes(); ++i)
    REQUIRE(sam->getEquation(i, 1) == i);
}


TEST_CASE("TestSAM.MixedBasisDirichlet1P")
{
  SIM2D sim({1,1});
  sim.read("src/LinAlg/Test/refdata/sam_2D_dir_1P.xinp");
  sim.preprocess();

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_mixedbasis_dir_1P.ref");

  int i, eq = 1;
  for (i = 1; i <= 4; ++i) {
    REQUIRE(sam->getEquation(4*(i-1)+1, 1) == eq++);
    REQUIRE(sam->getEquation(4*(i-1)+2, 1) == eq++);
    REQUIRE(sam->getEquation(4*(i-1)+3, 1) == eq++);
    REQUIRE(sam->getEquation(4*(i-1)+4, 1) == 0);
  }

  for (i = 17; i <= 25; ++i)
    REQUIRE(sam->getEquation(i, 1) == eq++);
}


TEST_CASE("TestSAM.SingleBasis2P")
{
  SIM2D sim(1);
  sim.read("src/LinAlg/Test/refdata/sam_2D_2P.xinp");
  sim.preprocess();

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_singlebasis_2P.ref");

  for (int i = 1; i <= sam->getNoNodes(); ++i)
    REQUIRE(sam->getEquation(i, 1) == i);
}


TEST_CASE("TestSAM.MixedBasis2P")
{
  SIM2D sim({1,1});
  sim.read("src/LinAlg/Test/refdata/sam_2D_2P.xinp");
  sim.preprocess();

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_mixedbasis_2P.ref");

  for (int i = 1; i <= sam->getNoNodes(); ++i)
    REQUIRE(sam->getEquation(i, 1) == i);
}


TEST_CASE("TestSAM.SingleBasisDirichlet2P")
{
  SIM2D sim(1);
  sim.read("src/LinAlg/Test/refdata/sam_2D_dir_2P.xinp");
  sim.preprocess();

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_singlebasis_dir_2P.ref");

  REQUIRE(sam->getEquation(1, 1) == 1);
  REQUIRE(sam->getEquation(2, 1) == 2);
  REQUIRE(sam->getEquation(3, 1) == 3);
  REQUIRE(sam->getEquation(4, 1) == 4);
  REQUIRE(sam->getEquation(5, 1) == 0);
  REQUIRE(sam->getEquation(6, 1) == 0);
}


TEST_CASE("TestSAM.MixedBasisDirichlet2P")
{
  SIM2D sim({1,1});
  sim.read("src/LinAlg/Test/refdata/sam_2D_dir_2P.xinp");
  sim.preprocess();

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_mixedbasis_dir_2P.ref");

  int i, eq = 1;
  for (i = 1; i <= 13; ++i)
    REQUIRE(sam->getEquation(i, 1) == eq++);

  const ASMbase* pch = sim.getPatch(2);
  for (i = 1; i <= 3; ++i) {
    REQUIRE(sam->getEquation(pch->getNodeID((i-1)*3+2), 1) == eq++);
    REQUIRE(sam->getEquation(pch->getNodeID((i-1)*3+3), 1) == 0);
  }
  REQUIRE(sam->getEquation(20, 1) == eq++);
  REQUIRE(sam->getEquation(21, 1) == eq++);
}
