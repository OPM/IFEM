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

#include "Catch2Support.h"

#include <fstream>
#include <numeric>


namespace
{
  using IntMat  = std::vector<IntVec>;
  using CharVec = std::vector<unsigned char>;


std::istream& operator>> (std::istream& is, IntVec& v)
{
  size_t nv;
  is >> nv;
  v.resize(nv);
  for (size_t i = 0; i < nv; i++)
    is >> v[i];

  return is;
}


IntMat readIntMatrix (size_t r, const std::string& file)
{
  IntMat result(r);
  std::ifstream f(file);
  for (size_t i = 0; i < r; i++)
    f >> result[i];

  return result;
}


void check_intmatrices_equal (const SAM* sam, const std::string& path)
{
  std::vector<IntSet> A;
  sam->getDofCouplings(A);

  size_t i, j;
  IntSet::const_iterator it;
  IntMat B = readIntMatrix(A.size(), path);
  for (i = 0; i < A.size(); i++)
    for (it = A[i].begin(), j = 0; it != A[i].end(); ++it, j++)
      CHECK(*it == B[i][j]);
};


  class MySIM : public SIM2D
  {
  public:
    template<class T> MySIM (const char* infile, const T& nf) : SIM2D(nf)
    {
      REQUIRE(this->read(infile));
      REQUIRE(this->preprocess());
    }
  };
}


TEST_CASE("TestSAM.SingleBasis1P")
{
  MySIM sim("src/LinAlg/Test/refdata/sam_2D_1P.xinp",1);

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_singlebasis_1P.ref");

  for (int i = 1; i <= sam->getNoNodes(); ++i)
    CHECK(sam->getEquation(i, 1) == i);
}


TEST_CASE("TestSAM.SingleBasisDirichlet1P")
{
  MySIM sim("src/LinAlg/Test/refdata/sam_2D_dir_1P.xinp",1);

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_singlebasis_dir_1P.ref");

  int i, eq = 1;
  for (i = 1; i <= 3; ++i) {
    CHECK(sam->getEquation(3*(i-1)+1, 1) == eq++);
    CHECK(sam->getEquation(3*(i-1)+2, 1) == eq++);
    CHECK(sam->getEquation(3*(i-1)+3, 1) == 0);
  }
}


TEST_CASE("TestSAM.MixedBasis1P")
{
  ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS1;
  ASMmxBase::itgBasis = 2;

  MySIM sim("src/LinAlg/Test/refdata/sam_2D_1P.xinp",CharVec{1,1});

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_mixedbasis_1P.ref");

  for (int i = 1; i <= sam->getNoNodes(); ++i)
    CHECK(sam->getEquation(i, 1) == i);
}


TEST_CASE("TestSAM.MixedBasisDirichlet1P")
{
  MySIM sim("src/LinAlg/Test/refdata/sam_2D_dir_1P.xinp",CharVec{1,1});

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_mixedbasis_dir_1P.ref");

  int i, eq = 1;
  for (i = 1; i <= 4; ++i) {
    CHECK(sam->getEquation(4*(i-1)+1, 1) == eq++);
    CHECK(sam->getEquation(4*(i-1)+2, 1) == eq++);
    CHECK(sam->getEquation(4*(i-1)+3, 1) == eq++);
    CHECK(sam->getEquation(4*(i-1)+4, 1) == 0);
  }

  for (i = 17; i <= 25; ++i)
    CHECK(sam->getEquation(i, 1) == eq++);
}


TEST_CASE("TestSAM.SingleBasis2P")
{
  MySIM sim("src/LinAlg/Test/refdata/sam_2D_2P.xinp",1);

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_singlebasis_2P.ref");

  for (int i = 1; i <= sam->getNoNodes(); ++i)
    CHECK(sam->getEquation(i, 1) == i);
}


TEST_CASE("TestSAM.MixedBasis2P")
{
  MySIM sim("src/LinAlg/Test/refdata/sam_2D_2P.xinp",CharVec{1,1});

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_mixedbasis_2P.ref");

  for (int i = 1; i <= sam->getNoNodes(); ++i)
    CHECK(sam->getEquation(i, 1) == i);
}


TEST_CASE("TestSAM.SingleBasisDirichlet2P")
{
  MySIM sim("src/LinAlg/Test/refdata/sam_2D_dir_2P.xinp",1);

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_singlebasis_dir_2P.ref");

  CHECK(sam->getEquation(1, 1) == 1);
  CHECK(sam->getEquation(2, 1) == 2);
  CHECK(sam->getEquation(3, 1) == 3);
  CHECK(sam->getEquation(4, 1) == 4);
  CHECK(sam->getEquation(5, 1) == 0);
  CHECK(sam->getEquation(6, 1) == 0);
}


TEST_CASE("TestSAM.MixedBasisDirichlet2P")
{
  MySIM sim("src/LinAlg/Test/refdata/sam_2D_dir_2P.xinp",CharVec{1,1});

  const SAM* sam = sim.getSAM();
  check_intmatrices_equal(sam, "src/LinAlg/Test/refdata/sam_2D_mixedbasis_dir_2P.ref");

  int i, eq = 1;
  for (i = 1; i <= 13; ++i)
    CHECK(sam->getEquation(i, 1) == eq++);

  const ASMbase* pch = sim.getPatch(2);
  for (i = 1; i <= 3; ++i) {
    CHECK(sam->getEquation(pch->getNodeID((i-1)*3+2), 1) == eq++);
    CHECK(sam->getEquation(pch->getNodeID((i-1)*3+3), 1) == 0);
  }
  CHECK(sam->getEquation(20, 1) == eq++);
  CHECK(sam->getEquation(21, 1) == eq++);
}


TEST_CASE("TestSAM.NodeAndLocalDof")
{
  for (unsigned char nf : {1, 2})
  {
    MySIM sim("src/LinAlg/Test/refdata/sam_2D_dir_1P.xinp",nf);

    const SAM* sam = sim.getSAM();
    int nEq = 0;
    for (int inod = 1; inod <= sam->getNoNodes(); inod++)
      for (int ldof = 1; ldof <= nf; ldof++)
        if (int ieq = sam->getEquation(inod,ldof); ieq > 0)
        {
          ++nEq;
          CHECK(sam->getNodeAndLocalDof(ieq,true) == std::make_pair(inod,ldof));
        }

    CHECK(nEq == sam->getNoEquations());
  }
}


TEST_CASE("TestSAM.Norms")
{
  MySIM sim("src/LinAlg/Test/refdata/sam_2D_1P.xinp",3);

  const SAM* sam = sim.getSAM();

  const size_t neq = sam->getNoEquations();
  const size_t nno = sam->getNoNodes();
  REQUIRE(neq == 3*nno);

  Vector x(neq),y(neq);
  std::iota(x.begin(),x.end(),1.0);

  size_t comp = 1;
  double refn = neq-2;
  double norm = sam->normInf(x,comp);
  CHECK(norm == refn);
  CHECK(comp == nno);

  comp = 2;
  refn = neq-1;
  norm = sam->normInf(x,comp);
  CHECK(norm == refn);
  CHECK(comp == nno);

  comp = 3;
  refn = neq;
  norm = sam->normInf(x,comp);
  CHECK(norm == refn);
  CHECK(comp == nno);

  comp = 999;
  refn = neq*neq + (neq-1)*(neq-1) + (neq-2)*(neq-2);
  norm = sam->normInf(x,comp);
  CHECK(norm*norm == refn);
  CHECK(comp == nno);

  comp = 2;
  norm = sam->normInf(y,comp);
  CHECK(norm == 0.0);
  CHECK(comp == 1);

  comp = 999;
  norm = sam->normInf(y,comp);
  CHECK(norm == 0.0);
  CHECK(comp == 1);
}
