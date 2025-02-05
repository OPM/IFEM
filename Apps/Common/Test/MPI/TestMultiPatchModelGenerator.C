//==============================================================================
//!
//! \file TestMultiPatchModelGenerator.C
//!
//! \date Dec 14 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for multi-patch model generators with subdivision.
//!
//==============================================================================

#include "SIMMultiPatchModelGen.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include <fstream>

#include "gtest/gtest.h"


static void check_intvectors_equal (const std::vector<int>& A,
                                    const std::string& Bfile)
{
  std::ifstream f(Bfile);

  size_t Bsize;
  f >> Bsize;
  if (!f) Bsize = 0;

  ASSERT_EQ(A.size(),Bsize);
  std::vector<int> B(Bsize);
  for (int& bi : B) f >> bi;
  for (size_t i = 0; i < A.size(); i++)
    EXPECT_EQ(A[i],B[i]);
};


template<class Dim> class SIMMultiPatch : public SIMMultiPatchModelGen<Dim>
{
public:
  explicit SIMMultiPatch(int n1) : SIMMultiPatchModelGen<Dim>(n1,false) {}
};


TEST(TestMultiPatchModelGenerator2D, SubdivisionsMPI)
{
  SIMMultiPatch<SIM2D> sim(1);
  ASSERT_TRUE(sim.read("Test/refdata/modelgen2d_subdivision.xinp"));
  ASSERT_TRUE(sim.preprocess());
  const ProcessAdm& adm = sim.getProcessAdm();
  std::stringstream str;
  str << "Test/refdata/modelgen2d_subdivision_nodes" << adm.getProcId() << ".ref";
  check_intvectors_equal(adm.dd.getMLGN(), str.str());
  str.str("");
  str << "Test/refdata/modelgen2d_subdivision_eqs" << adm.getProcId() << ".ref";
  check_intvectors_equal(adm.dd.getMLGEQ(), str.str());
}

/*
TES(TestMultiPatchModelGenerator3D, SubdivisionsMPI)
{
  SIMMultiPatch<SIM3D> sim(1);
  ASSERT_TRUE(sim.read("Test/refdata/modelgen3d_subdivision.xinp"));
  ASSERT_TRUE(sim.preprocess());
  const ProcessAdm& adm = sim.getProcessAdm();
  for (int i : adm.dd.getMLGN()) std::cout <<" "<< i;
  std::cout << std::endl;
//  std::stringstream str;
//  str << "Test/refdata/modelgen3d_subdivision_nodes" << adm.getProcId() << ".ref";
//  check_intvectors_equal(adm.dd.getMLGN(), str.str());
//  str.str("");
//  str << "Test/refdata/modelgen2d_subdivision_eqs" << adm.getProcId() << ".ref";
//  check_intvectors_equal(adm.dd.getMLGEQ(), str.str());
}
*/
