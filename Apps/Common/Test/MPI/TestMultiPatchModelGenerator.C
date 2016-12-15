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

#include "gtest/gtest.h"
#include "IFEM.h"
#include "IntegrandBase.h"
#include "SIMMultiPatchModelGen.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include <fstream>

typedef std::vector<int> IntVec;

template<class Dim>
class DummySIM : public SIMMultiPatchModelGen<Dim> {
public:
  class DummyIntegrand : public IntegrandBase {};
  DummySIM() : SIMMultiPatchModelGen<Dim>(1)
  { Dim::myProblem = new DummyIntegrand; }
};


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


TEST(TestMultiPatchModelGenerator2D, SubdivisionsMPI)
{
  DummySIM<SIM2D> sim;
  ASSERT_TRUE(sim.read("Test/refdata/modelgen2d_subdivision.xinp"));
  ASSERT_TRUE(sim.preprocess());
  const ProcessAdm& adm = sim.getProcessAdm();
  std::stringstream str;
  str << "Test/refdata/modelgen2d_subdivision_nodes" << adm.getProcId() << ".ref";
  IntVec B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGN(), B);
  str.str("");
  str << "Test/refdata/modelgen2d_subdivision_eqs" << adm.getProcId() << ".ref";
  B = readIntVector(str.str());
  check_intvectors_equal(adm.dd.getMLGEQ(), B);
}

/*
TES(TestMultiPatchModelGenerator3D, SubdivisionsMPI)
{
  DummySIM<SIM3D> sim;
  ASSERT_TRUE(sim.read("Test/refdata/modelgen3d_subdivision.xinp"));
  ASSERT_TRUE(sim.preprocess());
  const ProcessAdm& adm = sim.getProcessAdm();
  for (auto& it : adm.dd.getMLGN())
    IFEM::cout << it <<  " ";
  IFEM::cout << std::endl;
//  std::stringstream str;
//  str << "Test/refdata/modelgen3d_subdivision_nodes" << adm.getProcId() << ".ref";
//  IntVec B = readIntVector(str.str());
//  check_intvectors_equal(adm.dd.getMLGN(), B);
//  str.str("");
//  str << "Test/refdata/modelgen2d_subdivision_eqs" << adm.getProcId() << ".ref";
//  B = readIntVector(str.str());
//  check_intvectors_equal(adm.dd.getMLGEQ(), B);
}
*/
