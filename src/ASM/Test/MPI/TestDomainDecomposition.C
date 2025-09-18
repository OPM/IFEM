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
#include "readIntVec.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

namespace {

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
  REQUIRE(A.size() == B.size());
  for (size_t i = 0; i < A.size(); i++)
    REQUIRE(A[i] == B[i]);
}

}


TEST_CASE("TestDomainDecomposition.Corner2D")
{
  const int param = GENERATE(0,1);

  SECTION("Corner " + std::to_string(param)) {
    std::string file = "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_corner";
    if (param == 1)
      file += "_fail.xinp";
    else
      file += ".xinp";

    SIM2D sim;
    REQUIRE(sim.read(file.c_str()));

    if (param == 0) {
      REQUIRE(sim.preprocess());
      const ProcessAdm& adm = sim.getProcessAdm();
      std::stringstream str;
      str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_corner_nodes";
      str << adm.getProcId() << ".ref";
      IntVec B = readIntVector(str.str());
      check_intvectors_equal(adm.dd.getMLGN(), B);
    } else
      REQUIRE(!sim.preprocess());
  }
}


TEST_CASE("TestDomainDecomposition.Corner2DLR")
{
  const int param = GENERATE(0,1);

  SECTION("Corner " + std::to_string(param)) {
    std::stringstream str;
    std::string file = "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_corner";
    if (param == 1)
      file += "_fail.xinp";
    else
      file += ".xinp";

    SIM2D sim;
    sim.opt.discretization = ASM::LRSpline;
    REQUIRE(sim.read(file.c_str()));

    if (param == 0) {
      REQUIRE(sim.preprocess());
      const ProcessAdm& adm = sim.getProcessAdm();
      std::stringstream str;
      str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_cornerLR_nodes";
      str << adm.getProcId() << ".ref";
      IntVec B = readIntVector(str.str());
      check_intvectors_equal(adm.dd.getMLGN(), B);
    } else
      REQUIRE(!sim.preprocess());
  }
}


TEST_CASE("TestDomainDecomposition.SetupSingleBasis2D")
{
  const int param = GENERATE(0,1);

  SECTION("Orient " + std::to_string(param)) {
    SIM2D sim(2);
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
    str << param << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.preprocess());

    const ProcessAdm& adm = sim.getProcessAdm();
    const SAM* sam = sim.getSAM();
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
    str << param << "_nodes" << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGN(), B);
    REQUIRE(sam->getNoNodes() == 9);
    const std::vector<int> maxeqs {16, 26, 36, 44};
    REQUIRE(adm.dd.getMaxEq() == maxeqs[adm.getProcId()]);
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
    str << param << "_eqs" << adm.getProcId() << ".ref";
    B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGEQ(), B);
    REQUIRE(adm.dd.getMaxDOF() == 2*adm.dd.getMaxNode());
  }
}


TEST_CASE("TestDomainDecomposition.SetupSingleBasis2DLR")
{
  const int param = GENERATE(0,1);

  SECTION("Orient " + std::to_string(param)) {
    SIM2D sim(2);
    sim.opt.discretization = ASM::LRSpline;
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
    str << param << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.preprocess());

    const ProcessAdm& adm = sim.getProcessAdm();
    const SAM* sam = sim.getSAM();
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_LR_orient";
    str << param << "_nodes" << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGN(), B);
    REQUIRE(sam->getNoNodes() == 9);
    const std::vector<int> maxeqs {16, 26, 36, 44};
    REQUIRE(adm.dd.getMaxEq() == maxeqs[adm.getProcId()]);
  }
}


TEST_CASE("TestDomainDecomposition.SetupSingleBasisPeriodic2D")
{
  const int param = GENERATE(0,1);

  SECTION("Orient " + std::to_string(param)) {
    SIM2D sim(2);
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_periodic";
    str << param << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.preprocess());

    const ProcessAdm& adm = sim.getProcessAdm();
    const SAM* sam = sim.getSAM();
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_periodic";
    str << param << "_nodes" << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGN(), B);
    REQUIRE(sam->getNoNodes() == 9);
  }
}


TEST_CASE("TestDomainDecomposition.SetupSingleBasisGlobalLM2D")
{
  const int param = GENERATE(0,1);

  SECTION("Orient " + std::to_string(param)) {
    TestGlobalLMSIM<SIM2D> sim(2);
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
    str << param << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.preprocess());

    const ProcessAdm& adm = sim.getProcessAdm();
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
    str << param<< "_eqs" << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    if (adm.getProcId() == 0)
      B.push_back(17);
    else {
      for (int& it : B)
        if (it > 16)
          ++it;
      B.push_back(17);
    }
    check_intvectors_equal(adm.dd.getMLGEQ(), B);
  }
}


TEST_CASE("TestDomainDecomposition.SetupSingleBasisBlockEqsComponent2D")
{
  const int param = GENERATE(0,1);

  SECTION("Orient " + std::to_string(param)) {
    SIM2D sim(3);
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_blocks_components_orient";
    str << param << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.preprocess());

    const ProcessAdm& adm = sim.getProcessAdm();
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
    str << param << "_eqs" << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGEQ(1), B);
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_blocks_components_orient";
    str << param << "_block2_eqs" << adm.getProcId() << ".ref";
    B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGEQ(2), B);
  }
}


TEST_CASE("TestDomainDecomposition.SetupMixedBasis2D")
{
  const int param = GENERATE(0,1);

  SECTION("Orient " + std::to_string(param)) {
    SIM2D sim({2,2});
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
    str << param << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.preprocess());

    const ProcessAdm& adm = sim.getProcessAdm();
    const SAM* sam = sim.getSAM();
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_mixed_orient";
    str << param << "_nodes" << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGN(), B);
    REQUIRE(sam->getNoNodes() == 25);
    REQUIRE(sam->getNoDOFs() == 50);
    const std::vector<int> maxeqs {48, 81, 114, 140};
    REQUIRE(adm.dd.getMaxEq() == maxeqs[adm.getProcId()]);
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_mixed_orient";
    str << param << "_eqs" << adm.getProcId() << ".ref";
    B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGEQ(), B);
    REQUIRE(adm.dd.getMinDOF() == 2*(adm.dd.getMinNode()-1)+1);
    REQUIRE(adm.dd.getMaxDOF() == 2*adm.dd.getMaxNode());
  }
}


TEST_CASE("TestDomainDecompositionSetupMixedBasis2DLR")
{
  const int param = GENERATE(0,1);

  SECTION("Orient " + std::to_string(param)) {
    SIM2D sim({2,2});
    sim.opt.discretization = ASM::LRSpline;
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
    str << param << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.preprocess());

    const ProcessAdm& adm = sim.getProcessAdm();
    const SAM* sam = sim.getSAM();
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_mixed_LR_orient";
    str << param << "_nodes" << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGN(), B);
    REQUIRE(sam->getNoNodes() == 25);
    REQUIRE(sam->getNoDOFs() == 50);
    const std::vector<int> maxeqs {48, 81, 114, 140};
    REQUIRE(adm.dd.getMaxEq() == maxeqs[adm.getProcId()]);
  }
}


TEST_CASE("TestDomainDecomposition.SetupMixedBasisPeriodic2D")
{
  const int param = GENERATE(0,1);

  SECTION("Orient " + std::to_string(param)) {
    SIM2D sim({2,2});
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_periodic";
    str << param << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.preprocess());

    const ProcessAdm& adm = sim.getProcessAdm();
    const SAM* sam = sim.getSAM();
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_mixed_periodic";
    str << param << "_nodes" << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGN(), B);
    REQUIRE(sam->getNoNodes() == 25);
    REQUIRE(sam->getNoDOFs() == 50);
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_mixed_periodic";
    str << param << "_eqs" << adm.getProcId() << ".ref";
    B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGEQ(), B);
  }
}


TEST_CASE("TestDomainDecomposition.SetupMixedBasisBlockEqsBasis2D")
{
  const int param = GENERATE(0,1);

  SECTION("Orient " + std::to_string(param)) {
    SIM2D sim({2,2});
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_blocks_basis_orient";
    str << param << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.preprocess());

    const ProcessAdm& adm = sim.getProcessAdm();
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_mixed_orient";
    str << param << "_eqs" << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGEQ(0), B);
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_blocks_basis_orient";
    str << param << "_block1_eqs" << adm.getProcId() << ".ref";
    B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGEQ(1), B);
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_blocks_basis_orient";
    str << param << "_block2_eqs" << adm.getProcId() << ".ref";
    B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGEQ(2), B);
  }
}


TEST_CASE("TestDomainDecomposition.SetupMixedBasisBlockEqsBasisGlobalLM2D")
{
  const int param = GENERATE(0,1);

  SECTION("Orient " + std::to_string(param)) {
    TestGlobalLMSIM<SIM2D> sim({2,2});
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_blocks_basis_orient";
    str << param << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.preprocess());

    const ProcessAdm& adm = sim.getProcessAdm();
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_mixed_orient";
    str << param << "_eqs" << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    if (adm.getProcId() == 0)
      B.push_back(49);
    else {
      for (int& it : B)
        if (it > 48)
          ++it;
      B.push_back(49);
    }
    check_intvectors_equal(adm.dd.getMLGEQ(0), B);
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_blocks_basis_orient";
    str << param << "_block1_eqs" << adm.getProcId() << ".ref";
    B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGEQ(1), B);
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_blocks_basis_orient";
    str << param << "_block2_eqs" << adm.getProcId() << ".ref";
    B = readIntVector(str.str());
    if (adm.getProcId() == 0)
      B.push_back(19);
    else {
      for (int& it : B)
        if (it > 18)
          ++it;
      B.push_back(19);
    }
    check_intvectors_equal(adm.dd.getMLGEQ(2), B);
  }
}


TEST_CASE("TestDomainDecomposition.SetupSingleBasis3D")
{
  const int param = GENERATE(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);

  SECTION("Orient " + std::to_string(param)) {
    SIM3D sim(3);
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
    str << param << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.preprocess());

    const ProcessAdm& adm = sim.getProcessAdm();
    const SAM* sam = sim.getSAM();
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
    str << param << "_nodes" << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGN(), B);
    REQUIRE(sam->getNoNodes() == 45);
    const std::vector<int> maxeqs {111, 195, 254, 300};
    REQUIRE(adm.dd.getMaxEq() == maxeqs[adm.getProcId()]);
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
    str << param << "_eqs" << adm.getProcId() << ".ref";
    B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGEQ(), B);
    REQUIRE(adm.dd.getMaxDOF() == 3*adm.dd.getMaxNode());
  }
}


TEST_CASE("TestDomainDecomposition.SetupSingleBasisLR3D")
{
  const int param = GENERATE(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);

  SECTION("Orient " + std::to_string(param)) {
    SIM3D sim(3);
    sim.opt.discretization = ASM::LRSpline;
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
    str << param << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.preprocess());

    const ProcessAdm& adm = sim.getProcessAdm();
    const SAM* sam = sim.getSAM();
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_LR_orient";
    str << param << "_nodes" << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGN(), B);
    REQUIRE(sam->getNoNodes() == 45);
    const std::vector<int> maxeqs {111, 195, 254, 300};
    REQUIRE(adm.dd.getMaxEq() == maxeqs[adm.getProcId()]);
  }
}


TEST_CASE("TestDomainDecompositionSetupSingleBasisPeriodic3D")
{
  const int param = GENERATE(0,1,2);

  SECTION("Orient " + std::to_string(param)) {
    SIM3D sim(3);
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_periodic";
    str << param << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.preprocess());

    const ProcessAdm& adm = sim.getProcessAdm();
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_periodic";
    str << param << "_nodes" << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGN(), B);
  }
}


TEST_CASE("TestDomainDecomposition.SetupSingleBasisBlockEqsComponent3D")
{
  const int param = 0;

  SECTION("Orient " + std::to_string(param)) {
    SIM3D sim(4);
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_blocks_components_orient";
    str << param << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.preprocess());

    const ProcessAdm& adm = sim.getProcessAdm();
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
    str << param << "_eqs" << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGEQ(1), B);
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
    str << param << "_nodes" << adm.getProcId() << ".ref";
    B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGEQ(2), B);
  }
}


TEST_CASE("TestDomainDecomposition.SetupMixedBasis3D")
{
  const int param = 0;

  SECTION("Orient " + std::to_string(param)) {
    SIM3D sim({3,1});
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
    str << param << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.preprocess());

    const ProcessAdm& adm = sim.getProcessAdm();
    const SAM* sam = sim.getSAM();
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_mixed_orient";
    str << param << "_nodes" << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGN(), B);
    REQUIRE(sam->getNoNodes() == 157);
    const std::vector<int> maxeqs {337, 607, 828, 1007};
    REQUIRE(adm.dd.getMaxEq() == maxeqs[adm.getProcId()]);
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_mixed_orient";
    str << param << "_eqs" << adm.getProcId() << ".ref";
    B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGEQ(), B);
  }
}


TEST_CASE("TestDomainDecomposition.SetupMixedBasisPeriodic3D")
{
  const int param = GENERATE(0,1,2);

  SECTION("Orient " + std::to_string(param)) {
    SIM3D sim({1,1});
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_periodic";
    str << param << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.preprocess());

    const ProcessAdm& adm = sim.getProcessAdm();
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_mixed_periodic";
    str << param << "_nodes" << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());

    check_intvectors_equal(adm.dd.getMLGN(), B);
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_mixed_periodic";
    str << param << "_eqs" << adm.getProcId() << ".ref";
    B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGEQ(), B);
  }
}


TEST_CASE("TestDomainDecomposition.SetupMixedBasisPeriodicLM3D")
{
  const int param = 0;

  SECTION("Orient " + std::to_string(param)) {
    TestGlobalLMSIM<SIM3D> sim({1,1});
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_periodiclm";
    str << param << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));
    REQUIRE(sim.preprocess());

    const ProcessAdm& adm = sim.getProcessAdm();
    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_mixed_periodiclm";
    str << param << "_nodes" << adm.getProcId() << ".ref";
    IntVec B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGN(), B);

    str.str("");
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_mixed_periodiclm";
    str << param << "_eqs" << adm.getProcId() << ".ref";
    B = readIntVector(str.str());
    check_intvectors_equal(adm.dd.getMLGEQ(), B);
  }
}


TEST_CASE("TestDomainDecomposition.Corner3D")
{
  const int param = GENERATE(0,1);

  SECTION(param == 0 ? "Success" : "Failure") {
    SIM3D sim;
    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_corner";
    if (param > 0)
      str << "_fail";
    str << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));

    if (param < 1) {
      REQUIRE(sim.preprocess());
      const ProcessAdm& adm = sim.getProcessAdm();
      str.str("");
      str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_corner_nodes";
      str << adm.getProcId() << ".ref";
      IntVec B = readIntVector(str.str());
      check_intvectors_equal(adm.dd.getMLGN(), B);
    } else
      REQUIRE(!sim.preprocess());
  }
}


TEST_CASE("TestDomainDecomposition.Corner3DLR")
{
  const int param = GENERATE(0,1);

  SECTION(param == 0 ? "Success" : "Failure") {
    SIM3D sim;
    sim.opt.discretization = ASM::LRSpline;

    std::stringstream str;
    str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_corner";
    if (param > 0)
      str << "_fail";
    str << ".xinp";
    REQUIRE(sim.read(str.str().c_str()));

    if (param < 1) {
      REQUIRE(sim.preprocess());
      const ProcessAdm& adm = sim.getProcessAdm();
      str.str("");
      str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_cornerLR_nodes";
      str << adm.getProcId() << ".ref";
      IntVec B = readIntVector(str.str());
      check_intvectors_equal(adm.dd.getMLGN(), B);
    } else
      REQUIRE(!sim.preprocess());
  }
}
