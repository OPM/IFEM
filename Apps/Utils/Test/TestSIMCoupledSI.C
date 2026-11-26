//==============================================================================
//!
//! \file TestSIMCoupledSI.C
//!
//! \date Nov 26 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for the semi-implicit coupling of two simulators.
//!
//==============================================================================

#include "SIMCoupledSI.h"
#include "ProcessAdm.h"

#include "Catch2Support.h"

#include <tinyxml2.h>
#include <vector>


namespace {

//! \brief A simulator of the coupling, whose iterate of each cycle is -1.
class MockSolver
{
public:
  //! \brief Solves a cycle, giving the iterate.
  SIM::ConvStatus solveIteration(TimeStep&)
  {
    sol(1) = -1.0;
    return SIM::OK;
  }
  //! \brief Post-processes the solution of the step.
  bool postSolve(const TimeStep&) { return true; }
  //! \brief Returns the parallel process administrator.
  const ProcessAdm& getProcessAdm() const { return adm; }
  //! \brief Returns no equation system administrator.
  const SAM* getSAM() const { return nullptr; }

  // Stubs of the methods of the coupling
  bool advanceStep(TimeStep&) { return true; } //!< Stub
  bool solveStep(TimeStep&) { return true; } //!< Stub
  bool saveStep(const TimeStep&, int&) { return true; } //!< Stub
  bool saveModel(char*, int&, int&) { return true; } //!< Stub
  bool init(const TimeStep&) { return true; } //!< Stub
  void setVTF(VTF*) {} //!< Stub
  VTF* getVTF() const { return nullptr; } //!< Stub
  //! \brief Stub.
  void registerDependency(SIMdependency*, const std::string&, short int,
                          const std::vector<ASMbase*>&, char, int) {}
  //! \brief Stub.
  void registerDependency(SIMdependency*, const std::string&, short int) {}

  Vector sol = Vector(1); //!< The solution
  ProcessAdm adm;    //!< The process administrator
};


//! \brief A coupling relaxing the solution of the second simulator, which
//! clamps the relaxed solution at zero, as a model may clamp its unknowns.
class Coupling : public SIMCoupledSI<MockSolver,MockSolver>
{
public:
  //! \brief The constructor sets the relaxation weight and the cycles.
  Coupling(MockSolver& s1, MockSolver& s2, int cycles) :
    SIMCoupledSI<MockSolver,MockSolver>(s1,s2), nCycles(cycles)
  {
    tinyxml2::XMLDocument doc;
    doc.Parse("<subiterations relax='0.5'/>");
    this->parseIterations(doc.RootElement());
  }

  //! \brief Converges after the given number of cycles.
  SIM::ConvStatus checkConvergence(const TimeStep& tp, SIM::ConvStatus,
                                   SIM::ConvStatus) override
  {
    return tp.iter+1 < nCycles ? SIM::OK : SIM::CONVERGED;
  }

  //! \brief Returns the solution of the second simulator.
  const Vector& getRelaxationVector() const override { return S2.sol; }

  //! \brief Sets the relaxed solution, clamped at zero.
  void setRelaxedSolution(const Vector& sol) override
  {
    relaxed.push_back(sol(1));
    S2.sol(1) = std::max(sol(1),0.0);
  }

  int nCycles; //!< The number of cycles
  std::vector<double> relaxed; //!< The relaxed solutions of the cycles
};

}


TEST_CASE("TestSIMCoupledSI.ClampedRelaxation")
{
  // The relaxed solution of the second cycle, the mean of -1 and -1, is
  // clamped at zero by the simulator, and the next cycle relaxes from that
  MockSolver s1, s2;
  Coupling coupling(s1, s2, 4);
  TimeStep tp;
  REQUIRE(coupling.solveStep(tp));

  // The first cycle is not relaxed, and the converged one is the solution
  REQUIRE(coupling.relaxed.size() == 2);
  REQUIRE(coupling.relaxed[0] == -1.0);
  REQUIRE(coupling.relaxed[1] == -0.5);
}
