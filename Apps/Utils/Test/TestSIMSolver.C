// $Id$
//==============================================================================
//!
//! \file TestSIMSolver.C
//!
//! \date Oct 8 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for the time stepping loop of the SIM solver template.
//!
//==============================================================================

#include "SIMSolver.h"
#include "TimeStep.h"

#include "Catch2Support.h"
#include "tinyxml2.h"

#include <cstring>
#include <string>
#include <vector>


namespace {

/*!
  \brief A solver integrating du/dt = 1, which fails for time steps larger
  than 0.1 ending in the interval (0.3,0.6), after corrupting its state.
*/

class FailingSolver
{
public:
  bool saveModel(char*, int&, int&) { return true; }
  bool saveStep(const TimeStep& tp, int&) { times.push_back(tp.time.t); return true; }
  bool advanceStep(TimeStep&) { ++advances; return true; }

  bool solveStep(TimeStep& tp)
  {
    const double t = tp.time.t, dt = tp.time.dt;
    u += dt;
    if (dt <= 0.1 || t <= 0.3 || t - dt >= 0.6)
      return true;

    u += 1000.0; // The failing step leaves a corrupt state
    ++failures;
    return false;
  }

  bool serialize(std::map<std::string,std::string>& data) const
  {
    if (!canSerialize)
      return false;
    data["u"] = std::string(reinterpret_cast<const char*>(&u),sizeof(u));
    return true;
  }

  bool deSerialize(const std::map<std::string,std::string>& data)
  {
    std::memcpy(&u,data.at("u").data(),sizeof(u));
    return true;
  }

  double u = 0.0;           //!< The solution
  int advances = 0;         //!< Number of advanced steps
  int failures = 0;         //!< Number of failed steps
  bool canSerialize = true; //!< If \e false, the state can not be serialized
  std::vector<double> times; //!< The times of the saved steps
};


//! \brief A SIM solver parsing its time stepping from an XML string.
class TestSolver : public SIMSolver<FailingSolver>
{
public:
  explicit TestSolver(FailingSolver& s) : SIMSolver<FailingSolver>(s) {}

  void parseXML(const char* xml)
  {
    tinyxml2::XMLDocument doc;
    REQUIRE(doc.Parse(xml) == tinyxml2::XML_SUCCESS);
    REQUIRE(this->parse(doc.RootElement()));
  }
};

const char* withCutback =
  "<timestepping maxCutbacks=\"3\">"
  "  <step start=\"0\" end=\"1.0\">0.25</step>"
  "</timestepping>";

}


TEST_CASE("TestSIMSolver.Cutback")
{
  FailingSolver solver;
  TestSolver sim(solver);
  sim.parseXML(withCutback);
  REQUIRE(sim.solveProblem(nullptr) == 0);

  // The failing steps are restarted from the state before them,
  // and the initial state is saved in addition to each step
  REQUIRE(solver.failures > 0);
  REQUIRE(solver.advances == static_cast<int>(solver.times.size()) - 1 +
                             solver.failures);
  REQUIRE_THAT(solver.u, WithinRel(1.0));
  REQUIRE_THAT(solver.times.back(), WithinRel(1.0));
  for (size_t i = 1; i < solver.times.size(); i++)
    REQUIRE(solver.times[i] > solver.times[i-1]);
}


TEST_CASE("TestSIMSolver.NoCutback")
{
  FailingSolver solver;
  TestSolver sim(solver);
  sim.parseXML("<timestepping>"
               "  <step start=\"0\" end=\"1.0\">0.25</step>"
               "</timestepping>");
  REQUIRE(sim.solveProblem(nullptr) == 3);
  REQUIRE(solver.failures == 1);
}


TEST_CASE("TestSIMSolver.NoSerialization")
{
  // Without a serialized state the time steps are not cut back
  FailingSolver solver;
  solver.canSerialize = false;
  TestSolver sim(solver);
  sim.parseXML(withCutback);
  REQUIRE(sim.solveProblem(nullptr) == 3);
  REQUIRE(solver.failures == 1);
}


TEST_CASE("TestSIMSolver.MaxCutbacks")
{
  // One cut-back by f2 = 0.5 is not enough to get below 0.1
  FailingSolver solver;
  TestSolver sim(solver);
  sim.parseXML("<timestepping maxCutbacks=\"1\" f2=\"0.5\">"
               "  <step start=\"0\" end=\"1.0\">0.25</step>"
               "</timestepping>");
  REQUIRE(sim.solveProblem(nullptr) == 3);
  REQUIRE(solver.failures == 2);
}
