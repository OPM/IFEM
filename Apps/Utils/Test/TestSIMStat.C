//==============================================================================
//!
//! \file TestSIMStat.C
//!
//! \date Oct 16 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for the stationary driver.
//!
//==============================================================================

#include "SIMStat.h"

#include "SAM.h"
#include "SIM2D.h"
#include "SIMadmin.h"
#include "SIMsolution.h"
#include "TimeStep.h"

#include "Catch2Support.h"

#include <string>
#include <vector>


namespace {

//! \brief A time stepping simulator whose step is an affine contraction.
//! \details Shaped like the simulators the driver is meant for, a Dim with
//! the solution vectors mixed in. The step is a map with a fixed point of
//! its own rather than a scaling, since a scaling towards zero never
//! settles in the relative measure the driver uses. The two fields approach
//! that point at very different rates, which is what a saddle point problem
//! does and what the driver has to be told about.
class StepSIM : public SIM2D, public SIMsolution
{
public:
  //! \brief The constructor sets up the solution vectors and the two fields.
  //! \details No mesh is preprocessed here, so the assembly management data
  //! the driver takes its norms through is an empty one. Its inner product
  //! is that of the vectors, there being no shared degrees of freedom.
  StepSIM()
  {
    mySam = new SAM();
    this->initSolution(nDofs,2);
    for (int d = 1; d <= nDofs; d++)
      (d <= nDofs/2 ? fastDofs : slowDofs).push_back(d);
  }

  //! \brief Takes one step of the contraction.
  bool solveTimeStep(const TimeStep&)
  {
    ++nSteps;
    Vector& x = solution.front();
    for (int d : fastDofs) x(d) = 0.5*x(d) + 0.5;
    for (int d : slowDofs) x(d) = 0.9*x(d) + 0.1;
    return true;
  }

  //! \brief Puts the two solution vectors at different points.
  //! \details What a driver leaves in index one is the step before this
  //! call, which is not in general the point the iteration starts at.
  void startFrom(double now, double before)
  {
    solution[0].fill(now);
    solution[1].fill(before);
  }

  //! \brief Counts the post-processing, which is wanted once per solve.
  bool postSolve(const TimeStep&) { ++nPost; return true; }

  //! \brief Pushes this step's solution onto the one before it.
  bool advanceStep(const TimeStep&) { solution[1] = solution[0]; return true; }

  //! \brief Nothing of the model needs initializing here.
  bool init(const TimeStep&) { return true; }

  //! \brief Returns the two fields, the fast one first.
  std::vector<ConvGroup> convergenceGroups() const
  { return {{"fast", fastDofs}, {"slow", slowDofs}}; }

  //! \brief No model is set up, and the summary is not what is tested.
  void printSolutionSummary(const Vector&, int, const char*,
                            std::streamsize) override {}

  //! \brief The tag this simulator would read its own options from.
  std::string inputContext = "steps";

  int nSteps = 0; //!< Steps taken since the simulator was made
  int nPost = 0;  //!< Times the post-processing has run

  static constexpr int nDofs = 8; //!< Degrees of freedom, split over the two

private:
  std::vector<int> fastDofs; //!< The field which settles quickly
  std::vector<int> slowDofs; //!< The field which does not
};


//! \brief Keeps the convergence history off the console while a test runs.
//! \details The level is shared, so it is put back on the way out.
struct Quiet
{
  //! \brief The constructor silences the console.
  Quiet() : was(SIMadmin::msgLevel) { SIMadmin::msgLevel = -1; }
  //! \brief The destructor puts the level back.
  ~Quiet() { SIMadmin::msgLevel = was; }
  int was; //!< The level to restore
};


//! \brief The driver, with the knobs the input would otherwise turn.
class Stat : public SIMStat<StepSIM>
{
public:
  //! \brief Sets which field the change is measured on.
  void measureOn(const std::string& field) { cNorm = field; }
  //! \brief Sets the change accepted as stationary.
  void tolerance(double tol) { rTol = tol; }
  //! \brief Sets the iteration cap.
  void cap(int n) { maxIter = n; }
  //! \brief Turns the acceleration off, so the rate is the scheme's own.
  void unaccelerated() { anderson.setDepth(0); }
};

}


TEST_CASE("TestSIMStat.MeasuresTheFieldItIsTold")
{
  Quiet quiet;
  // Both fields reach the same point; the slow one takes far longer to get
  // there, so which is measured decides when the iteration stops
  TimeStep tp;
  Stat fast;
  fast.unaccelerated();
  fast.tolerance(1.0e-8);
  fast.measureOn("fast");
  REQUIRE(fast.solveStep(tp));

  Stat slow;
  slow.unaccelerated();
  slow.tolerance(1.0e-8);
  slow.measureOn("slow");
  REQUIRE(slow.solveStep(tp));

  INFO("measuring the slow field has to cost more iterations than the fast");
  REQUIRE(slow.nSteps > fast.nSteps);

  // and both arrive, the fixed point being one on every degree of freedom
  for (int d = 1; d <= StepSIM::nDofs; d++)
    REQUIRE_THAT(slow.SIMsolution::getSolution()(d), WithinAbs(1.0, 1.0e-6));
}


TEST_CASE("TestSIMStat.StartsFromThePointItIsAt")
{
  Quiet quiet;
  // The first iteration does not advance, so index one still holds the step
  // before this call rather than the point the iteration starts at. The step
  // this simulator takes does not look at index one, so what is left there
  // must not reach the iteration at all
  TimeStep tp;
  Stat same;
  same.startFrom(0.0,0.0);
  same.tolerance(1.0e-8);
  REQUIRE(same.solveStep(tp));

  Stat apart;
  apart.startFrom(0.0,-20.0);
  apart.tolerance(1.0e-8);
  REQUIRE(apart.solveStep(tp));

  INFO("what index one held must not change the iteration");
  REQUIRE(apart.nSteps == same.nSteps);

  for (int d = 1; d <= StepSIM::nDofs; d++)
    REQUIRE_THAT(apart.SIMsolution::getSolution()(d), WithinAbs(1.0, 1.0e-6));
}


TEST_CASE("TestSIMStat.AllWaitsForTheSlowestField")
{
  Quiet quiet;
  // Taking every field at once is as slow as its slowest part, which is the
  // reason the driver lets a simulator name them separately at all
  TimeStep tp;
  Stat all;
  all.unaccelerated();
  all.tolerance(1.0e-8);
  all.measureOn("all");
  REQUIRE(all.solveStep(tp));

  Stat slow;
  slow.unaccelerated();
  slow.tolerance(1.0e-8);
  slow.measureOn("slow");
  REQUIRE(slow.solveStep(tp));

  REQUIRE_THAT(static_cast<double>(all.nSteps),
               WithinRel(static_cast<double>(slow.nSteps), 0.1));
}


TEST_CASE("TestSIMStat.MeasuresThePrimaryFieldByDefault")
{
  Quiet quiet;
  // Told nothing, it measures the one the simulator puts first
  TimeStep tp;
  Stat byDefault;
  byDefault.unaccelerated();
  byDefault.tolerance(1.0e-8);
  REQUIRE(byDefault.solveStep(tp));

  Stat named;
  named.unaccelerated();
  named.tolerance(1.0e-8);
  named.measureOn("fast");
  REQUIRE(named.solveStep(tp));

  REQUIRE(byDefault.nSteps == named.nSteps);
}


TEST_CASE("TestSIMStat.RefusesAFieldTheSimulatorDoesNotHave")
{
  Quiet quiet;
  // Rather than measuring something else and converging on it
  TimeStep tp;
  Stat sim;
  sim.measureOn("temperature");
  REQUIRE_FALSE(sim.solveStep(tp));
  REQUIRE(sim.nSteps == 0);
}


TEST_CASE("TestSIMStat.PostProcessesOnceRatherThanPerIteration")
{
  Quiet quiet;
  // The iterations are of the pseudo time stepping and are not steps anyone
  // asked to be post-processed; only the state converged to is
  TimeStep tp;
  Stat sim;
  sim.unaccelerated();
  sim.tolerance(1.0e-8);
  REQUIRE(sim.solveStep(tp));

  REQUIRE(sim.nSteps > 1);
  REQUIRE(sim.nPost == 1);
}


TEST_CASE("TestSIMStat.TheAccelerationEarnsItsKeep")
{
  Quiet quiet;
  // The point of the acceleration: the same fixed point, in fewer steps
  TimeStep tp;
  Stat plain;
  plain.unaccelerated();
  plain.tolerance(1.0e-8);
  plain.measureOn("slow");
  REQUIRE(plain.solveStep(tp));

  Stat accelerated;
  accelerated.tolerance(1.0e-8);
  accelerated.measureOn("slow");
  REQUIRE(accelerated.solveStep(tp));

  INFO("the accelerated iteration is not shorter than the plain one");
  REQUIRE(accelerated.nSteps < plain.nSteps);

  for (int d = 1; d <= StepSIM::nDofs; d++)
    REQUIRE_THAT(accelerated.SIMsolution::getSolution()(d), WithinAbs(1.0, 1.0e-6));
}


TEST_CASE("TestSIMStat.StopsAtTheIterationCap")
{
  Quiet quiet;
  // A tolerance it cannot reach must end the loop rather than run forever
  TimeStep tp;
  Stat sim;
  sim.unaccelerated();
  sim.tolerance(1.0e-30);
  sim.cap(7);
  REQUIRE(sim.solveStep(tp));

  REQUIRE(sim.nSteps == 8); // the cap is on the counter, which starts at zero
}
