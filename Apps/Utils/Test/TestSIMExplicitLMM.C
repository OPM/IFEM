// $Id$
//==============================================================================
//!
//! \file TestSIMExplicitLMM.C
//!
//! \date Nov 24 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for the explicit linear multistep time stepping of SIMs.
//!
//==============================================================================

#include "SIMExplicitLMM.h"
#include "ProcessAdm.h"
#include "TimeDomain.h"

#include "Catch2Support.h"

#include <cstring>
#include <map>
#include <string>


namespace {

/*!
  \brief A solver of du/dt = u, whose next solve can be made to fail after
  the right-hand side has been assembled.
*/

class ODESolver
{
public:
  static int msgLevel; //!< Output level, set by the time stepping

  const ProcessAdm& getProcessAdm() const { return adm; }
  bool initDirichlet(double) { return true; }
  bool updateDirichlet(double, const Vector*, bool = false) { return true; }
  bool applyDirichlet(Vector&) const { return true; }
  bool hasIC(const std::string&) const { return false; }
  void setMode(SIM::SolutionMode) {}
  void printSolutionSummary(const Vector&, int, const char*) {}

  //! \brief The problem, which only names the solution.
  struct Problem
  {
    //! \brief Returns the name of the solution field.
    std::string getField1Name(size_t) const { return "u"; }
  };
  //! \brief Returns the problem.
  const Problem* getProblem() const { return &problem; }

  Vector& getSolution(int i = 0) { return sol[i]; }

  bool assembleSystem(const TimeDomain&, const Vectors& u, bool = true)
  {
    rhs = StdVector(std::vector<double>{u.front()(1)});
    return true;
  }

  SystemVector* getRHSvector(size_t, bool copy)
  {
    return copy ? rhs.copy() : &rhs;
  }

  void addToRHSvector(size_t, const SystemVector& vec, double scale)
  {
    rhs.add(vec,scale);
  }

  bool solveSystem(Vector& x)
  {
    if (failNext) {
      failNext = false;
      return false;
    }
    x = rhs;
    return true;
  }

  bool advanceStep(TimeStep&)
  {
    sol[1] = sol[0];
    return true;
  }

  bool serialize(std::map<std::string,std::string>& data) const
  {
    for (int i = 0; i < 2; i++)
      data["u" + std::to_string(i)] =
        std::string(reinterpret_cast<const char*>(sol[i].ptr()),sizeof(double));
    return true;
  }

  bool deSerialize(const std::map<std::string,std::string>& data)
  {
    for (int i = 0; i < 2; i++)
      std::memcpy(sol[i].ptr(),data.at("u" + std::to_string(i)).data(),
                  sizeof(double));
    return true;
  }

  Vector sol[2] = { Vector(std::vector<double>{1.0}),
                    Vector(std::vector<double>{1.0}) }; //!< The solutions
  bool failNext = false; //!< If \e true, the next solve fails

private:
  ProcessAdm adm;  //!< Process administrator
  Problem problem; //!< The problem
  StdVector  rhs; //!< Right-hand side vector
};

int ODESolver::msgLevel = 0;


//! \brief Integrates five steps of size 0.1, failing once in step \a fail,
//! which is retried with half the step size.
//! \param[in] method The time integration method
//! \param[in] fail The failing step, if any
//! \param[out] exact The forward Euler solution of the actual step sizes
//! \return The solution, or -1 if the failed step could not be retried
double integrate(TimeIntegration::Method method, int fail, double& exact)
{
  ODESolver ode;
  TimeIntegration::SIMExplicitLMM<ODESolver> lmm(ode,method,false);

  TimeStep tp;
  exact = 1.0;
  std::map<std::string,std::string> state;
  REQUIRE(lmm.serialize(state));
  for (tp.step = 1; tp.step <= 5; tp.step++) {
    tp.time.dt = 0.1;
    tp.time.t += tp.time.dt;
    REQUIRE(lmm.advanceStep(tp));
    ode.failNext = tp.step == fail;
    if (!lmm.solveStep(tp)) {
      // Restart the step with half the step size from the state after the
      // previous one
      tp.time.t -= 0.5*tp.time.dt;
      tp.time.dt *= 0.5;
      REQUIRE(lmm.deSerialize(state));
      if (!lmm.advanceStep(tp))
        return -1.0;
      REQUIRE(lmm.solveStep(tp));
    }
    exact *= 1.0 + tp.time.dt;
    state.clear();
    REQUIRE(lmm.serialize(state));
  }

  return ode.getSolution()(1);
}

}


TEST_CASE("TestSIMExplicitLMM.Retry")
{
  // A step of forward Euler retried with a reduced step size, after its
  // failure, continues from the state after the previous step
  double exact = 0.0;
  const int fail = GENERATE(1, 3);
  const double u = integrate(TimeIntegration::AB1,fail,exact);
  REQUIRE_THAT(u, WithinRel(exact, 1.0e-14));

  // The coefficients of higher orders assume the spacing of the history,
  // which a step size reduced by a cut-back, and its recovery, does not
  // have, also when the failing step starts with a lower order
  const auto method = GENERATE(TimeIntegration::AB2, TimeIntegration::AB3);
  REQUIRE(integrate(method,fail,exact) == -1.0);
}
