// $Id$
//==============================================================================
//!
//! \file TestNonLinSIM.C
//!
//! \date Nov 25 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Tests the nonlinear quasi-static solution driver.
//!
//==============================================================================

#include "SAM.h"
#include "SIMgeneric.h"
#include "SIMdummy.h"

#include "NonLinSIM.h"
#include "ElmMats.h"
#include "AlgEqSystem.h"
#include "TimeStep.h"

#include "gtest/gtest.h"
#include <numeric>


// SAM class representing a single-DOF system.
class SAM1DOF : public SAM
{
public:
  SAM1DOF()
  {
    nmmnpc = nel = nnod = ndof = neq = 1;
    mmnpc  = new int[1]; mmnpc[0] = 1;
    mpmnpc = new int[2]; std::iota(mpmnpc,mpmnpc+2,1);
    madof  = new int[2]; std::iota(madof,madof+2,1);
    msc    = new int[1]; msc[0] = 1;
    EXPECT_TRUE(this->initSystemEquations());
  }
  virtual ~SAM1DOF() {}
};


// Simulator class for a single-DOF skew bar.
class Bar1DOF : public SIMdummy<SIMgeneric>
{
public:
  Bar1DOF() { mySam = new SAM1DOF(); }
  virtual ~Bar1DOF() {}
  virtual bool assembleSystem(const TimeDomain&, const Vectors& prevSol,
                              bool newLHSmatrix, bool)
  {
    const double x = 5.0;             // Initial X-coordinate of end point
    const double y = 2.0;             // Initial Y-coordinate of end point
    const double K = 26925.824035673; // Axial stiffness
    const double F = -200.0;          // External load (constant)

    double u   = prevSol.front()[0]; // Current deflection
    double L   = hypot(x,y+u); // Current length of the bar
    double L0  = hypot(x,y);   // Initial length of the bar
    double eps = L/L0 - 1.0;   // Axial strain
    double N   = K*eps;        // Axial force
    double Km  = K/L0;         // Nominal material stiffness
    double Kg  = N/L;          // Nominal geometric stiffness

    Vec3 X(x/L,(y+u)/L,0.0); // Beam axis
    Vec3 Y(-X.y,X.x,0.0);    // Normal axis

    ElmMats elm;
    elm.resize(1,1);
    elm.redim(1);
    elm.vec.resize(1);
    elm.A.front().fill(Km*X.y*X.y + Kg*Y.y*Y.y); // Stiffness matrix
    elm.b.front().fill(-N*X.y);                  // Internal forces
    elm.vec.front() = prevSol.front();
    myEqSys->initialize(newLHSmatrix);
    if (!myEqSys->assemble(&elm,1))
      return false;

    // Add in the external load
    if (!mySam->assembleSystem(*myEqSys->getVector(),&F,1))
      return false;

    return myEqSys->finalize(newLHSmatrix);
  }
};


// Nonlinear simulation driver with line search.
class TestNonLinSIM : public NonLinSIM
{
public:
  TestNonLinSIM(SIMbase& sim, double lsTol = 0.0) : NonLinSIM(sim)
  {
    eta = lsTol;
    rTol = 1.0e-16;
  }
  virtual ~TestNonLinSIM() {}
};


static void runSingleDof (NonLinSIM& solver, int& n, double& s)
{
  TimeStep tp;
  ASSERT_TRUE(solver.initSol());
  ASSERT_TRUE(solver.advanceStep(tp));
  ASSERT_EQ(solver.solveStep(tp),SIM::CONVERGED);
  s = solver.getSolution()[0]; n = tp.iter;
  std::cout <<"  Solution "<< s <<" obtained in "<< n <<" iterations.\n";
}


TEST(TestNonLinSIM, SingleDOF)
{
  Bar1DOF simulator;
  ASSERT_TRUE(simulator.initSystem(LinAlg::DENSE));

  TestNonLinSIM integrator1(simulator);
  TestNonLinSIM integrator2(simulator,0.0001);

  int    n1, n2;
  double s1, s2;
  runSingleDof(integrator1,n1,s1);
  runSingleDof(integrator2,n2,s2);

  EXPECT_FLOAT_EQ(s1,s2);
  EXPECT_EQ(n1,5);
  EXPECT_EQ(n2,3);
}
