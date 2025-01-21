// $Id$
//==============================================================================
//!
//! \file TestNewmark.C
//!
//! \date May 14 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Tests the Newmark time integrator.
//!
//==============================================================================

#include "SAM.h"
#include "IntegrandBase.h"
#include "SIMgeneric.h"
#include "SIMdummy.h"

#include "GenAlphaSIM.h"
#include "NewmarkNLSIM.h"
#include "HHTSIM.h"
#include "HHTMats.h"
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


// SAM class representing a 2-DOF system.
class SAM2DOF : public SAM
{
public:
  SAM2DOF()
  {
    ndof   = neq = 2;
    nmmnpc = nel = nnod = 1;
    mmnpc  = new int[1]; mmnpc[0] = 1;
    mpmnpc = new int[2]; std::iota(mpmnpc,mpmnpc+2,1);
    madof  = new int[2]; madof[0] = 1; madof[1] = 3;
    msc    = new int[2]; msc[0] = msc[1] = 1;
    EXPECT_TRUE(this->initSystemEquations());
  }
  virtual ~SAM2DOF() {}
};


// SAM class representing a 2-DOF system where the 2nd DOF is prescribed.
class SAM2DOFprescr : public SAM
{
public:
  SAM2DOFprescr()
  {
    nmmnpc = nnod = ndof = 2;
    nmmceq = nceq = nel = neq = 1;
    mmnpc  = new int[2]; std::iota(mmnpc,mmnpc+2,1);
    mpmnpc = new int[2]; mpmnpc[0] = 1; mpmnpc[1] = 3;
    madof  = new int[3]; std::iota(madof,madof+3,1);
    msc    = new int[2]; msc[0] = 1; msc[1] = 0;
    mmceq  = new int[1]; mmceq[0] = 2;
    mpmceq = new int[2]; std::iota(mpmceq,mpmceq+2,1);
    ttcc   = new double[1]; ttcc[0] = 0.0;
    EXPECT_TRUE(this->initSystemEquations());
  }
  virtual ~SAM2DOFprescr() {}

  bool updateContrainEqs(double time, const Vector* prevSol)
  {
    if (!prevSol)
      ttcc[0] = 0.0;
    else if (prevSol->size() >= 2)
      ttcc[0] = sin(10.0*time) - (*prevSol)(2);
    else
      ttcc[0] = sin(10.0*time);
    return true;
  }
};


// Dummy integrand class, for time integration scheme testing.
class Problem : public IntegrandBase
{
public:
  Problem() : IntegrandBase(1) { memset(prm,0,sizeof(prm)); }
  virtual ~Problem() {}

  virtual void setIntegrationPrm(unsigned short int i, double p) { prm[i] = p; }
  virtual double getIntegrationPrm(unsigned short int i) const
  { return m_mode == SIM::DYNAMIC ? prm[i] : 0.0; }
  const double* getIntPrm() { return prm; }

private:
  double prm[5];
};


// Common base class for the test simulators.
class TestSIMbase : public SIMdummy<SIMgeneric>
{
public:
  TestSIMbase() : SIMdummy<SIMgeneric>(new Problem) {}
  virtual ~TestSIMbase() {}

protected:
  bool assembleMass(double M, size_t ndof = 1)
  {
    ElmMats elm;
    elm.resize(1,0);
    elm.redim(ndof);
    elm.A.front().diag(M);
    return myEqSys->assemble(&elm,1);
  }
};


// Simulator class for a single-DOF oscillator.
class SIM1DOF : public TestSIMbase
{
  const double M = 10.0;   // Mass of the oscillator
  const double K = 1000.0; // Stiffness of the oscillator
  const double F = 1.0;    // External load (constant)

public:
  SIM1DOF() { mySam = new SAM1DOF(); }
  virtual ~SIM1DOF() {}

  double getInitAcc() const { return F/M; }

  virtual bool assembleSystem(const TimeDomain& time,
                              const Vectors& prevSol,
                              bool newLHSmatrix, bool)
  {
    myEqSys->initialize(newLHSmatrix);

    bool ok;
    if (myProblem->getMode() == SIM::MASS_ONLY)
      ok = this->assembleMass(M);
    else {
      NewmarkMats* elm;
      const double* intPrm = static_cast<Problem*>(myProblem)->getIntPrm();
      bool useHHT = intPrm[4] == 1.0;
      if (useHHT)
        elm = new HHTMats(intPrm[2],intPrm[0],intPrm[1]);
      else
        elm = new NewmarkMats(intPrm[0],intPrm[1],intPrm[2],intPrm[3]);
      elm->rhsOnly = !newLHSmatrix;
      elm->resize(3, useHHT ? 2 : 1);
      elm->redim(1);
      elm->setStepSize(time.dt,time.it);
      elm->A[1].diag(M); // Mass matrix
      elm->A[2].diag(K); // Stiffness matrix
      elm->b[0] = -K*prevSol.front(); // Elastic forces
      elm->vec = prevSol;
      ok = myEqSys->assemble(elm,1);
    }

    // Add in the external load
    ok &= mySam->assembleSystem(*myEqSys->getVector(),&F,1);

    return ok && myEqSys->finalize(newLHSmatrix);
  }
};


// Simulator class for a 2-DOF damped oscillator.
class SIM2DOFdmp : public TestSIMbase
{
  const double M  = 10.0;  // Mass of the oscillator
  const double Kx = 200.0; // Stiffness of the oscillator in x-direction
  const double Ky = 100.0; // Stiffness of the oscillator in y-direction
  const double Cx = 100.0; // Velocity-proportional damping coefficient
  const double Fx = 100.0; // External load in x-direction (constant)
  const double Fy = 1.0;   // External load in y-direction (constant)

public:
  SIM2DOFdmp() { mySam = new SAM2DOF(); }
  virtual ~SIM2DOFdmp() {}

  double getInitAcc(size_t i = 0) const { return (i == 0 ? Fx : Fy)/M; }

  virtual bool assembleSystem(const TimeDomain& time,
                              const Vectors& prevSol,
                              bool newLHSmatrix, bool)
  {
    myEqSys->initialize(newLHSmatrix);

    bool ok;
    if (myProblem->getMode() == SIM::MASS_ONLY)
      ok = this->assembleMass(M,2);
    else {
      NewmarkMats* elm;
      double v = prevSol[prevSol.size()-2][0];
      const double* intPrm = static_cast<Problem*>(myProblem)->getIntPrm();
      bool nlDyn = intPrm[3] <= 0.0;
      if (nlDyn)
        elm = new HHTMats(intPrm[2],intPrm[0],intPrm[1],intPrm[4]!=1.0);
      else
        elm = new NewmarkMats(intPrm[0],intPrm[1],intPrm[2],intPrm[3]);
      elm->resize(nlDyn ? 5 : 4, nlDyn ? 4 : 2);
      elm->redim(2);
      elm->setStepSize(time.dt,time.it);
      elm->A[1].diag(M); // Mass matrix
      elm->A[2](1,1) = Kx; // Stiffness matrix
      elm->A[2](2,2) = Ky; // Stiffness matrix
      if (nlDyn) elm->A[3].clear(); // No geometric stiffness
      elm->A[nlDyn ? 4 : 3](1,1) = Cx*v; // Damping matrix
      elm->b[0](1) = -Kx*prevSol[0][0]; // Elastic force in X-direction
      elm->b[0](2) = -Ky*prevSol[0][1]; // Elastic force in Y-direction
      elm->b[nlDyn ? 3 : 1](1) = 0.5*Cx*v*v; // Damping force, integral of Cx*v
      elm->vec = prevSol;
      ok = myEqSys->assemble(elm,1);
      delete elm;
    }

    // Add in the external load
    double F[2] = { Fx, Fy };
    ok &= mySam->assembleSystem(*myEqSys->getVector(),F,1);

    return ok && myEqSys->finalize(newLHSmatrix);
  }
};


// Simulator class for a two-DOF oscillator with prescribed motion.
class SIM2DOFprescr : public TestSIMbase
{
  const double M = 1.0;    // Mass of the oscillator
  const double K = 1000.0; // Stiffness of the oscillator

public:
  SIM2DOFprescr() { mySam = new SAM2DOFprescr(); }
  virtual ~SIM2DOFprescr() {}

  virtual bool assembleSystem(const TimeDomain& time,
                              const Vectors& prevSol,
                              bool newLHSmatrix, bool)
  {
    myEqSys->initialize(newLHSmatrix);

    const double* intPrm = static_cast<Problem*>(myProblem)->getIntPrm();
    NewmarkMats elm(intPrm[0],intPrm[1],intPrm[2],intPrm[3]);
    elm.resize(3,1); elm.redim(2);
    elm.setStepSize(time.dt,time.it);
    elm.A[1].diag(M); // Mass matrix
    elm.A[2].diag(K); // Stiffness matrix
    elm.A[2](2,1) = elm.A[2](1,2) = -K;
    elm.b[0] = elm.A[2]*prevSol[0]; // Elastic forces
    elm.b[0] *= -1.0;
    elm.vec = prevSol;

    return myEqSys->assemble(&elm,1) && myEqSys->finalize(newLHSmatrix);
  }

  virtual bool updateDirichlet (double time, const Vector* prevSol)
  {
    return static_cast<SAM2DOFprescr*>(mySam)->updateContrainEqs(time,prevSol);
  }
};


// Newmark time integrator with numerical damping (alpha_H = -0.1).
class Newmark : public NewmarkSIM
{
public:
  Newmark(SIMbase& sim, bool useDispl, int nupd = 20) : NewmarkSIM(sim)
  {
    beta = 0.3025; gamma = 0.6;
    solveDisp = useDispl;
    predictor = useDispl ? 'd' : 'a';
    nupdat = nupd;
    this->initPrm();
    this->initSol(3);
  }
  virtual ~Newmark() {}
};


// Generalized-alpha time integrator with numerical damping (alpha_H = -0.1).
class GenAlpha : public GenAlphaSIM
{
public:
  GenAlpha(SIMbase& sim, bool useDispl) : GenAlphaSIM(sim)
  {
    solveDisp = useDispl;
    predictor = useDispl ? 'd' : 'a';
    this->initPrm();
    this->initSol(3);
  }
  virtual ~GenAlpha() {}
};


void runSingleDof (NewmarkSIM& solver, double rtol = 0.5e-11)
{
  TimeStep tp;
  tp.time.dt = 0.01;
  tp.stopTime = 0.65;
  solver.opt.solver = LinAlg::DENSE;

  ASSERT_TRUE(solver.initEqSystem());
  ASSERT_TRUE(solver.initAcc());
  EXPECT_FLOAT_EQ(solver.getAcceleration()[0],0.1);

  //               at t=0.1           at t=0.25         at t=0.5
  double u[3] = { 0.000457484252515, 0.00178698471292, 0.000732016593476 };
  double v[3] = { 0.008368445734720, 0.00592764975245,-0.009365075630580 };
  double a[3] = { 0.054251574748500,-0.07869847129160, 0.026798340652400 };

  while (solver.advanceStep(tp))
  {
    ASSERT_TRUE(solver.solveStep(tp) == SIM::CONVERGED);
    double dis = solver.getSolution()[0];
    double vel = solver.getVelocity()[0];
    double acc = solver.getAcceleration()[0];

    // Check the response at three randomly selected time steps
    if (tp.step == 10) {
      EXPECT_NEAR(dis,u[0], (dis+u[0])*rtol);
      EXPECT_NEAR(vel,v[0], (vel+v[0])*rtol);
      EXPECT_NEAR(acc,a[0], (acc+a[0])*rtol);
    }
    else if (tp.step == 25) {
      EXPECT_NEAR(dis,u[1], (dis+u[1])*rtol);
      EXPECT_NEAR(vel,v[1], (vel+v[1])*rtol);
      EXPECT_NEAR(acc,a[1],-(acc+a[1])*rtol);
    }
    else if (tp.step == 50) {
      EXPECT_NEAR(dis,u[2], (dis+u[2])*rtol);
      EXPECT_NEAR(vel,v[2],-(vel+v[2])*rtol);
      EXPECT_NEAR(acc,a[2], (acc+a[2])*rtol);
    }
  }
}


static void printVec (const char* name, const Vector& vec)
{
  std::ios::fmtflags stdFlags = std::cout.flags(std::ios::scientific);
  std::streamsize stdPrec = std::cout.precision(6);
  std::cout <<"  "<< name <<" =";
  for (double v : vec) std::cout <<" "<< v;
  std::cout << std::endl;
  std::cout.flags(stdFlags);
  std::cout.precision(stdPrec);
}


void runTwoDof (NewmarkSIM& solver, double refAcc, double rtol = 0.5e-11)
{
  TimeStep tp;
  tp.time.dt = 0.01;
  tp.stopTime = 0.65;
  solver.opt.solver = LinAlg::DENSE;
  ASSERT_TRUE(solver.initEqSystem());
  ASSERT_TRUE(solver.initAcc());
  EXPECT_FLOAT_EQ(solver.getAcceleration()[0],refAcc);

  while (solver.advanceStep(tp))
  {
    ASSERT_TRUE(solver.solveStep(tp) == SIM::CONVERGED);
    printVec("u",solver.getSolution());
    printVec("v",solver.getVelocity());
    printVec("a",solver.getAcceleration());
  }
}


void runPrescribed (NewmarkSIM& solver, double rtol = 0.5e-11)
{
  TimeStep tp;
  tp.time.dt = 0.01;
  tp.stopTime = 0.65;
  solver.opt.solver = LinAlg::DENSE;
  ASSERT_TRUE(solver.initEqSystem());

  //              at t=0.1         at t=0.25         at t=0.5
  double u[3] = {0.9312639435267, 0.3547915343361, -1.075289543029 };
  double v[3] = { 16.61378563136, -8.930847244282,  11.85822968077 };
  double a[3] = {-89.79295871881,  243.6806097678,  116.3652683660 };

  while (solver.advanceStep(tp))
  {
    ASSERT_TRUE(solver.solveStep(tp) == SIM::CONVERGED);
    double dis = solver.getSolution()[0];
    double vel = solver.getVelocity()[0];
    double acc = solver.getAcceleration()[0];

    // Check the response at three randomly selected time steps
    if (tp.step == 10) {
      EXPECT_NEAR(dis,u[0], (dis+u[0])*rtol);
      EXPECT_NEAR(vel,v[0], (vel+v[0])*rtol);
      EXPECT_NEAR(acc,a[0],-(acc+a[0])*rtol);
    }
    else if (tp.step == 25) {
      EXPECT_NEAR(dis,u[1], (dis+u[1])*rtol);
      EXPECT_NEAR(vel,v[1],-(vel+v[1])*rtol);
      EXPECT_NEAR(acc,a[1], (acc+a[1])*rtol);
    }
    else if (tp.step == 50) {
      EXPECT_NEAR(dis,u[2],-(dis+u[2])*rtol);
      EXPECT_NEAR(vel,v[2], (vel+v[2])*rtol);
      EXPECT_NEAR(acc,a[2], (acc+a[2])*rtol);
    }
  }
}


TEST(TestNewmark, SingleDOFa)
{
  SIM1DOF simulator;
  Newmark integrator(simulator,false);
  runSingleDof(integrator);
}

TEST(TestGenAlpha, SingleDOFa)
{
  SIM1DOF simulator;
  GenAlpha integrator(simulator,false);
  runSingleDof(integrator,0.02);
}

TEST(TestNewmark, SingleDOFu)
{
  SIM1DOF simulator;
  Newmark integrator(simulator,true,0);
  runSingleDof(integrator);
}

TEST(TestHHT, SingleDOFu)
{
  SIM1DOF simulator;
  HHTSIM integrator(simulator);
  integrator.initPrm();
  integrator.initSol();
  runSingleDof(integrator,0.9);
}

TEST(TestNewmark, Damped)
{
  SIM2DOFdmp simulator;
  Newmark integrator(simulator,true);
  runTwoDof(integrator,simulator.getInitAcc());
}

TEST(TestNewmarkNL, Damped)
{
  SIM2DOFdmp simulator;
  NewmarkNLSIM integrator(simulator);
  integrator.initPrm();
  integrator.initSol();
  runTwoDof(integrator,simulator.getInitAcc());
}

TEST(TestHHT, Damped)
{
  SIM2DOFdmp simulator;
  HHTSIM integrator(simulator);
  integrator.initPrm();
  integrator.initSol();
  runTwoDof(integrator,simulator.getInitAcc());
}

TEST(TestNewmark, Prescribed)
{
  SIM2DOFprescr simulator;
  Newmark integrator(simulator,true);
  runPrescribed(integrator);
}

/* does not work, yet
TEST(TestGenAlpha, SingleDOFu)
{
  SIM1DOF simulator;
  GenAlpha integrator(simulator,true);
  runSingleDof(integrator,0.02);
}

TEST(TestGenAlpha, Prescribed)
{
  SIM2DOFprescr simulator;
  GenAlpha integrator(simulator,true);
  runPrescribed(integrator,0.02);
}
*/
