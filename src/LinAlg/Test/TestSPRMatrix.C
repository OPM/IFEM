// $Id$
//==============================================================================
//!
//! \file TestSPRMatrix.C
//!
//! \date Sep 12 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Unit tests for the SPR equation solver.
//!
//==============================================================================

#include "SAM.h"
#include "SIMgeneric.h"
#include "SIMdummy.h"

#include "SPRMatrix.h"
#include "ElmMats.h"
#include "AlgEqSystem.h"

#include "gtest/gtest.h"

namespace // To avoid conflicting names with other tests
{

// SAM class representing a single-DOF system.
class SAM1DOF : public SAM
{
public:
  SAM1DOF()
  {
    nmmnpc = nel = nnod = ndof = neq = 1;
    mmnpc  = new int[1]; mmnpc[0] = 1;
    mpmnpc = new int[2]; mpmnpc[0] = 1; mpmnpc[1] = 2;
    madof  = new int[2]; madof[0] = 1; madof[1] = 2;
    msc    = new int[1]; msc[0] = 1;
    mpmceq = new int[1]; mpmceq[1] = 0;
    EXPECT_TRUE(this->initSystemEquations());
  }
  virtual ~SAM1DOF() {}
};


// SAM class representing a 2-DOF system consisting of 2 elements and 3 nodes.
class SAM2DOF : public SAM
{
public:
  SAM2DOF()
  {
    ndof   = 6;
    neq    = 2;
    nmmnpc = 4;
    nel    = 2;
    nnod   = 3;
    mmnpc  = new int[4]; mmnpc[0] = 1; mmnpc[1] = 2; mmnpc[2] = 2; mmnpc[3] = 3;
    mpmnpc = new int[3]; mpmnpc[0] = 1; mpmnpc[1] = 3; mpmnpc[2] = 5;
    madof  = new int[4]; madof[0] = 1; madof[1] = 3; madof[2] = 5; madof[3] = 7;
    msc    = new int[6]; msc[0]=msc[1]=msc[4]=msc[5] = 0; msc[2]=msc[3] = 1;
    mpmceq = new int[1]; mpmceq[1] = 0;
    EXPECT_TRUE(this->initSystemEquations());
  }
  virtual ~SAM2DOF() {}
};


// Simulator class for a single-DOF skew bar.
class Bar1DOF : public SIMdummy<SIMgeneric>
{
public:
  Bar1DOF() { mySam = new SAM1DOF(); }
  virtual ~Bar1DOF() {}

  using SIMdummy<SIMgeneric>::assembleSystem;
  virtual bool assembleSystem(const TimeDomain&, const Vectors&, bool, bool)
  {
    const double x = 5.0;             // X-coordinate of end point
    const double y = 2.0;             // Y-coordinate of end point
    const double K = 26925.824035673; // Axial stiffness
    const double F = -200.0;          // External load (constant)

    double L  = hypot(x,y); // Length of the bar
    double Km = K/L;        // Nominal material stiffness
    Vec3   X(x/L,y/L,0.0);  // Beam axis

    myEqSys->initialize(true);
    ElmMats elm;
    elm.resize(1,1);
    elm.redim(1);

    elm.A.front().fill(Km*X.y*X.y); // Stiffness matrix
    if (!myEqSys->assemble(&elm,1))
      return false;

    // Add in the external load (local DOF 1 in node 1)
    if (!mySam->assembleSystem(*myEqSys->getVector(),F,{1,1}))
      return false;

    return myEqSys->finalize(true);
  }
};


// Simulator class for the two-DOF bar system.
class Bar2DOF : public SIMdummy<SIMgeneric>
{
public:
  Bar2DOF() { mySam = new SAM2DOF(); }
  virtual ~Bar2DOF() {}

  using SIMdummy<SIMgeneric>::assembleSystem;
  virtual bool assembleSystem(const TimeDomain&, const Vectors&, bool, bool)
  {
    const double x = 5.0;             // X-coordinate of top point
    const double y = 2.0;             // Y-coordinate of top point
    const double z = 8.0;             // X-coordinate of end point
    const double K = 26925.824035673; // Axial stiffness
    const double F = -200.0;          // External load (constant)

    double x2 = z-x;
    double L1 = hypot(x,y);  // Length of first bar
    double L2 = hypot(x2,y); // Length of second bar
    double K1 = K/L1;        // Nominal material stiffness, element 1
    double K2 = K/L2;        // Nominal material stiffness, element 2

    Vec3 X1( x/L1, y/L1,0.0); // Beam axis, element 1
    Vec3 X2(x2/L2,-y/L2,0.0); // Beam axis, element 2

    // Lambda function setting up the 4x4 element stiffness matrix
    auto&& barElement = [](Matrix& Km, const Vec3& Xaxis, double Knom)
    {
      for (int j = 1; j <= 2; j++)
        for (int i = 1; i <= 2; i++)
        {
          Km(  i,  j) = Knom*Xaxis[i-1]*Xaxis[j-1];
          Km(2+i,2+j) = Km(i  ,j);
          Km(  i,2+j) = Km(2+i,j) = -Km(i,j);
        }
    };

    // Assemble the stiffness matrix
    myEqSys->initialize(true);
    ElmMats elm;
    elm.resize(1,1);
    elm.redim(4);

    barElement(elm.A.front(),X1,K1);
    if (!myEqSys->assemble(&elm,1))
      return false;
    barElement(elm.A.front(),X2,K2);
    if (!myEqSys->assemble(&elm,2))
      return false;

    // Add in the external load (local DOF 2 in node 2)
    if (!mySam->assembleSystem(*myEqSys->getVector(),F,{2,2}))
      return false;

    return myEqSys->finalize(true);
  }
};
}


TEST(TestSPRMatrix, SingleDOF)
{
  Bar1DOF simulator;
  Vectors solution(1);
#ifdef HAS_SPR
  ASSERT_TRUE(simulator.initSystem(LinAlg::SPR));
#else
  ASSERT_TRUE(simulator.initSystem(LinAlg::DENSE));
#endif
  ASSERT_TRUE(simulator.assembleSystem());
  ASSERT_TRUE(simulator.solveSystem(solution,1));
  std::cout <<"Solution vector:" << solution.front();
}


TEST(TestSPRMatrix, TwoDOF)
{
  Bar2DOF simulator;
  Vectors solution(1);
#ifdef HAS_SPR
  ASSERT_TRUE(simulator.initSystem(LinAlg::SPR));
#else
  ASSERT_TRUE(simulator.initSystem(LinAlg::DENSE));
#endif
  ASSERT_TRUE(simulator.assembleSystem());
  ASSERT_TRUE(simulator.solveSystem(solution,1));
  std::cout <<"Solution vector:" << solution.front();
}
