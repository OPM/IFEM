//==============================================================================
//!
//! \file TestInitialConditions.C
//!
//! \date Apr 6 2016
//!
//! \author Timo van Opstal / NTNU
//!
//! \brief Tests initial condition handling.
//!
//==============================================================================

#include "SIM2D.h"
#include "Function.h"

#include "gtest/gtest.h"


class TestSIM : public SIM2D
{
public:
  TestSIM() : SIM2D(4) {}
  virtual ~TestSIM() {}
  const InitialCondVec& getIC() const { return myICs.begin()->second; }
};


TEST(TestInitialConditions, Parse)
{
  const char* geo = "<geometry><topologysets>"
    "<set name='left' type='edge'><item patch='1'>1</item></set>"
    "<set name='right' type='edge'><item patch='1'>2</item></set>"
    "<set name='bottom' type='edge'><item patch='1'>3</item></set>"
    "</topologysets></geometry>";

  const char* dbc = "<boundaryconditions>"
    "<dirichlet set='left' comp='1'>1</dirichlet>"
    "<dirichlet set='right' component='2'>2</dirichlet>"
    "<dirichlet set='bottom' comp='3' component='-3'>3</dirichlet>"
    "<dirichlet set='bottom' component='-4' comp='4'>4</dirichlet>"
    "</boundaryconditions>";

  TestSIM sim;
  ASSERT_TRUE(sim.loadXML(geo));
  ASSERT_TRUE(sim.loadXML(dbc));
  ASSERT_TRUE(sim.loadXML("<initialcondition field='solution' type='expression' comp='1'>1</initialcondition>"));
  ASSERT_TRUE(sim.loadXML("<initialcondition field='solution' type='expression' component='2'>2</initialcondition>"));
  ASSERT_TRUE(sim.loadXML("<initialcondition field='solution' type='expression' comp='3' component='-3'>3</initialcondition>"));
  ASSERT_TRUE(sim.loadXML("<initialcondition field='solution' type='expression' component='-4' comp='4'>4</initialcondition>"));

  // Recognize both comp and component attributes and correct priority
  // Boundary conditions
  for (int i = 1; i < 5; i++)
    EXPECT_FLOAT_EQ((double)i,(*sim.getSclFunc(i))(Vec3()));

  // Initial conditions
  for (const SIMinput::ICInfo& info : sim.getIC())
    switch (info.component) {
      case 1:
        EXPECT_EQ(info.function,"1");
        break;
      case 2:
        EXPECT_EQ(info.function,"2");
        break;
      case 3:
        EXPECT_EQ(info.function,"3");
        break;
      case 4:
        EXPECT_EQ(info.function,"4");
        break;
      default:
        EXPECT_TRUE(false);
      }
}
