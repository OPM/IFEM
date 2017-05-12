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
  TestSIM() : SIM2D(4)
  {
    EXPECT_TRUE(this->read("src/SIM/Test/refdata/input.xinp"));
  }
  virtual ~TestSIM() {}
  const InitialCondVec* getIC() const
  {
    return myICs.empty() ? nullptr : &myICs.begin()->second;
  }
};


TEST(TestInitialConditions, Parse)
{
  TestSIM sim;

  // Recognize both comp and component attributes and correct priority
  // Boundary conditions
  for (int i = 1; i < 5; i++)
    ASSERT_FLOAT_EQ((float)i,(*sim.getSclFunc(i))(Vec3()));

  // Initial conditions
  const SIMinput::InitialCondVec* ic = sim.getIC();
  ASSERT_TRUE(ic != nullptr);
  for (auto info = ic->begin(); info != ic->end(); ++info)
    switch (info->component) {
      case 1:
        ASSERT_EQ(info->function, "1");
        break;
      case 2:
        ASSERT_EQ(info->function, "2");
        break;
      case 3:
        ASSERT_EQ(info->function, "3");
        break;
      case 4:
        ASSERT_EQ(info->function, "4");
        break;
      default:
        ASSERT_TRUE(false);
      }
}
