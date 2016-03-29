//==============================================================================
//!
//! \file TestScopedLogger.C
//!
//! \date Oct 10 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for scoped logger.
//!
//==============================================================================

#include "ScopedLogger.h"

#include "gtest/gtest.h"

void MockFunction(std::ostream& stream)
{
  ScopedLogger log("MockFunction", stream);
}

TEST(TestScopedLogger, General)
{
  std::stringstream str;
  MockFunction(str);

  char tmp[1024];
  str.getline(tmp, 1024);
#ifdef HAVE_MPI
  ASSERT_STREQ(tmp, "[0]: Entering \"MockFunction\"");
#else
  ASSERT_STREQ(tmp, "Entering \"MockFunction\"");
#endif
  str.getline(tmp, 1024);
#ifdef HAVE_MPI
  ASSERT_STREQ(tmp, "[0]: Exiting \"MockFunction\"");
#else
  ASSERT_STREQ(tmp, "Exiting \"MockFunction\"");
#endif


}
