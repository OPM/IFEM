//==============================================================================
//!
//! \file TestASMsupel.C
//!
//! \date Mar 30 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Tests for superelement patches.
//!
//==============================================================================

#include "ASMsupel.h"
#include <fstream>

#include "gtest/gtest.h"


TEST(TestASMsupel, Read)
{
  ASMsupel pch;
  ASMbase::resetNumbering();
  std::ifstream is("src/ASM/Test/refdata/Supel.dat");
  ASSERT_TRUE(pch.read(is));
  ASSERT_FALSE(pch.empty());
  ASSERT_TRUE(pch.generateFEMTopology());
  EXPECT_EQ(pch.getNoNodes(),2U);
}
