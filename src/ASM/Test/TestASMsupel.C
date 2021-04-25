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


struct TestCase
{
  const char* file;
  size_t      nsup;
};


class TestASMsup : public testing::Test,
                   public testing::WithParamInterface<TestCase>
{
};


TEST_P(TestASMsup, Read)
{
  ASMsupel pch;
  ASMbase::resetNumbering();
  std::cout <<"Checking "<< GetParam().file << std::endl;
  std::ifstream is(GetParam().file);
  ASSERT_TRUE(pch.read(is));
  ASSERT_FALSE(pch.empty());
  ASSERT_TRUE(pch.generateFEMTopology());
  EXPECT_EQ(pch.getNoNodes(),GetParam().nsup);
}


const std::vector<TestCase> testFiles = {
  { "src/ASM/Test/refdata/Supel.dat", 2U },
  { "src/ASM/Test/refdata/kjoint.dat", 4U }};


INSTANTIATE_TEST_CASE_P(TestASMsup, TestASMsup, testing::ValuesIn(testFiles));
