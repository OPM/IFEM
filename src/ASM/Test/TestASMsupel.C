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

#include "Catch2Support.h"


namespace {

struct TestCase
{
  const char* file;
  size_t      nsup;
};

}


TEST_CASE("TestASMsup.Read")
{
  const TestCase param = GENERATE(
    TestCase{"src/ASM/Test/refdata/Supel.dat", 2U},
    TestCase{"src/ASM/Test/refdata/kjoint.dat", 4U}
  );

  SECTION(param.file) {
    ASMsupel pch;
    ASMbase::resetNumbering();
    std::cout <<"Checking "<< param.file << std::endl;
    std::ifstream is(param.file);
    REQUIRE(pch.read(is));
    REQUIRE(!pch.empty());
    REQUIRE(pch.generateFEMTopology());
    REQUIRE(pch.getNoNodes() == param.nsup);
  }
}
