//==============================================================================
//!
//! \file TestSIMargsBase.C
//!
//! \date Oct 15 2026
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Tests for base class for command-line argument parsing.
//!
//==============================================================================

#include "SIMargsBase.h"

#include "Catch2Support.h"
#include <cstring>


TEST_CASE("TestSIMargsBase.preParse")
{
  constexpr int argc = 7;
  const char* testarg[argc] = {
    "myProgram",
    "-2D",
    "-adap",
    "-foo","22",
    "-bar",
    "with_include.xml"
  };

  char* argv[argc];
  for (int i = 0; i < argc; i++)
    argv[i] = const_cast<char*>(testarg[i]);

  int foo = 0;
  bool bar = false;

  auto external_args = [argv,&foo,&bar](int& i)
  {
    if (!strcmp(argv[i],"-foo") && i+1 < argc)
      foo = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-bar"))
      bar = true;
    else
      return false;
    return true;
  };

  SIMargsBase args("Test");
  char* infile = args.preParse(argc,argv,external_args);
  REQUIRE(infile != nullptr);
  CHECK(!strcmp(infile,"with_include.xml"));
  CHECK(args.dim == 2);
  CHECK(args.adap);
  CHECK(bar);
}
