//==============================================================================
//!
//! \file TestElementBlock.C
//!
//! \date Oct 7 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for storage of a standard FE grid block of uniform element type.
//!
//==============================================================================

#include "ElementBlock.h"

#include "Catch2Support.h"


TEST_CASE("TestElementBlock.Resize")
{
  ElementBlock block(8);
  block.resize(3, 3, 3);
  REQUIRE(block.getNoNodes() == 27);

  block.unStructResize(3, 8);
  REQUIRE(block.getNoNodes() == 8);
}
