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

#include "gtest/gtest.h"

TEST(TestElementBlock, Resize)
{
  ElementBlock block(8);
  block.resize(3, 3, 3);
  EXPECT_EQ(block.getNoNodes(), 27U);

  block.unStructResize(3, 8);
  EXPECT_EQ(block.getNoNodes(), 8U);
}
