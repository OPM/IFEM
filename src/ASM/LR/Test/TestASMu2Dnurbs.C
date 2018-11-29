//==============================================================================
//!
//! \file TestASMu2Dnurbs.C
//!
//! \date Feb 8 2019
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for driver for assembly of unstructured 2D NURBS FE models.
//!
//==============================================================================

#include "ASMs2D.h"
#include "ASMu2Dnurbs.h"
#include "Vec3.h"

#include <LRSpline/LRSplineSurface.h>
#include "gtest/gtest.h"
#include <fstream>


TEST(TestASMu2Dnurbs, EvalPoint)
{
  ASMu2Dnurbs u2Dpch(2, 1);
  std::ifstream g2file("src/ASM/LR/Test/refdata/hole2D.g2");
  std::ifstream g2file2("src/ASM/LR/Test/refdata/hole2D.g2");
  ASMs2D s2Dpch(2,1);
  ASSERT_TRUE(u2Dpch.read(g2file));
  ASSERT_TRUE(s2Dpch.read(g2file2));
  ASSERT_TRUE(u2Dpch.generateFEMTopology());
  ASSERT_TRUE(s2Dpch.generateFEMTopology());

  double xi[2] = {0.1, 0.1};
  double param1[2], param2[2];
  Vec3 x1, x2;
  s2Dpch.evalPoint(xi,param1,x1);
  u2Dpch.evalPoint(xi,param2,x2);
  EXPECT_FLOAT_EQ(param1[0], param2[0]);
  EXPECT_FLOAT_EQ(param1[1], param2[1]);
  EXPECT_FLOAT_EQ(x1[0], x2[0]);
  EXPECT_FLOAT_EQ(x1[1], x2[1]);
}
