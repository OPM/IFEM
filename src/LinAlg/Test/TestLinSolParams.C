//==============================================================================
//!
//! \file TestLinSolParams.C
//!
//! \date Nov 2 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for LinSolParams
//!
//==============================================================================

#include "LinSolParams.h"
#include "tinyxml.h"

#include "gtest/gtest.h"

TEST(TestLinSolParams, ParseBlocks)
{
  TiXmlDocument doc;
  doc.LoadFile("src/LinAlg/Test/refdata/linsolver_blocks.xml");

  LinSolParams params(2);
  params.read(doc.RootElement());

  // global settings
  ASSERT_EQ(params.getDimension(), 2u);
  ASSERT_EQ(params.getNoBlocks(), 2u);
  ASSERT_STREQ(params.getPreconditioner(), "fieldsplit");
  ASSERT_EQ(params.getMethod(), "bcgs");
  ASSERT_FLOAT_EQ(params.getAbsTolerance(), 1e-12);
  ASSERT_FLOAT_EQ(params.getRelTolerance(), 1e-4);
  ASSERT_FLOAT_EQ(params.getDivTolerance(), 1e6);
  ASSERT_EQ(params.getMaxIterations(), 123);

  // first block
  ASSERT_EQ(params.getBlock(0).basis, 1u);
  ASSERT_EQ(params.getBlock(0).comps, 12u);
  ASSERT_EQ(params.getBlock(0).prec, "sor");

  // second block
  ASSERT_EQ(params.getBlock(1).basis, 2u);
  ASSERT_EQ(params.getBlock(1).comps, 3u);
  ASSERT_EQ(params.getBlock(1).prec, "ilu");
  ASSERT_EQ(params.getBlock(1).ilu_fill_level, 1);
}


TEST(TestLinSolParams, ParseGAMG)
{
  TiXmlDocument doc;
  doc.LoadFile("src/LinAlg/Test/refdata/linsolver_gamg.xml");

  LinSolParams params(2);
  params.read(doc.RootElement());

  ASSERT_EQ(params.getBlock(0).gamg.type, "agg");
  ASSERT_EQ(params.getBlock(0).gamg.procEqLimit, 1000);
  ASSERT_EQ(params.getBlock(0).gamg.repartition, 1);
  ASSERT_EQ(params.getBlock(0).gamg.useAggGasm, 1);
  ASSERT_EQ(params.getBlock(0).gamg.reuseInterp, 1);
  ASSERT_FLOAT_EQ(params.getBlock(0).gamg.threshold, 0.1);
  ASSERT_EQ(params.getBlock(0).finesmoother, "compositedir");
  ASSERT_EQ(params.getBlock(0).dirSmoother.size(), 1u);
  ASSERT_EQ(params.getBlock(0).dirSmoother[0].type, "ilu");
  ASSERT_EQ(params.getBlock(0).dirSmoother[0].order, 12);
}


TEST(TestLinSolParams, ParseML)
{
  TiXmlDocument doc;
  doc.LoadFile("src/LinAlg/Test/refdata/linsolver_ml.xml");

  LinSolParams params(2);
  params.read(doc.RootElement());

  ASSERT_EQ(params.getBlock(0).ml.coarsePackage, "petsc");
  ASSERT_EQ(params.getBlock(0).ml.coarseSolver, "lu");
  ASSERT_EQ(params.getBlock(0).ml.coarsenScheme, "agg");
  ASSERT_EQ(params.getBlock(0).ml.symmetrize, 1);
  ASSERT_EQ(params.getBlock(0).ml.blockScaling, 1);
  ASSERT_EQ(params.getBlock(0).ml.putOnSingleProc, 1000);
  ASSERT_EQ(params.getBlock(0).ml.reuseInterp, 1);
  ASSERT_EQ(params.getBlock(0).ml.reusable, 1);
  ASSERT_EQ(params.getBlock(0).ml.keepAggInfo, 1);
  ASSERT_EQ(params.getBlock(0).ml.aux, 1);
  ASSERT_FLOAT_EQ(params.getBlock(0).ml.threshold, 0.1);
  ASSERT_FLOAT_EQ(params.getBlock(0).ml.auxThreshold, 0.05);
}


TEST(TestLinSolParams, ParseHypre)
{
  TiXmlDocument doc;
  doc.LoadFile("src/LinAlg/Test/refdata/linsolver_hypre.xml");

  LinSolParams params(2);
  params.read(doc.RootElement());

  ASSERT_EQ(params.getBlock(0).hypre.type, "boom");
  ASSERT_EQ(params.getBlock(0).hypre.noAggCoarse, 1);
  ASSERT_EQ(params.getBlock(0).hypre.noPathAggCoarse, 2);
  ASSERT_FLOAT_EQ(params.getBlock(0).hypre.truncation, 0.005);
  ASSERT_FLOAT_EQ(params.getBlock(0).hypre.threshold, 0.01);
  ASSERT_EQ(params.getBlock(0).hypre.coarsenScheme, "agg");
}
