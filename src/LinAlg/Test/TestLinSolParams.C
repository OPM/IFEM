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

  LinSolParams params;
  params.read(doc.RootElement());

  // global settings
  ASSERT_EQ(params.getNoBlocks(), 2u);
  ASSERT_EQ(params.getStringValue("pc"), "fieldsplit");
  ASSERT_EQ(params.getStringValue("type"), "bcgs");
  ASSERT_FLOAT_EQ(params.getDoubleValue("atol"), 1e-12);
  ASSERT_FLOAT_EQ(params.getDoubleValue("rtol"), 1e-4);
  ASSERT_FLOAT_EQ(params.getDoubleValue("dtol"), 1e6);
  ASSERT_EQ(params.getIntValue("maxits"), 123);

  // first block
  ASSERT_EQ(params.getBlock(0).basis, 1u);
  ASSERT_EQ(params.getBlock(0).comps, 12u);
  ASSERT_EQ(params.getBlock(0).getStringValue("pc"), "sor");

  // second block
  ASSERT_EQ(params.getBlock(1).basis, 2u);
  ASSERT_EQ(params.getBlock(1).comps, 3u);
  ASSERT_EQ(params.getBlock(1).getStringValue("pc"), "ilu");
  ASSERT_EQ(params.getBlock(1).getIntValue("ilu_fill_level"), 1);
}


TEST(TestLinSolParams, ParseGAMG)
{
  TiXmlDocument doc;
  doc.LoadFile("src/LinAlg/Test/refdata/linsolver_gamg.xml");

  LinSolParams params;
  params.read(doc.RootElement());

  ASSERT_EQ(params.getBlock(0).getStringValue("gamg_type"), "agg");
  ASSERT_EQ(params.getBlock(0).getIntValue("gamg_proc_eq_limit"), 1000);
  ASSERT_EQ(params.getBlock(0).getIntValue("gamg_repartition"), 1);
  ASSERT_EQ(params.getBlock(0).getIntValue("gamg_use_agg_gasm"), 1);
  ASSERT_EQ(params.getBlock(0).getIntValue("gamg_reuse_interpolation"), 1);
  ASSERT_FLOAT_EQ(params.getBlock(0).getDoubleValue("gamg_threshold"), 0.1);
  ASSERT_EQ(params.getBlock(0).getStringValue("multigrid_finesmoother"), "compositedir");
  ASSERT_EQ(params.getBlock(0).dirSmoother.size(), 1u);
  ASSERT_EQ(params.getBlock(0).dirSmoother[0].type, "ilu");
  ASSERT_EQ(params.getBlock(0).dirSmoother[0].order, 12);
}


TEST(TestLinSolParams, ParseML)
{
  TiXmlDocument doc;
  doc.LoadFile("src/LinAlg/Test/refdata/linsolver_ml.xml");

  LinSolParams params;
  params.read(doc.RootElement());

  ASSERT_EQ(params.getBlock(0).getStringValue("ml_coarse_package"), "petsc");
  ASSERT_EQ(params.getBlock(0).getStringValue("ml_coarse_solver"), "lu");
  ASSERT_EQ(params.getBlock(0).getStringValue("ml_coarsen_scheme"), "agg");
  ASSERT_EQ(params.getBlock(0).getIntValue("ml_symmetrize"), 1);
  ASSERT_EQ(params.getBlock(0).getIntValue("ml_block_scaling"), 1);
  ASSERT_EQ(params.getBlock(0).getIntValue("ml_put_on_single_proc"), 1000);
  ASSERT_EQ(params.getBlock(0).getIntValue("ml_reuse_interpolation"), 1);
  ASSERT_EQ(params.getBlock(0).getIntValue("ml_reusable"), 1);
  ASSERT_EQ(params.getBlock(0).getIntValue("ml_keep_agg_info"), 1);
  ASSERT_EQ(params.getBlock(0).getIntValue("ml_aux"), 1);
  ASSERT_FLOAT_EQ(params.getBlock(0).getDoubleValue("ml_threshold"), 0.1);
  ASSERT_FLOAT_EQ(params.getBlock(0).getDoubleValue("ml_aux_threshold"), 0.05);
}


TEST(TestLinSolParams, ParseHypre)
{
  TiXmlDocument doc;
  doc.LoadFile("src/LinAlg/Test/refdata/linsolver_hypre.xml");

  LinSolParams params;
  params.read(doc.RootElement());

  ASSERT_EQ(params.getBlock(0).getStringValue("hypre_type"), "boom");
  ASSERT_EQ(params.getBlock(0).getIntValue("hypre_no_agg_coarse"), 1);
  ASSERT_EQ(params.getBlock(0).getIntValue("hypre_no_path_agg_coarse"), 2);
  ASSERT_FLOAT_EQ(params.getBlock(0).getDoubleValue("hypre_truncation"), 0.005);
  ASSERT_FLOAT_EQ(params.getBlock(0).getDoubleValue("hypre_threshold"), 0.01);
  ASSERT_EQ(params.getBlock(0).getStringValue("hypre_coarsen_scheme"), "agg");
}
