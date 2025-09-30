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
#include "tinyxml2.h"

#include "Catch2Support.h"

using Catch::Matchers::WithinRel;


TEST_CASE("TestLinSolParams.ParseBlocks")
{
  tinyxml2::XMLDocument doc;
  doc.LoadFile("src/LinAlg/Test/refdata/linsolver_blocks.xml");

  LinSolParams params;
  params.read(doc.RootElement());

  // global settings
  REQUIRE(params.getNoBlocks() == 2);
  REQUIRE(params.getStringValue("pc") == "fieldsplit");
  REQUIRE(params.getStringValue("type") == "bcgs");
  REQUIRE_THAT(params.getDoubleValue("atol"), WithinRel(1e-12));
  REQUIRE_THAT(params.getDoubleValue("rtol"), WithinRel(1e-4));
  REQUIRE_THAT(params.getDoubleValue("dtol"), WithinRel(1e6));
  REQUIRE(params.getIntValue("maxits") == 123);

  // first block
  REQUIRE(params.getBlock(0).basis == 1);
  REQUIRE(params.getBlock(0).comps == 12);
  REQUIRE(params.getBlock(0).getStringValue("pc") == "sor");

  // second block
  REQUIRE(params.getBlock(1).basis == 2);
  REQUIRE(params.getBlock(1).comps == 3);
  REQUIRE(params.getBlock(1).getStringValue("pc") == "ilu");
  REQUIRE(params.getBlock(1).getIntValue("ilu_fill_level") == 1);
}


TEST_CASE("TestLinSolParams.ParseGAMG")
{
  tinyxml2::XMLDocument doc;
  doc.LoadFile("src/LinAlg/Test/refdata/linsolver_gamg.xml");

  LinSolParams params;
  params.read(doc.RootElement());

  REQUIRE(params.getBlock(0).getStringValue("gamg_type") == "agg");
  REQUIRE(params.getBlock(0).getIntValue("gamg_proc_eq_limit") == 1000);
  REQUIRE(params.getBlock(0).getIntValue("gamg_repartition") == 1);
  REQUIRE(params.getBlock(0).getIntValue("gamg_use_agg_gasm") == 1);
  REQUIRE(params.getBlock(0).getIntValue("gamg_reuse_interpolation") == 1);
  REQUIRE_THAT(params.getBlock(0).getDoubleValue("gamg_threshold"), WithinRel(0.1));
  REQUIRE(params.getBlock(0).getStringValue("multigrid_finesmoother") == "compositedir");
  REQUIRE(params.getBlock(0).dirSmoother.size() == 1);
  REQUIRE(params.getBlock(0).dirSmoother[0].type == "ilu");
  REQUIRE(params.getBlock(0).dirSmoother[0].order == 12);
}


TEST_CASE("TestLinSolParams.ParseML")
{
  tinyxml2::XMLDocument doc;
  doc.LoadFile("src/LinAlg/Test/refdata/linsolver_ml.xml");

  LinSolParams params;
  params.read(doc.RootElement());

  REQUIRE(params.getBlock(0).getStringValue("ml_coarse_package") == "petsc");
  REQUIRE(params.getBlock(0).getStringValue("ml_coarse_solver") == "lu");
  REQUIRE(params.getBlock(0).getStringValue("ml_coarsen_scheme") == "agg");
  REQUIRE(params.getBlock(0).getIntValue("ml_symmetrize") == 1);
  REQUIRE(params.getBlock(0).getIntValue("ml_block_scaling") == 1);
  REQUIRE(params.getBlock(0).getIntValue("ml_put_on_single_proc") == 1000);
  REQUIRE(params.getBlock(0).getIntValue("ml_reuse_interpolation") == 1);
  REQUIRE(params.getBlock(0).getIntValue("ml_reusable") == 1);
  REQUIRE(params.getBlock(0).getIntValue("ml_keep_agg_info") == 1);
  REQUIRE(params.getBlock(0).getIntValue("ml_aux") == 1);
  REQUIRE_THAT(params.getBlock(0).getDoubleValue("ml_threshold"), WithinRel(0.1));
  REQUIRE_THAT(params.getBlock(0).getDoubleValue("ml_aux_threshold"), WithinRel(0.05));
}


TEST_CASE("TestLinSolParams.ParseHypre")
{
  tinyxml2::XMLDocument doc;
  doc.LoadFile("src/LinAlg/Test/refdata/linsolver_hypre.xml");

  LinSolParams params;
  params.read(doc.RootElement());

  REQUIRE(params.getBlock(0).getStringValue("hypre_type") == "boom");
  REQUIRE(params.getBlock(0).getIntValue("hypre_no_agg_coarse") == 1);
  REQUIRE(params.getBlock(0).getIntValue("hypre_no_path_agg_coarse") == 2);
  REQUIRE_THAT(params.getBlock(0).getDoubleValue("hypre_truncation"), WithinRel(0.005));
  REQUIRE_THAT(params.getBlock(0).getDoubleValue("hypre_threshold"), WithinRel(0.01));
  REQUIRE(params.getBlock(0).getStringValue("hypre_coarsen_scheme") == "agg");
}
