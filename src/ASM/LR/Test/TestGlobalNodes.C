//==============================================================================
//!
//! \file TestGlobalNodes.C
//!
//! \date Jun 21 2024
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for simple global node establishment for unstructured FE models.
//!
//==============================================================================

#include "LR/GlobalNodes.h"
#include "ASMuCube.h"
#include "ASMuSquare.h"

#include "LRSpline/LRSplineVolume.h"

#include "Catch2Support.h"


TEST_CASE("TestGlobalNodes.2D")
{
  ASMuSquare pch1;
  pch1.generateFEMTopology();
  ASMuSquare pch2(2, 1.0, 0.0);
  pch2.generateFEMTopology();
  ASMuSquare pch3(2, 1.0, 1.0);
  pch3.generateFEMTopology();

  std::vector<const LR::LRSpline*> splines{pch1.getBasis(),
                                           pch2.getBasis(),
                                           pch3.getBasis()};
  std::vector<ASM::Interface> ifs;
  ifs.push_back(ASM::Interface{1, 2, 2, 1, 0, 1, 1, 0});
  ifs.push_back(ASM::Interface{2, 3, 4, 3, 0, 1, 1, 0});

  auto nodes = GlobalNodes::calcGlobalNodes(splines, ifs);
  const auto ref = std::vector{
    GlobalNodes::IntVec{0, 1, 2, 3},
    GlobalNodes::IntVec{1, 4, 3, 5},
    GlobalNodes::IntVec{3, 5, 6, 7},
  };

  REQUIRE(nodes == ref);
}


TEST_CASE("TestGlobalNodes.3D")
{
  ASMuCube pch1;
  pch1.generateFEMTopology();
  ASMuCube pch2(2, 1.0, 0.0, 0.0);
  pch2.generateFEMTopology();
  ASMuCube pch3(2, 1.0, 1.0, 1.0);
  pch3.generateFEMTopology();

  std::vector<const LR::LRSpline*> splines{pch1.getBasis(),
                                           pch2.getBasis(),
                                           pch3.getBasis()};
  std::vector<ASM::Interface> ifs;
  ifs.push_back(ASM::Interface{1, 2, 2, 1, 0, 2, 1, 0});
  ifs.push_back(ASM::Interface{2, 3, 6, 5, 0, 2, 1, 0});

  auto nodes = GlobalNodes::calcGlobalNodes(splines, ifs);

  const auto ref = std::vector{
    GlobalNodes::IntVec{0,  1, 2,  3,  4,  5,  6,  7},
    GlobalNodes::IntVec{1,  8, 3,  9,  5, 10,  7, 11},
    GlobalNodes::IntVec{5, 10, 7, 11, 12, 13, 14, 15},
  };

  REQUIRE(nodes == ref);
}
