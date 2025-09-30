//==============================================================================
//!
//! \file TestASMu1D.C
//!
//! \date Jun 1 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Tests for unstructured 1D FE models.
//!
//==============================================================================

#include "ASMbase.h"
#include "SIM1D.h"

#include "Catch2Support.h"


TEST_CASE("TestASMu1D.Read")
{
  std::string xml("<geometry dim='3'><patchfile type='xml'>"
                  "src/ASM/Test/refdata/bridge.xinp"
                  "</patchfile></geometry>");

  SIM1D sim;
  REQUIRE(sim.loadXML(xml.c_str()));
  REQUIRE(sim.createFEMmodel());

  ASMbase* pch1 = sim.getPatch(1);
  REQUIRE(pch1->getNoNodes() == 58);
  REQUIRE(pch1->getNoElms() == 40);

  Vec3 Xn = pch1->getCoord(56);
  REQUIRE(Xn.x == 4.25);
  REQUIRE(Xn.y == 6.25);
  REQUIRE(Xn.z ==2.0);

  Matrix Xe;
  REQUIRE(pch1->getElementCoordinates(Xe,2));
  REQUIRE(Xe.rows() == 3);
  REQUIRE(Xe.cols() == 2);
  REQUIRE(Xe(1,1) == 5.0);
  REQUIRE(Xe(2,1) == 10.0);
  REQUIRE(Xe(3,1) == 0.0);
  REQUIRE(Xe(1,2) == 2.5);
  REQUIRE(Xe(2,2) == 10.0);
  REQUIRE(Xe(3,2) == 0.0);

  IntVec s3 = {19,20,21,22,23};
  IntVec S3 = pch1->getNodeSet(pch1->getNodeSetIdx("S03"));
  REQUIRE(S3 == s3);
}
