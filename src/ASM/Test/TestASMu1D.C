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

#include "gtest/gtest.h"


TEST(TestASMu1D, Read)
{
  std::string xml("<geometry dim='3'><patchfile type='xml'>"
                  "  src/ASM/Test/refdata/bridge.xinp"
                  "</patchfile></geometry>");

  SIM1D sim;
  ASSERT_TRUE(sim.loadXML(xml.c_str()));
  ASSERT_TRUE(sim.createFEMmodel());

  ASMbase* pch1 = sim.getPatch(1);
  EXPECT_EQ(pch1->getNoNodes(),58u);
  EXPECT_EQ(pch1->getNoElms(),40u);

  Vec3 Xn = pch1->getCoord(56);
  EXPECT_EQ(Xn.x,4.25);
  EXPECT_EQ(Xn.y,6.25);
  EXPECT_EQ(Xn.z,2.0);

  Matrix Xe;
  EXPECT_TRUE(pch1->getElementCoordinates(Xe,2));
  EXPECT_EQ(Xe.rows(),3U);
  EXPECT_EQ(Xe.cols(),2U);
  EXPECT_EQ(Xe(1,1),5.0);
  EXPECT_EQ(Xe(2,1),10.0);
  EXPECT_EQ(Xe(3,1),0.0);
  EXPECT_EQ(Xe(1,2),2.5);
  EXPECT_EQ(Xe(2,2),10.0);
  EXPECT_EQ(Xe(3,2),0.0);

  IntVec s3 = {19,20,21,22,23};
  IntVec S3 = pch1->getNodeSet(pch1->getNodeSetIdx("S03"));
  EXPECT_TRUE(S3 == s3);
}
