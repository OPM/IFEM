//==============================================================================
//!
//! \file TestASMs3Dmx.C
//!
//! \date Aug 25 2023
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for structured 3D mixed spline FE models.
//!
//==============================================================================

#include "ASMCube.h"

#include "gtest/gtest.h"


namespace {

const char* cubeTH_1 =
  "700 1 0 0\n"
  "3 0\n"
  "5 3\n"
  "0 0 0 0.5 0.5 1 1 1\n"
  "5 3\n"
  "0 0 0 0.5 0.5 1 1 1\n"
  "5 3\n"
  "0 0 0 0.5 0.5 1 1 1\n"
  "0 0 0\n"
  "0.25 0 0\n"
  "0.5 0 0\n"
  "0.75 0 0\n"
  "1 0 0\n"
  "0 0.25 0\n"
  "0.25 0.25 0\n"
  "0.5 0.25 0\n"
  "0.75 0.25 0\n"
  "1 0.25 0\n"
  "0 0.5 0\n"
  "0.25 0.5 0\n"
  "0.5 0.5 0\n"
  "0.75 0.5 0\n"
  "1 0.5 0\n"
  "0 0.75 0\n"
  "0.25 0.75 0\n"
  "0.5 0.75 0\n"
  "0.75 0.75 0\n"
  "1 0.75 0\n"
  "0 1 0\n"
  "0.25 1 0\n"
  "0.5 1 0\n"
  "0.75 1 0\n"
  "1 1 0\n"
  "0 0 0.25\n"
  "0.25 0 0.25\n"
  "0.5 0 0.25\n"
  "0.75 0 0.25\n"
  "1 0 0.25\n"
  "0 0.25 0.25\n"
  "0.25 0.25 0.25\n"
  "0.5 0.25 0.25\n"
  "0.75 0.25 0.25\n"
  "1 0.25 0.25\n"
  "0 0.5 0.25\n"
  "0.25 0.5 0.25\n"
  "0.5 0.5 0.25\n"
  "0.75 0.5 0.25\n"
  "1 0.5 0.25\n"
  "0 0.75 0.25\n"
  "0.25 0.75 0.25\n"
  "0.5 0.75 0.25\n"
  "0.75 0.75 0.25\n"
  "1 0.75 0.25\n"
  "0 1 0.25\n"
  "0.25 1 0.25\n"
  "0.5 1 0.25\n"
  "0.75 1 0.25\n"
  "1 1 0.25\n"
  "0 0 0.5\n"
  "0.25 0 0.5\n"
  "0.5 0 0.5\n"
  "0.75 0 0.5\n"
  "1 0 0.5\n"
  "0 0.25 0.5\n"
  "0.25 0.25 0.5\n"
  "0.5 0.25 0.5\n"
  "0.75 0.25 0.5\n"
  "1 0.25 0.5\n"
  "0 0.5 0.5\n"
  "0.25 0.5 0.5\n"
  "0.5 0.5 0.5\n"
  "0.75 0.5 0.5\n"
  "1 0.5 0.5\n"
  "0 0.75 0.5\n"
  "0.25 0.75 0.5\n"
  "0.5 0.75 0.5\n"
  "0.75 0.75 0.5\n"
  "1 0.75 0.5\n"
  "0 1 0.5\n"
  "0.25 1 0.5\n"
  "0.5 1 0.5\n"
  "0.75 1 0.5\n"
  "1 1 0.5\n"
  "0 0 0.75\n"
  "0.25 0 0.75\n"
  "0.5 0 0.75\n"
  "0.75 0 0.75\n"
  "1 0 0.75\n"
  "0 0.25 0.75\n"
  "0.25 0.25 0.75\n"
  "0.5 0.25 0.75\n"
  "0.75 0.25 0.75\n"
  "1 0.25 0.75\n"
  "0 0.5 0.75\n"
  "0.25 0.5 0.75\n"
  "0.5 0.5 0.75\n"
  "0.75 0.5 0.75\n"
  "1 0.5 0.75\n"
  "0 0.75 0.75\n"
  "0.25 0.75 0.75\n"
  "0.5 0.75 0.75\n"
  "0.75 0.75 0.75\n"
  "1 0.75 0.75\n"
  "0 1 0.75\n"
  "0.25 1 0.75\n"
  "0.5 1 0.75\n"
  "0.75 1 0.75\n"
  "1 1 0.75\n"
  "0 0 1\n"
  "0.25 0 1\n"
  "0.5 0 1\n"
  "0.75 0 1\n"
  "1 0 1\n"
  "0 0.25 1\n"
  "0.25 0.25 1\n"
  "0.5 0.25 1\n"
  "0.75 0.25 1\n"
  "1 0.25 1\n"
  "0 0.5 1\n"
  "0.25 0.5 1\n"
  "0.5 0.5 1\n"
  "0.75 0.5 1\n"
  "1 0.5 1\n"
  "0 0.75 1\n"
  "0.25 0.75 1\n"
  "0.5 0.75 1\n"
  "0.75 0.75 1\n"
  "1 0.75 1\n"
  "0 1 1\n"
  "0.25 1 1\n"
  "0.5 1 1\n"
  "0.75 1 1\n"
  "1 1 1\n\n";

const char* cubeTH_p =
  "700 1 0 0\n"
  "3 0\n"
  "4 3\n"
  "0 0 0 0.5 1 1 1\n"
  "4 3\n"
  "0 0 0 0.5 1 1 1\n"
  "4 3\n"
  "0 0 0 0.5 1 1 1\n"
  "0 0 0\n"
  "0.25 0 0\n"
  "0.75 0 0\n"
  "1 0 0\n"
  "0 0.25 0\n"
  "0.25 0.25 0\n"
  "0.75 0.25 0\n"
  "1 0.25 0\n"
  "0 0.75 0\n"
  "0.25 0.75 0\n"
  "0.75 0.75 0\n"
  "1 0.75 0\n"
  "0 1 0\n"
  "0.25 1 0\n"
  "0.75 1 0\n"
  "1 1 0\n"
  "0 0 0.25\n"
  "0.25 0 0.25\n"
  "0.75 0 0.25\n"
  "1 0 0.25\n"
  "0 0.25 0.25\n"
  "0.25 0.25 0.25\n"
  "0.75 0.25 0.25\n"
  "1 0.25 0.25\n"
  "0 0.75 0.25\n"
  "0.25 0.75 0.25\n"
  "0.75 0.75 0.25\n"
  "1 0.75 0.25\n"
  "0 1 0.25\n"
  "0.25 1 0.25\n"
  "0.75 1 0.25\n"
  "1 1 0.25\n"
  "0 0 0.75\n"
  "0.25 0 0.75\n"
  "0.75 0 0.75\n"
  "1 0 0.75\n"
  "0 0.25 0.75\n"
  "0.25 0.25 0.75\n"
  "0.75 0.25 0.75\n"
  "1 0.25 0.75\n"
  "0 0.75 0.75\n"
  "0.25 0.75 0.75\n"
  "0.75 0.75 0.75\n"
  "1 0.75 0.75\n"
  "0 1 0.75\n"
  "0.25 1 0.75\n"
  "0.75 1 0.75\n"
  "1 1 0.75\n"
  "0 0 1\n"
  "0.25 0 1\n"
  "0.75 0 1\n"
  "1 0 1\n"
  "0 0.25 1\n"
  "0.25 0.25 1\n"
  "0.75 0.25 1\n"
  "1 0.25 1\n"
  "0 0.75 1\n"
  "0.25 0.75 1\n"
  "0.75 0.75 1\n"
  "1 0.75 1\n"
  "0 1 1\n"
  "0.25 1 1\n"
  "0.75 1 1\n"
  "1 1 1\n\n";

const char* cubeTH_2 =
  "700 1 0 0\n"
  "3 0\n"
  "3 2\n"
  "0 0 0.5 1 1\n"
  "3 2\n"
  "0 0 0.5 1 1\n"
  "3 2\n"
  "0 0 0.5 1 1\n"
  "0 0 0\n"
  "0.5 0 0\n"
  "1 0 0\n"
  "0 0.5 0\n"
  "0.5 0.5 0\n"
  "1 0.5 0\n"
  "0 1 0\n"
  "0.5 1 0\n"
  "1 1 0\n"
  "0 0 0.5\n"
  "0.5 0 0.5\n"
  "1 0 0.5\n"
  "0 0.5 0.5\n"
  "0.5 0.5 0.5\n"
  "1 0.5 0.5\n"
  "0 1 0.5\n"
  "0.5 1 0.5\n"
  "1 1 0.5\n"
  "0 0 1\n"
  "0.5 0 1\n"
  "1 0 1\n"
  "0 0.5 1\n"
  "0.5 0.5 1\n"
  "1 0.5 1\n"
  "0 1 1\n"
  "0.5 1 1\n"
  "1 1 1\n\n";

const char* cubeFRTH_1 =
  "700 1 0 0\n"
  "3 0\n"
  "3 3\n"
  "0 0 0 1 1 1\n"
  "3 3\n"
  "0 0 0 1 1 1\n"
  "3 3\n"
  "0 0 0 1 1 1\n"
  "0 0 0\n"
  "0.5 0 0\n"
  "1 0 0\n"
  "0 0.5 0\n"
  "0.5 0.5 0\n"
  "1 0.5 0\n"
  "0 1 0\n"
  "0.5 1 0\n"
  "1 1 0\n"
  "0 0 0.5\n"
  "0.5 0 0.5\n"
  "1 0 0.5\n"
  "0 0.5 0.5\n"
  "0.5 0.5 0.5\n"
  "1 0.5 0.5\n"
  "0 1 0.5\n"
  "0.5 1 0.5\n"
  "1 1 0.5\n"
  "0 0 1\n"
  "0.5 0 1\n"
  "1 0 1\n"
  "0 0.5 1\n"
  "0.5 0.5 1\n"
  "1 0.5 1\n"
  "0 1 1\n"
  "0.5 1 1\n"
  "1 1 1\n\n";

const char* cubeRT_1 =
  "700 1 0 0\n"
  "3 0\n"
  "3 3\n"
  "0 0 0 1 1 1\n"
  "2 2\n"
  "0 0 1 1\n"
  "2 2\n"
  "0 0 1 1\n"
  "0 0 0\n"
  "0.5 0 0\n"
  "1 0 0\n"
  "0 1 0\n"
  "0.5 1 0\n"
  "1 1 0\n"
  "0 0 1\n"
  "0.5 0 1\n"
  "1 0 1\n"
  "0 1 1\n"
  "0.5 1 1\n"
  "1 1 1\n\n";

const char* cubeRT_2 =
  "700 1 0 0\n"
  "3 0\n"
  "2 2\n"
  "0 0 1 1\n"
  "3 3\n"
  "0 0 0 1 1 1\n"
  "2 2\n"
  "0 0 1 1\n"
  "0 0 0\n"
  "1 0 0\n"
  "0 0.5 0\n"
  "1 0.5 0\n"
  "0 1 0\n"
  "1 1 0\n"
  "0 0 1\n"
  "1 0 1\n"
  "0 0.5 1\n"
  "1 0.5 1\n"
  "0 1 1\n"
  "1 1 1\n\n";

const char* cubeRT_3 =
  "700 1 0 0\n"
  "3 0\n"
  "2 2\n"
  "0 0 1 1\n"
  "2 2\n"
  "0 0 1 1\n"
  "3 3\n"
  "0 0 0 1 1 1\n"
  "0 0 0\n"
  "1 0 0\n"
  "0 1 0\n"
  "1 1 0\n"
  "0 0 0.5\n"
  "1 0 0.5\n"
  "0 1 0.5\n"
  "1 1 0.5\n"
  "0 0 1\n"
  "1 0 1\n"
  "0 1 1\n"
  "1 1 1\n\n";

}


class TestASMs3Dmx : public testing::Test,
                     public testing::WithParamInterface<int>
{
};


TEST_P(TestASMs3Dmx, WriteFRTH)
{
  ASMmxBase::Type = GetParam() == 0 ? ASMmxBase::FULL_CONT_RAISE_BASIS1
                                    : ASMmxBase::FULL_CONT_RAISE_BASIS2;
  ASMbase::resetNumbering();
  ASMmxCube pch1({1,1});
  EXPECT_TRUE(pch1.generateFEMTopology());

  std::stringstream str;
  EXPECT_TRUE(pch1.write(str, 1));
  EXPECT_EQ(str.str(), GetParam() == 0 ? cubeFRTH_1 : ASMCube::cube);

  str.str("");
  EXPECT_TRUE(pch1.write(str, 2));
  EXPECT_EQ(str.str(), GetParam() == 1 ? cubeFRTH_1 : ASMCube::cube);

  EXPECT_FALSE(pch1.write(str, 3));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::GEOMETRY_BASIS));
  EXPECT_EQ(str.str(), ASMCube::cube);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::PROJECTION_BASIS));
  EXPECT_EQ(str.str(), cubeFRTH_1);

  EXPECT_FALSE(pch1.write(str, ASM::PROJECTION_BASIS_2));
  EXPECT_FALSE(pch1.write(str, ASM::REFINEMENT_BASIS));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::INTEGRATION_BASIS));
  EXPECT_EQ(str.str(), ASMCube::cube);
}


TEST(TestASMs3Dmx, WriteRT)
{
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  ASMbase::resetNumbering();
  ASMmxCube pch1({1,1,1});
  EXPECT_TRUE(pch1.generateFEMTopology());

  std::stringstream str;
  EXPECT_TRUE(pch1.write(str, 1));
  EXPECT_EQ(str.str(), cubeRT_1);

  str.str("");
  EXPECT_TRUE(pch1.write(str, 2));
  EXPECT_EQ(str.str(), cubeRT_2);

  str.str("");
  EXPECT_TRUE(pch1.write(str, 3));
  EXPECT_EQ(str.str(), cubeRT_3);

  str.str("");
  EXPECT_TRUE(pch1.write(str, 4));
  EXPECT_EQ(str.str(), ASMCube::cube);

  EXPECT_FALSE(pch1.write(str, 5));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::GEOMETRY_BASIS));
  EXPECT_EQ(str.str(), ASMCube::cube);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::PROJECTION_BASIS));
  EXPECT_EQ(str.str(), cubeFRTH_1);

  EXPECT_FALSE(pch1.write(str, ASM::PROJECTION_BASIS_2));
  EXPECT_FALSE(pch1.write(str, ASM::REFINEMENT_BASIS));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::INTEGRATION_BASIS));
  EXPECT_EQ(str.str(), ASMCube::cube);
}


TEST(TestASMs3Dmx, WriteSG)
{
  ASMmxBase::Type = ASMmxBase::SUBGRID;
  ASMbase::resetNumbering();
  ASMmxCube pch1({1,1});
  EXPECT_TRUE(pch1.generateFEMTopology());

  std::stringstream str;
  EXPECT_TRUE(pch1.write(str, 1));
  EXPECT_EQ(str.str(), cubeTH_p);

  str.str("");
  EXPECT_TRUE(pch1.write(str, 2));
  EXPECT_EQ(str.str(), ASMCube::cube);

  EXPECT_FALSE(pch1.write(str, 3));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::GEOMETRY_BASIS));
  EXPECT_EQ(str.str(), ASMCube::cube);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::PROJECTION_BASIS));
  EXPECT_EQ(str.str(), cubeTH_p);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::PROJECTION_BASIS_2));
  EXPECT_EQ(str.str(), cubeFRTH_1);

  EXPECT_FALSE(pch1.write(str, ASM::REFINEMENT_BASIS));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::INTEGRATION_BASIS));
  EXPECT_EQ(str.str(), cubeTH_p);
}


TEST_P(TestASMs3Dmx, WriteTH)
{
  ASMmxBase::Type = GetParam() == 0 ? ASMmxBase::REDUCED_CONT_RAISE_BASIS1
                                    : ASMmxBase::REDUCED_CONT_RAISE_BASIS2;
  ASMbase::resetNumbering();
  ASMmxCube pch1({1,1});
  EXPECT_TRUE(pch1.uniformRefine(0, 1));
  EXPECT_TRUE(pch1.uniformRefine(1, 1));
  EXPECT_TRUE(pch1.uniformRefine(2, 1));
  EXPECT_TRUE(pch1.generateFEMTopology());

  std::stringstream str;
  EXPECT_TRUE(pch1.write(str, 1));
  EXPECT_EQ(str.str(), GetParam() == 0 ? cubeTH_1 : cubeTH_2);

  str.str("");
  EXPECT_TRUE(pch1.write(str, 2));
  EXPECT_EQ(str.str(), GetParam() == 1 ? cubeTH_1 : cubeTH_2);

  EXPECT_FALSE(pch1.write(str, 3));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::GEOMETRY_BASIS));
  EXPECT_EQ(str.str(), cubeTH_2);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::PROJECTION_BASIS));
  EXPECT_EQ(str.str(), cubeTH_p);

  EXPECT_FALSE(pch1.write(str, ASM::PROJECTION_BASIS_2));
  EXPECT_FALSE(pch1.write(str, ASM::REFINEMENT_BASIS));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::INTEGRATION_BASIS));
  EXPECT_EQ(str.str(), cubeTH_2);
}

INSTANTIATE_TEST_SUITE_P(TestASMs3Dmx, TestASMs3Dmx, testing::Values(0,1));
