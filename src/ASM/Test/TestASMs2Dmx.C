//==============================================================================
//!
//! \file TestASMs2Dmx.C
//!
//! \date Aug 25 2023
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for structured 2D mixed spline FE models.
//!
//==============================================================================

#include "ASMmxBase.h"
#include "ASMSquare.h"

#include "gtest/gtest.h"


namespace {

// formatting matches write routine, do not change
const char* squareTH_1 =
  "200 1 0 0\n"
  "2 0\n"
  "5 3\n"
  "0 0 0 0.5 0.5 1 1 1\n"
  "5 3\n"
  "0 0 0 0.5 0.5 1 1 1\n"
  "0 0\n"
  "0.25 0\n"
  "0.5 0\n"
  "0.75 0\n"
  "1 0\n"
  "0 0.25\n"
  "0.25 0.25\n"
  "0.5 0.25\n"
  "0.75 0.25\n"
  "1 0.25\n"
  "0 0.5\n"
  "0.25 0.5\n"
  "0.5 0.5\n"
  "0.75 0.5\n"
  "1 0.5\n"
  "0 0.75\n"
  "0.25 0.75\n"
  "0.5 0.75\n"
  "0.75 0.75\n"
  "1 0.75\n"
  "0 1\n"
  "0.25 1\n"
  "0.5 1\n"
  "0.75 1\n"
  "1 1\n\n";

// formatting matches write routine, do not change
const char* squareTH_2 =
  "200 1 0 0\n"
  "2 0\n"
  "3 2\n"
  "0 0 0.5 1 1\n"
  "3 2\n"
  "0 0 0.5 1 1\n"
  "0 0\n"
  "0.5 0\n"
  "1 0\n"
  "0 0.5\n"
  "0.5 0.5\n"
  "1 0.5\n"
  "0 1\n"
  "0.5 1\n"
  "1 1\n\n";

const char* squareTH_p =
  "200 1 0 0\n"
  "2 0\n"
  "4 3\n"
  "0 0 0 0.5 1 1 1\n"
  "4 3\n"
  "0 0 0 0.5 1 1 1\n"
  "0 0\n"
  "0.25 0\n"
  "0.75 0\n"
  "1 0\n"
  "0 0.25\n"
  "0.25 0.25\n"
  "0.75 0.25\n"
  "1 0.25\n"
  "0 0.75\n"
  "0.25 0.75\n"
  "0.75 0.75\n"
  "1 0.75\n"
  "0 1\n"
  "0.25 1\n"
  "0.75 1\n"
  "1 1\n\n";

const char* squareFRTH_1 =
  "200 1 0 0\n"
  "2 0\n"
  "3 3\n"
  "0 0 0 1 1 1\n"
  "3 3\n"
  "0 0 0 1 1 1\n"
  "0 0\n"
  "0.5 0\n"
  "1 0\n"
  "0 0.5\n"
  "0.5 0.5\n"
  "1 0.5\n"
  "0 1\n"
  "0.5 1\n"
  "1 1\n\n";

const char* squareRT_1 =
  "200 1 0 0\n"
  "2 0\n"
  "3 3\n"
  "0 0 0 1 1 1\n"
  "2 2\n"
  "0 0 1 1\n"
  "0 0\n"
  "0.5 0\n"
  "1 0\n"
  "0 1\n"
  "0.5 1\n"
  "1 1\n\n";

const char* squareRT_2 =
  "200 1 0 0\n"
  "2 0\n"
  "2 2\n"
  "0 0 1 1\n"
  "3 3\n"
  "0 0 0 1 1 1\n"
  "0 0\n"
  "1 0\n"
  "0 0.5\n"
  "1 0.5\n"
  "0 1\n"
  "1 1\n\n";

}


class TestASMs2Dmx : public testing::Test,
                     public testing::WithParamInterface<int>
{
};


TEST_P(TestASMs2Dmx, WriteFRTH)
{
  ASMmxBase::Type = GetParam() == 0 ? ASMmxBase::FULL_CONT_RAISE_BASIS1
                                    : ASMmxBase::FULL_CONT_RAISE_BASIS2;
  ASMbase::resetNumbering();
  ASMmxSquare pch1({1,1});
  EXPECT_TRUE(pch1.generateFEMTopology());

  std::stringstream str;
  EXPECT_TRUE(pch1.write(str, 1));
  EXPECT_EQ(str.str(), GetParam() == 0 ? squareFRTH_1 : ASMSquare::square);

  str.str("");
  EXPECT_TRUE(pch1.write(str, 2));
  EXPECT_EQ(str.str(), GetParam() == 1 ? squareFRTH_1 : ASMSquare::square);

  EXPECT_FALSE(pch1.write(str, 3));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::GEOMETRY_BASIS));
  EXPECT_EQ(str.str(), ASMSquare::square);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::PROJECTION_BASIS));
  EXPECT_EQ(str.str(), squareFRTH_1);

  EXPECT_FALSE(pch1.write(str, ASM::PROJECTION_BASIS_2));
  EXPECT_FALSE(pch1.write(str, ASM::REFINEMENT_BASIS));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::INTEGRATION_BASIS));
  EXPECT_EQ(str.str(), ASMSquare::square);
}


TEST(TestASMs2Dmx, WriteRT)
{
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  ASMbase::resetNumbering();
  ASMmxSquare pch1({1,1,1});
  EXPECT_TRUE(pch1.generateFEMTopology());

  std::stringstream str;
  EXPECT_TRUE(pch1.write(str, 1));
  EXPECT_EQ(str.str(), squareRT_1);

  str.str("");
  EXPECT_TRUE(pch1.write(str, 2));
  EXPECT_EQ(str.str(), squareRT_2);

  str.str("");
  EXPECT_TRUE(pch1.write(str, 3));
  EXPECT_EQ(str.str(), ASMSquare::square);

  EXPECT_FALSE(pch1.write(str, 4));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::GEOMETRY_BASIS));
  EXPECT_EQ(str.str(), squareFRTH_1);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::PROJECTION_BASIS));
  EXPECT_EQ(str.str(), squareFRTH_1);

  EXPECT_FALSE(pch1.write(str, ASM::PROJECTION_BASIS_2));
  EXPECT_FALSE(pch1.write(str, ASM::REFINEMENT_BASIS));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::INTEGRATION_BASIS));
  EXPECT_EQ(str.str(), ASMSquare::square);
}


TEST(TestASMs2Dmx, WriteSG)
{
  ASMmxBase::Type = ASMmxBase::SUBGRID;
  ASMbase::resetNumbering();
  ASMmxSquare pch1({1,1});
  EXPECT_TRUE(pch1.generateFEMTopology());

  std::stringstream str;
  EXPECT_TRUE(pch1.write(str, 1));
  EXPECT_EQ(str.str(), squareTH_p);

  str.str("");
  EXPECT_TRUE(pch1.write(str, 2));
  EXPECT_EQ(str.str(), ASMSquare::square);

  EXPECT_FALSE(pch1.write(str, 3));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::GEOMETRY_BASIS));
  EXPECT_EQ(str.str(), ASMSquare::square);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::PROJECTION_BASIS));
  EXPECT_EQ(str.str(), squareTH_p);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::PROJECTION_BASIS_2));
  EXPECT_EQ(str.str(), squareFRTH_1);

  EXPECT_FALSE(pch1.write(str, ASM::REFINEMENT_BASIS));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::INTEGRATION_BASIS));
  EXPECT_EQ(str.str(), squareTH_p);
}


TEST_P(TestASMs2Dmx, WriteTH)
{
  ASMmxBase::Type = GetParam() == 0 ? ASMmxBase::REDUCED_CONT_RAISE_BASIS1
                                    : ASMmxBase::REDUCED_CONT_RAISE_BASIS2;
  ASMbase::resetNumbering();
  ASMmxSquare pch1({1,1});
  EXPECT_TRUE(pch1.uniformRefine(0, 1));
  EXPECT_TRUE(pch1.uniformRefine(1, 1));
  EXPECT_TRUE(pch1.generateFEMTopology());

  std::stringstream str;
  EXPECT_TRUE(pch1.write(str, 1));
  EXPECT_EQ(str.str(), GetParam() == 0 ? squareTH_1 : squareTH_2);

  str.str("");
  EXPECT_TRUE(pch1.write(str, 2));
  EXPECT_EQ(str.str(), GetParam() == 1 ? squareTH_1 : squareTH_2);

  EXPECT_FALSE(pch1.write(str, 3));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::GEOMETRY_BASIS));
  EXPECT_EQ(str.str(), squareTH_2);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::PROJECTION_BASIS));
  EXPECT_EQ(str.str(), squareTH_p);

  EXPECT_FALSE(pch1.write(str, ASM::PROJECTION_BASIS_2));
  EXPECT_FALSE(pch1.write(str, ASM::REFINEMENT_BASIS));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::INTEGRATION_BASIS));
  EXPECT_EQ(str.str(), squareTH_2);
}

INSTANTIATE_TEST_SUITE_P(TestASMs2Dmx, TestASMs2Dmx, testing::Values(0,1));
