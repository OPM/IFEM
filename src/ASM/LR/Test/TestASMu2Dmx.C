//==============================================================================
//!
//! \file TestASMu2Dmx.C
//!
//! \date Aug 25 2023
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for structured 2D mixed LR spline FE models.
//!
//==============================================================================

#include "ASMmxBase.h"
#include "ASMuSquare.h"

#include "gtest/gtest.h"


namespace {

// formatting matches write routine, do not change
const char* squareTH_1 =
  "# LRSPLINE SURFACE\n"
  "#\tp1\tp2\tNbasis\tNline\tNel\tdim\trat\n"
  "\t3\t3\t25\t6\t4\t2\t0\n"
  "# Basis functions:\n"
  "0: [0 0 0 0.5 ] x [0 0 0 0.5 ] 0 0 (1)\n"
  "1: [0.5 1 1 1 ] x [0 0 0 0.5 ] 1 0 (1)\n"
  "2: [0.5 0.5 1 1 ] x [0 0 0 0.5 ] 0.75 0 (1)\n"
  "3: [0 0 0.5 0.5 ] x [0 0 0 0.5 ] 0.25 0 (1)\n"
  "4: [0 0.5 0.5 1 ] x [0 0 0 0.5 ] 0.5 0 (1)\n"
  "5: [0 0 0 0.5 ] x [0.5 1 1 1 ] 0 1 (1)\n"
  "6: [0.5 1 1 1 ] x [0.5 1 1 1 ] 1 1 (1)\n"
  "7: [0.5 0.5 1 1 ] x [0.5 1 1 1 ] 0.75 1 (1)\n"
  "8: [0 0 0.5 0.5 ] x [0.5 1 1 1 ] 0.25 1 (1)\n"
  "9: [0 0.5 0.5 1 ] x [0.5 1 1 1 ] 0.5 1 (1)\n"
  "10: [0 0 0 0.5 ] x [0.5 0.5 1 1 ] 0 0.75 (1)\n"
  "11: [0.5 1 1 1 ] x [0.5 0.5 1 1 ] 1 0.75 (1)\n"
  "12: [0.5 0.5 1 1 ] x [0.5 0.5 1 1 ] 0.75 0.75 (1)\n"
  "13: [0 0 0.5 0.5 ] x [0.5 0.5 1 1 ] 0.25 0.75 (1)\n"
  "14: [0 0.5 0.5 1 ] x [0.5 0.5 1 1 ] 0.5 0.75 (1)\n"
  "15: [0 0 0 0.5 ] x [0 0 0.5 0.5 ] 0 0.25 (1)\n"
  "16: [0.5 1 1 1 ] x [0 0 0.5 0.5 ] 1 0.25 (1)\n"
  "17: [0.5 0.5 1 1 ] x [0 0 0.5 0.5 ] 0.75 0.25 (1)\n"
  "18: [0 0 0.5 0.5 ] x [0 0 0.5 0.5 ] 0.25 0.25 (1)\n"
  "19: [0 0.5 0.5 1 ] x [0 0 0.5 0.5 ] 0.5 0.25 (1)\n"
  "20: [0 0 0 0.5 ] x [0 0.5 0.5 1 ] 0 0.5 (1)\n"
  "21: [0.5 1 1 1 ] x [0 0.5 0.5 1 ] 1 0.5 (1)\n"
  "22: [0.5 0.5 1 1 ] x [0 0.5 0.5 1 ] 0.75 0.5 (1)\n"
  "23: [0 0 0.5 0.5 ] x [0 0.5 0.5 1 ] 0.25 0.5 (1)\n"
  "24: [0 0.5 0.5 1 ] x [0 0.5 0.5 1 ] 0.5 0.5 (1)\n"
  "# Mesh lines:\n"
  "0 x [0, 1] (3)\n"
  "0.5 x [0, 1] (2)\n"
  "1 x [0, 1] (3)\n"
  "[0, 1] x 0 (3)\n"
  "[0, 1] x 0.5 (2)\n"
  "[0, 1] x 1 (3)\n"
  "# Elements:\n"
  "0 [2] : (0, 0) x (0.5, 0.5)    {0, 3, 4, 15, 18, 19, 20, 23, 24}\n"
  "1 [2] : (0.5, 0) x (1, 0.5)    {1, 2, 4, 16, 17, 19, 21, 22, 24}\n"
  "2 [2] : (0, 0.5) x (0.5, 1)    {5, 8, 9, 10, 13, 14, 20, 23, 24}\n"
  "3 [2] : (0.5, 0.5) x (1, 1)    {6, 7, 9, 11, 12, 14, 21, 22, 24}\n";

// formatting matches write routine, do not change
const char* squareTH_2 =
  "# LRSPLINE SURFACE\n"
  "#\tp1\tp2\tNbasis\tNline\tNel\tdim\trat\n"
  "\t2\t2\t9\t6\t4\t2\t0\n"
  "# Basis functions:\n"
  "0: [0 0 0.5 ] x [0 0 0.5 ] 0 0 (1)\n"
  "1: [0.5 1 1 ] x [0 0 0.5 ] 1 0 (1)\n"
  "2: [0 0.5 1 ] x [0 0 0.5 ] 0.5 0 (1)\n"
  "3: [0 0 0.5 ] x [0.5 1 1 ] 0 1 (1)\n"
  "4: [0.5 1 1 ] x [0.5 1 1 ] 1 1 (1)\n"
  "5: [0 0.5 1 ] x [0.5 1 1 ] 0.5 1 (1)\n"
  "6: [0 0 0.5 ] x [0 0.5 1 ] 0 0.5 (1)\n"
  "7: [0.5 1 1 ] x [0 0.5 1 ] 1 0.5 (1)\n"
  "8: [0 0.5 1 ] x [0 0.5 1 ] 0.5 0.5 (1)\n"
  "# Mesh lines:\n"
  "0 x [0, 1] (2)\n"
  "0.5 x [0, 1] (1)\n"
  "1 x [0, 1] (2)\n"
  "[0, 1] x 0 (2)\n"
  "[0, 1] x 0.5 (1)\n"
  "[0, 1] x 1 (2)\n"
  "# Elements:\n"
  "0 [2] : (0, 0) x (0.5, 0.5)    {0, 2, 6, 8}\n"
  "1 [2] : (0.5, 0) x (1, 0.5)    {1, 2, 7, 8}\n"
  "2 [2] : (0, 0.5) x (0.5, 1)    {3, 5, 6, 8}\n"
  "3 [2] : (0.5, 0.5) x (1, 1)    {4, 5, 7, 8}\n";

const char* squareTH_p =
  "# LRSPLINE SURFACE\n"
  "#\tp1\tp2\tNbasis\tNline\tNel\tdim\trat\n"
  "\t3\t3\t16\t6\t4\t2\t0\n"
  "# Basis functions:\n"
  "0: [0 0 0 0.5 ] x [0 0 0 0.5 ] 0 0 (1)\n"
  "1: [0.5 1 1 1 ] x [0 0 0 0.5 ] 1 0 (1)\n"
  "2: [0 0.5 1 1 ] x [0 0 0 0.5 ] 0.75 0 (1)\n"
  "3: [0 0 0.5 1 ] x [0 0 0 0.5 ] 0.25 0 (1)\n"
  "4: [0 0 0 0.5 ] x [0.5 1 1 1 ] 0 1 (1)\n"
  "5: [0.5 1 1 1 ] x [0.5 1 1 1 ] 1 1 (1)\n"
  "6: [0 0.5 1 1 ] x [0.5 1 1 1 ] 0.75 1 (1)\n"
  "7: [0 0 0.5 1 ] x [0.5 1 1 1 ] 0.25 1 (1)\n"
  "8: [0 0 0 0.5 ] x [0 0.5 1 1 ] 0 0.75 (1)\n"
  "9: [0.5 1 1 1 ] x [0 0.5 1 1 ] 1 0.75 (1)\n"
  "10: [0 0.5 1 1 ] x [0 0.5 1 1 ] 0.75 0.75 (1)\n"
  "11: [0 0 0.5 1 ] x [0 0.5 1 1 ] 0.25 0.75 (1)\n"
  "12: [0 0 0 0.5 ] x [0 0 0.5 1 ] 0 0.25 (1)\n"
  "13: [0.5 1 1 1 ] x [0 0 0.5 1 ] 1 0.25 (1)\n"
  "14: [0 0.5 1 1 ] x [0 0 0.5 1 ] 0.75 0.25 (1)\n"
  "15: [0 0 0.5 1 ] x [0 0 0.5 1 ] 0.25 0.25 (1)\n"
  "# Mesh lines:\n"
  "0 x [0, 1] (3)\n"
  "0.5 x [0, 1] (1)\n"
  "1 x [0, 1] (3)\n"
  "[0, 1] x 0 (3)\n"
  "[0, 1] x 0.5 (1)\n"
  "[0, 1] x 1 (3)\n"
  "# Elements:\n"
  "0 [2] : (0, 0) x (0.5, 0.5)    {0, 2, 3, 8, 10, 11, 12, 14, 15}\n"
  "1 [2] : (0.5, 0) x (1, 0.5)    {1, 2, 3, 9, 10, 11, 13, 14, 15}\n"
  "2 [2] : (0, 0.5) x (0.5, 1)    {4, 6, 7, 8, 10, 11, 12, 14, 15}\n"
  "3 [2] : (0.5, 0.5) x (1, 1)    {5, 6, 7, 9, 10, 11, 13, 14, 15}\n";

const char* squareFRTH_1 =
  "# LRSPLINE SURFACE\n"
  "#\tp1\tp2\tNbasis\tNline\tNel\tdim\trat\n"
  "\t3\t3\t9\t4\t1\t2\t0\n"
  "# Basis functions:\n"
  "0: [0 0 0 1 ] x [0 0 0 1 ] 0 0 (1)\n"
  "1: [0 0 1 1 ] x [0 0 0 1 ] 0.5 0 (1)\n"
  "2: [0 1 1 1 ] x [0 0 0 1 ] 1 0 (1)\n"
  "3: [0 0 0 1 ] x [0 0 1 1 ] 0 0.5 (1)\n"
  "4: [0 0 1 1 ] x [0 0 1 1 ] 0.5 0.5 (1)\n"
  "5: [0 1 1 1 ] x [0 0 1 1 ] 1 0.5 (1)\n"
  "6: [0 0 0 1 ] x [0 1 1 1 ] 0 1 (1)\n"
  "7: [0 0 1 1 ] x [0 1 1 1 ] 0.5 1 (1)\n"
  "8: [0 1 1 1 ] x [0 1 1 1 ] 1 1 (1)\n"
  "# Mesh lines:\n0 x [0, 1] (3)\n"
  "1 x [0, 1] (3)\n"
  "[0, 1] x 0 (3)\n"
  "[0, 1] x 1 (3)\n"
  "# Elements:\n"
  "0 [2] : (0, 0) x (1, 1)    {0, 1, 2, 3, 4, 5, 6, 7, 8}\n";

const char* squareRT_1 =
  "# LRSPLINE SURFACE\n"
  "#\tp1\tp2\tNbasis\tNline\tNel\tdim\trat\n"
  "\t3\t2\t6\t4\t1\t2\t0\n"
  "# Basis functions:\n"
  "0: [0 0 0 1 ] x [0 0 1 ] 0 0 (1)\n"
  "1: [0 0 1 1 ] x [0 0 1 ] 0.5 0 (1)\n"
  "2: [0 1 1 1 ] x [0 0 1 ] 1 0 (1)\n"
  "3: [0 0 0 1 ] x [0 1 1 ] 0 1 (1)\n"
  "4: [0 0 1 1 ] x [0 1 1 ] 0.5 1 (1)\n"
  "5: [0 1 1 1 ] x [0 1 1 ] 1 1 (1)\n"
  "# Mesh lines:\n"
  "0 x [0, 1] (3)\n"
  "1 x [0, 1] (3)\n"
  "[0, 1] x 0 (2)\n"
  "[0, 1] x 1 (2)\n"
  "# Elements:\n"
  "0 [2] : (0, 0) x (1, 1)    {0, 1, 2, 3, 4, 5}\n";

const char* squareRT_2 =
  "# LRSPLINE SURFACE\n"
  "#\tp1\tp2\tNbasis\tNline\tNel\tdim\trat\n"
  "\t2\t3\t6\t4\t1\t2\t0\n"
  "# Basis functions:\n"
  "0: [0 0 1 ] x [0 0 0 1 ] 0 0 (1)\n"
  "1: [0 1 1 ] x [0 0 0 1 ] 1 0 (1)\n"
  "2: [0 0 1 ] x [0 0 1 1 ] 0 0.5 (1)\n"
  "3: [0 1 1 ] x [0 0 1 1 ] 1 0.5 (1)\n"
  "4: [0 0 1 ] x [0 1 1 1 ] 0 1 (1)\n"
  "5: [0 1 1 ] x [0 1 1 1 ] 1 1 (1)\n"
  "# Mesh lines:\n0 x [0, 1] (2)\n"
  "1 x [0, 1] (2)\n"
  "[0, 1] x 0 (3)\n"
  "[0, 1] x 1 (3)\n"
  "# Elements:\n"
  "0 [2] : (0, 0) x (1, 1)    {0, 1, 2, 3, 4, 5}\n";

}


class TestASMu2Dmx : public testing::Test,
                     public testing::WithParamInterface<int>
{
};


TEST_P(TestASMu2Dmx, WriteFRTH)
{
    ASMmxBase::Type = GetParam() == 0 ? ASMmxBase::FULL_CONT_RAISE_BASIS1
                                      : ASMmxBase::FULL_CONT_RAISE_BASIS2;
  ASMbase::resetNumbering();
  ASMmxuSquare pch1({1,1});
  EXPECT_TRUE(pch1.generateFEMTopology());

  std::stringstream str;
  EXPECT_TRUE(pch1.write(str, 1));
  EXPECT_EQ(str.str(), GetParam() == 0 ? squareFRTH_1 : ASMuSquare::square);

  str.str("");
  EXPECT_TRUE(pch1.write(str, 2));
  EXPECT_EQ(str.str(), GetParam() == 1 ? squareFRTH_1 : ASMuSquare::square);

  EXPECT_FALSE(pch1.write(str, 3));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::GEOMETRY_BASIS));
  EXPECT_EQ(str.str(), ASMuSquare::square);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::PROJECTION_BASIS));
  EXPECT_EQ(str.str(), squareFRTH_1);

  EXPECT_FALSE(pch1.write(str, ASM::PROJECTION_BASIS_2));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::REFINEMENT_BASIS));
  EXPECT_EQ(str.str(), squareFRTH_1);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::INTEGRATION_BASIS));
  EXPECT_EQ(str.str(), ASMuSquare::square);
}

TEST(TestASMu2Dmx, WriteRT)
{
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  ASMbase::resetNumbering();
  ASMmxuSquare pch1({1,1,1});
  EXPECT_TRUE(pch1.generateFEMTopology());

  std::stringstream str;
  EXPECT_TRUE(pch1.write(str, 1));
  EXPECT_EQ(str.str(), squareRT_1);

  str.str("");
  EXPECT_TRUE(pch1.write(str, 2));
  EXPECT_EQ(str.str(), squareRT_2);

  str.str("");
  EXPECT_TRUE(pch1.write(str, 3));
  EXPECT_EQ(str.str(), ASMuSquare::square);

  EXPECT_FALSE(pch1.write(str, 4));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::GEOMETRY_BASIS));
  EXPECT_EQ(str.str(), ASMuSquare::square);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::PROJECTION_BASIS));
  EXPECT_EQ(str.str(), squareFRTH_1);

  EXPECT_FALSE(pch1.write(str, ASM::PROJECTION_BASIS_2));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::REFINEMENT_BASIS));
  EXPECT_EQ(str.str(), squareFRTH_1);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::INTEGRATION_BASIS));
  EXPECT_EQ(str.str(), ASMuSquare::square);
}

TEST(TestASMu2Dmx, WriteSG)
{
  ASMmxBase::Type = ASMmxBase::SUBGRID;
  ASMbase::resetNumbering();
  ASMmxuSquare pch1({1,1});
  EXPECT_TRUE(pch1.generateFEMTopology());

  std::stringstream str;
  EXPECT_TRUE(pch1.write(str, 1));
  EXPECT_EQ(str.str(), squareTH_p);

  str.str("");
  EXPECT_TRUE(pch1.write(str, 2));
  EXPECT_EQ(str.str(), ASMuSquare::square);

  EXPECT_FALSE(pch1.write(str, 3));

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::GEOMETRY_BASIS));
  EXPECT_EQ(str.str(), ASMuSquare::square);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::PROJECTION_BASIS));
  EXPECT_EQ(str.str(), squareTH_p);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::PROJECTION_BASIS_2));
  EXPECT_EQ(str.str(), squareFRTH_1);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::REFINEMENT_BASIS));
  EXPECT_EQ(str.str(), squareFRTH_1);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::INTEGRATION_BASIS));
  EXPECT_EQ(str.str(), squareTH_p);
}


TEST_P(TestASMu2Dmx, WriteTH)
{
  ASMmxBase::Type = GetParam() == 0 ? ASMmxBase::REDUCED_CONT_RAISE_BASIS1
                                    : ASMmxBase::REDUCED_CONT_RAISE_BASIS2;
  ASMbase::resetNumbering();
  ASMmxuSquare pch1({1,1});
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

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::REFINEMENT_BASIS));
  EXPECT_EQ(str.str(), squareTH_p);

  str.str("");
  EXPECT_TRUE(pch1.write(str, ASM::INTEGRATION_BASIS));
  EXPECT_EQ(str.str(), squareTH_2);
}


INSTANTIATE_TEST_SUITE_P(TestASMu2Dmx, TestASMu2Dmx, testing::Values(0,1));
