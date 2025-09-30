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

#include "Catch2Support.h"


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


TEST_CASE("TestASMs3Dmx.WriteFRTH")
{
  const ASMmxBase::MixedType param = GENERATE(
    ASMmxBase::FULL_CONT_RAISE_BASIS1,
    ASMmxBase::FULL_CONT_RAISE_BASIS2
  );
  const int basis = param == ASMmxBase::FULL_CONT_RAISE_BASIS1 ? 1 : 2;
  SECTION(basis == 1 ? "Raise basis 1" : "Raise basis 2") {
    ASMmxBase::Type = param;

    ASMbase::resetNumbering();
    ASMmxCube pch1({1,1});
    REQUIRE(pch1.generateFEMTopology());

    std::stringstream str;
    REQUIRE(pch1.write(str, 1));
    REQUIRE(str.str() == (basis == 1 ? cubeFRTH_1 : ASMCube::cube));

    str.str("");
    REQUIRE(pch1.write(str, 2));
    REQUIRE(str.str() ==  (basis == 2 ? cubeFRTH_1 : ASMCube::cube));

    REQUIRE(!pch1.write(str, 3));

    str.str("");
    REQUIRE(pch1.write(str, ASM::GEOMETRY_BASIS));
    REQUIRE(str.str() == ASMCube::cube);

    str.str("");
    REQUIRE(pch1.write(str, ASM::PROJECTION_BASIS));
    REQUIRE(str.str() == cubeFRTH_1);

    REQUIRE(!pch1.write(str, ASM::PROJECTION_BASIS_2));
    REQUIRE(!pch1.write(str, ASM::REFINEMENT_BASIS));

    str.str("");
    REQUIRE(pch1.write(str, ASM::INTEGRATION_BASIS));
    REQUIRE(str.str() == ASMCube::cube);
  }
}


TEST_CASE("TestASMs3Dmx.WriteRT")
{
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  ASMbase::resetNumbering();
  ASMmxCube pch1({1,1,1});
  pch1.raiseOrder(1,1,1);
  REQUIRE(pch1.generateFEMTopology());

  std::stringstream str;
  REQUIRE(pch1.write(str, 1));
  REQUIRE(str.str() == cubeRT_1);

  str.str("");
  REQUIRE(pch1.write(str, 2));
  REQUIRE(str.str() == cubeRT_2);

  str.str("");
  REQUIRE(pch1.write(str, 3));
  REQUIRE(str.str() == cubeRT_3);

  str.str("");
  REQUIRE(pch1.write(str, 4));
  REQUIRE(str.str() == ASMCube::cube);

  REQUIRE(!pch1.write(str, 5));

  str.str("");
  REQUIRE(pch1.write(str, ASM::GEOMETRY_BASIS));
  REQUIRE(str.str() == cubeFRTH_1);

  str.str("");
  REQUIRE(pch1.write(str, ASM::PROJECTION_BASIS));
  REQUIRE(str.str() == cubeFRTH_1);

  REQUIRE(!pch1.write(str, ASM::PROJECTION_BASIS_2));
  REQUIRE(!pch1.write(str, ASM::REFINEMENT_BASIS));

  str.str("");
  REQUIRE(pch1.write(str, ASM::INTEGRATION_BASIS));
  REQUIRE(str.str() == ASMCube::cube);
}


TEST_CASE("TestASMs3Dmx.WriteSG")
{
  ASMmxBase::Type = ASMmxBase::SUBGRID;
  ASMbase::resetNumbering();
  ASMmxCube pch1({1,1});
  REQUIRE(pch1.generateFEMTopology());

  std::stringstream str;
  REQUIRE(pch1.write(str, 1));
  REQUIRE(str.str() == cubeTH_p);

  str.str("");
  REQUIRE(pch1.write(str, 2));
  REQUIRE(str.str() == ASMCube::cube);

  REQUIRE(!pch1.write(str, 3));

  str.str("");
  REQUIRE(pch1.write(str, ASM::GEOMETRY_BASIS));
  REQUIRE(str.str() == ASMCube::cube);

  str.str("");
  REQUIRE(pch1.write(str, ASM::PROJECTION_BASIS));
  REQUIRE(str.str() == cubeTH_p);

  str.str("");
  REQUIRE(pch1.write(str, ASM::PROJECTION_BASIS_2));
  REQUIRE(str.str() == cubeFRTH_1);

  REQUIRE(!pch1.write(str, ASM::REFINEMENT_BASIS));

  str.str("");
  REQUIRE(pch1.write(str, ASM::INTEGRATION_BASIS));
  REQUIRE(str.str() == cubeTH_p);
}


TEST_CASE("TestASMs3Dmx.WriteTH")
{
  const ASMmxBase::MixedType param = GENERATE(
    ASMmxBase::REDUCED_CONT_RAISE_BASIS1,
    ASMmxBase::REDUCED_CONT_RAISE_BASIS2
  );
  const int basis = param == ASMmxBase::REDUCED_CONT_RAISE_BASIS1 ? 1 : 2;
  SECTION(basis == 1 ? "Raise basis 1" : "Raise basis 2") {
    ASMmxBase::Type = param;
    ASMbase::resetNumbering();
    ASMmxCube pch1({1,1});
    REQUIRE(pch1.uniformRefine(0, 1));
    REQUIRE(pch1.uniformRefine(1, 1));
    REQUIRE(pch1.uniformRefine(2, 1));
    REQUIRE(pch1.generateFEMTopology());

    std::stringstream str;
    REQUIRE(pch1.write(str, 1));
    REQUIRE(str.str() == (basis == 1 ? cubeTH_1 : cubeTH_2));

    str.str("");
    REQUIRE(pch1.write(str, 2));
    REQUIRE(str.str() == (basis == 2 ? cubeTH_1 : cubeTH_2));

    REQUIRE(!pch1.write(str, 3));

    str.str("");
    REQUIRE(pch1.write(str, ASM::GEOMETRY_BASIS));
    REQUIRE(str.str() == cubeTH_2);

    str.str("");
    REQUIRE(pch1.write(str, ASM::PROJECTION_BASIS));
    REQUIRE(str.str() == cubeTH_p);

    REQUIRE(!pch1.write(str, ASM::PROJECTION_BASIS_2));
    REQUIRE(!pch1.write(str, ASM::REFINEMENT_BASIS));

    str.str("");
    REQUIRE(pch1.write(str, ASM::INTEGRATION_BASIS));
    REQUIRE(str.str() == cubeTH_2);
  }
}
