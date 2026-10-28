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

#include "Fields.h"
#include "ItgPoint.h"

#include <memory>

#include "Catch2Support.h"


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


TEST_CASE("TestASMs2Dmx.WriteFRTH")
{
  const ASMmxBase::MixedType param = GENERATE(
    ASMmxBase::FULL_CONT_RAISE_BASIS1,
    ASMmxBase::FULL_CONT_RAISE_BASIS2
  );

  const int basis = param == ASMmxBase::FULL_CONT_RAISE_BASIS1 ? 1 : 2;
  SECTION(basis == 1 ? "Raise basis 1" : "Raise basis 2") {
    ASMmxBase::Type = param;
    ASMbase::resetNumbering();
    ASMmxSquare pch1({1,1});
    REQUIRE(pch1.generateFEMTopology());

    std::stringstream str;
    REQUIRE(pch1.write(str, 1));
    REQUIRE(str.str() == (basis == 1 ? squareFRTH_1 : ASMSquare::square));

    str.str("");
    REQUIRE(pch1.write(str, 2));
    REQUIRE(str.str() == (basis == 2 ? squareFRTH_1 : ASMSquare::square));

    REQUIRE(!pch1.write(str, 3));

    str.str("");
    REQUIRE(pch1.write(str, ASM::GEOMETRY_BASIS));
    REQUIRE(str.str() == ASMSquare::square);

    str.str("");
    REQUIRE(pch1.write(str, ASM::PROJECTION_BASIS));
    REQUIRE(str.str() == squareFRTH_1);

    REQUIRE(!pch1.write(str, ASM::PROJECTION_BASIS_2));
    REQUIRE(!pch1.write(str, ASM::REFINEMENT_BASIS));

    str.str("");
    REQUIRE(pch1.write(str, ASM::INTEGRATION_BASIS));
    REQUIRE(str.str() == ASMSquare::square);
  }
}


TEST_CASE("TestASMs2Dmx.WriteRT")
{
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  ASMbase::resetNumbering();
  ASMmxSquare pch1({1,1,1});
  pch1.raiseOrder(1,1);
  REQUIRE(pch1.generateFEMTopology());

  std::stringstream str;
  REQUIRE(pch1.write(str, 1));
  REQUIRE(str.str() == squareRT_1);

  str.str("");
  REQUIRE(pch1.write(str, 2));
  REQUIRE(str.str() == squareRT_2);

  str.str("");
  REQUIRE(pch1.write(str, 3));
  REQUIRE(str.str() == ASMSquare::square);

  REQUIRE(!pch1.write(str, 4));

  str.str("");
  REQUIRE(pch1.write(str, ASM::GEOMETRY_BASIS));
  REQUIRE(str.str() == squareFRTH_1);

  str.str("");
  REQUIRE(pch1.write(str, ASM::PROJECTION_BASIS));
  REQUIRE(str.str() == squareFRTH_1);

  REQUIRE(!pch1.write(str, ASM::PROJECTION_BASIS_2));
  REQUIRE(!pch1.write(str, ASM::REFINEMENT_BASIS));

  str.str("");
  REQUIRE(pch1.write(str, ASM::INTEGRATION_BASIS));
  REQUIRE(str.str() == ASMSquare::square);
}


TEST_CASE("TestASMs2Dmx.WriteNedelec")
{
  // The curl-compatible components are the div-compatible ones swapped
  ASMmxBase::Type = ASMmxBase::CURL_COMPATIBLE;
  ASMbase::resetNumbering();
  ASMmxSquare pch1({1,1,1});
  pch1.raiseOrder(1,1);
  REQUIRE(pch1.generateFEMTopology());

  std::stringstream str;
  REQUIRE(pch1.write(str, 1));
  REQUIRE(str.str() == squareRT_2);

  str.str("");
  REQUIRE(pch1.write(str, 2));
  REQUIRE(str.str() == squareRT_1);

  str.str("");
  REQUIRE(pch1.write(str, 3));
  REQUIRE(str.str() == ASMSquare::square);

  REQUIRE(!pch1.write(str, 4));

  str.str("");
  REQUIRE(pch1.write(str, ASM::GEOMETRY_BASIS));
  REQUIRE(str.str() == squareFRTH_1);
}


TEST_CASE("TestASMs2Dmx.WriteSG")
{
  ASMmxBase::Type = ASMmxBase::SUBGRID;
  ASMbase::resetNumbering();
  ASMmxSquare pch1({1,1});
  REQUIRE(pch1.generateFEMTopology());

  std::stringstream str;
  REQUIRE(pch1.write(str, 1));
  REQUIRE(str.str() == squareTH_p);

  str.str("");
  REQUIRE(pch1.write(str, 2));
  REQUIRE(str.str() == ASMSquare::square);

  REQUIRE(!pch1.write(str, 3));

  str.str("");
  REQUIRE(pch1.write(str, ASM::GEOMETRY_BASIS));
  REQUIRE(str.str() == ASMSquare::square);

  str.str("");
  REQUIRE(pch1.write(str, ASM::PROJECTION_BASIS));
  REQUIRE(str.str() == squareTH_p);

  str.str("");
  REQUIRE(pch1.write(str, ASM::PROJECTION_BASIS_2));
  REQUIRE(str.str() == squareFRTH_1);

  REQUIRE(!pch1.write(str, ASM::REFINEMENT_BASIS));

  str.str("");
  REQUIRE(pch1.write(str, ASM::INTEGRATION_BASIS));
  REQUIRE(str.str() == squareTH_p);
}


TEST_CASE("TestASMs2Dmx.WriteTH")
{
  const ASMmxBase::MixedType param = GENERATE(
    ASMmxBase::REDUCED_CONT_RAISE_BASIS1,
    ASMmxBase::REDUCED_CONT_RAISE_BASIS2
  );

  const int basis = param == ASMmxBase::REDUCED_CONT_RAISE_BASIS1 ? 1 : 2;
  SECTION(basis == 1 ? "Raise basis 1" : "Raise basis 2") {
    ASMmxBase::Type = param;
    ASMbase::resetNumbering();
    ASMmxSquare pch1({1,1});
    REQUIRE(pch1.uniformRefine(0, 1));
    REQUIRE(pch1.uniformRefine(1, 1));
    REQUIRE(pch1.generateFEMTopology());

    std::stringstream str;
    REQUIRE(pch1.write(str, 1));
    REQUIRE(str.str() == (basis == 1 ? squareTH_1 : squareTH_2));

    str.str("");
    REQUIRE(pch1.write(str, 2));
    REQUIRE(str.str() == (basis == 2 ? squareTH_1 : squareTH_2));

    REQUIRE(!pch1.write(str, 3));

    str.str("");
    REQUIRE(pch1.write(str, ASM::GEOMETRY_BASIS));
    REQUIRE(str.str() == squareTH_2);

    str.str("");
    REQUIRE(pch1.write(str, ASM::PROJECTION_BASIS));
    REQUIRE(str.str() == squareTH_p);

    REQUIRE(!pch1.write(str, ASM::PROJECTION_BASIS_2));
    REQUIRE(!pch1.write(str, ASM::REFINEMENT_BASIS));

    str.str("");
    REQUIRE(pch1.write(str, ASM::INTEGRATION_BASIS));
    REQUIRE(str.str() == squareTH_2);
  }
}


TEST_CASE("TestASMs2Dmx.ProjectedFieldsSubgrid")
{
  // The projection basis of the subgrid spaces is the first basis, and the
  // projected fields are used by mixed multi-patch models
  const ASMmxBase::MixedType oldType = ASMmxBase::Type;
  ASMmxBase::Type = ASMmxBase::SUBGRID;
  ASMbase::resetNumbering();
  ASMmxSquare pch({2,1});
  REQUIRE(pch.uniformRefine(0,1));
  REQUIRE(pch.uniformRefine(1,1));
  REQUIRE(pch.generateFEMTopology());
  REQUIRE(pch.getNoProjectionNodes() == pch.getNoNodes(1));

  // The x-coordinates of the nodes of the first basis
  Vector coefs(pch.getNoNodes(1));
  for (size_t n = 1; n <= coefs.size(); n++)
    coefs(n) = pch.getCoord(n).x;

  std::unique_ptr<Fields> f(pch.getProjectedFields(coefs,1));
  REQUIRE(f);
  Vector vals;
  REQUIRE(f->valueFE(ItgPoint(0.3,0.7),vals));
  REQUIRE(vals.size() == 1);
  CHECK_THAT(vals(1), WithinAbs(0.3, 1e-12));

  ASMmxBase::Type = oldType;
}
