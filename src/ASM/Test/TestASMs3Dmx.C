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

#include "GoTools/trivariate/SplineVolume.h"

#include "Catch2Support.h"

#include "ASMs3Dmx.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "IntegrandBase.h"
#include "TimeDomain.h"

#include <array>
#include <map>
#include <sstream>


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

/*!
  \brief Records the element corners.
*/

class ElementCornersIntegrand : public IntegrandBase
{
public:
  //! \brief The constructor sets the number of space dimensions.
  explicit ElementCornersIntegrand(unsigned short int n) : IntegrandBase(n) {}

  //! \brief Returns the integrand type.
  int getIntegrandType() const override { return ELEMENT_CORNERS; }

  using IntegrandBase::evalIntMx;
  //! \brief Records the element corners.
  bool evalIntMx(LocalIntegral&, const MxFiniteElement& fe,
                 const Vec3&) const override
  {
#pragma omp critical
    corners[fe.iel] = fe.XC;
    return true;
  }

  mutable std::map<int,std::vector<Vec3>> corners; //!< Corners of each element
};

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


TEST_CASE("TestASMs3Dmx.Nedelec")
{
  ASMmxBase::Type = ASMmxBase::CURL_COMPATIBLE;
  ASMbase::resetNumbering();
  ASMmxCube pch1({1,1,1,0});
  pch1.raiseOrder(1,1,1);
  REQUIRE(pch1.generateFEMTopology());

  // Component b is one order lower in direction b,
  // and the last basis is one order lower in all directions
  const int b = GENERATE(1, 2, 3, 4);
  for (int d = 0; d < 3; ++d)
    REQUIRE(pch1.getBasis(b)->order(d) == (b == 4 || b == d+1 ? 2 : 3));

  // The geometry basis is left as it is
  for (int d = 0; d < 3; ++d)
    REQUIRE(pch1.getBasis(ASM::GEOMETRY_BASIS)->order(d) == 3);

  std::stringstream str;
  REQUIRE(pch1.write(str, 4));
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


TEST_CASE("TestASMs3Dmx.ContinuousBasis")
{
  const ASMmxBase::MixedType param = GENERATE(
    ASMmxBase::DIV_COMPATIBLE,
    ASMmxBase::CURL_COMPATIBLE
  );

  SECTION(param == ASMmxBase::DIV_COMPATIBLE ? "Div" : "Curl") {
    ASMmxBase::Type = param;
    ASMbase::resetNumbering();
    ASMmxCube pch1({1,1,1});
    pch1.raiseOrder(1,1,1);
    REQUIRE(pch1.generateFEMTopology());

    // The last basis of either pair is what the derivative of its vector
    // space lands in, taken patch by patch, so it is tied to nothing
    REQUIRE(pch1.isContinuousBasis(1));
    REQUIRE(pch1.isContinuousBasis(2));
    REQUIRE(pch1.isContinuousBasis(3));
    REQUIRE(!pch1.isContinuousBasis(4));

    // And either pair holds its vector DOFs as coefficients of a reference
    // basis, which change sign with the parameter they belong to
    pch1.setPiolaMapped(true);
    REQUIRE(pch1.dofsFollowParametrization(1));
    REQUIRE(pch1.dofsFollowParametrization(2));
    REQUIRE(pch1.dofsFollowParametrization(3));
    REQUIRE(!pch1.dofsFollowParametrization(4));
  }
}


TEST_CASE("TestASMs3Dmx.ContinuousOverBoundary")
{
  const ASMmxBase::MixedType param = GENERATE(
    ASMmxBase::DIV_COMPATIBLE,
    ASMmxBase::CURL_COMPATIBLE
  );

  const bool div = param == ASMmxBase::DIV_COMPATIBLE;
  SECTION(div ? "Div" : "Curl") {
    ASMmxBase::Type = param;
    ASMbase::resetNumbering();
    ASMmxCube pch1({1,1,1});
    pch1.raiseOrder(1,1,1);
    REQUIRE(pch1.generateFEMTopology());

    // Faces 1 and 2 hold u constant and the component through them is on
    // the first basis, 3 and 4 hold v, 5 and 6 hold w. An H(div) field
    // carries the component through a face over it and an H(curl) field
    // the two along it, so the two pairs answer the one way and the other
    for (int face = 1; face <= 6; face++)
    {
      const size_t through = (face+1)/2;
      for (size_t b = 1; b <= 3; b++)
        REQUIRE(pch1.isContinuousOverBoundary(face,b) ==
                (div ? b == through : b != through));

      // Neither carries the last basis over anything
      REQUIRE(!pch1.isContinuousOverBoundary(face,4));
    }
  }

  // A pair continuous in everything it has is continuous over every face
  ASMmxBase::Type = ASMmxBase::REDUCED_CONT_RAISE_BASIS1;
  ASMbase::resetNumbering();
  ASMmxCube pch2({1,1});
  REQUIRE(pch2.generateFEMTopology());

  for (int face = 1; face <= 6; face++)
    for (size_t b = 1; b <= 2; b++)
      REQUIRE(pch2.isContinuousOverBoundary(face,b));
}


TEST_CASE("TestASMs3Dmx.SubgridCorners")
{
  // A distorted triquadratic volume, which the refined basis does not
  // represent
  std::ostringstream volume;
  volume <<"700 1 0 0\n3 0\n4 3\n0 0 0 0.4 1 1 1\n"
         <<"3 3\n0 0 0 1 1 1\n3 3\n0 0 0 1 1 1\n";
  for (double w : {0.0, 0.5, 1.0})
    for (double v : {0.0, 0.5, 1.0})
      for (double u : {0.0, 0.2, 0.7, 1.0})
        volume << u + 0.05*v*w + 0.1*v*v <<' '<< v + 0.03*u - 0.1*w*w <<' '
               << w + 0.02*u*v + 0.1*u*u <<'\n';

  ASMmxBase::Type = ASMmxBase::SUBGRID;
  ASMs3Dmx pch({3,1});
  std::istringstream is(volume.str());
  REQUIRE(pch.read(is));
  ASMbase::resetNumbering();
  REQUIRE(pch.generateFEMTopology());

  ElementCornersIntegrand itg(3);
  GlobalIntegral dummy;
  pch.generateThreadGroups(itg,true,false);
  REQUIRE(pch.integrate(itg,dummy,TimeDomain()));

  // The corners are on the geometry, at the borders of the refined elements
  const Go::SplineVolume* geo = pch.getBasis(ASM::GEOMETRY_BASIS);
  std::array<RealArray,3> knots;
  for (int d = 0; d < 3; d++)
    pch.getBasis(1)->basis(d).knotsSimple(knots[d]);
  int iel = 0;
  for (size_t k = 0; k+1 < knots[2].size(); k++)
    for (size_t j = 0; j+1 < knots[1].size(); j++)
      for (size_t i = 0; i+1 < knots[0].size(); i++)
      {
        INFO("element " << iel+1);
        const std::vector<Vec3>& XC = itg.corners[++iel];
        REQUIRE(XC.size() == 8);
        for (size_t c = 0; c < 8; c++)
        {
          Go::Point X;
          geo->point(X,knots[0][i+c%2],knots[1][j+(c/2)%2],knots[2][k+c/4]);
          for (int d = 0; d < 3; d++)
            CHECK_THAT(XC[c][d], WithinAbs(X[d], 1e-14));
        }
      }
  CHECK(itg.corners.size() == static_cast<size_t>(iel));
}
