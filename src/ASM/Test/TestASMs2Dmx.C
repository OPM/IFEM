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

#include "Catch2Support.h"

#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "IntegrandBase.h"
#include "TimeDomain.h"
#include "GoTools/geometry/SplineSurface.h"

#include <map>
#include <sstream>


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


TEST_CASE("TestASMs2Dmx.ContinuousBasis")
{
  const ASMmxBase::MixedType param = GENERATE(
    ASMmxBase::DIV_COMPATIBLE,
    ASMmxBase::CURL_COMPATIBLE
  );

  SECTION(param == ASMmxBase::DIV_COMPATIBLE ? "Div" : "Curl") {
    ASMmxBase::Type = param;
    ASMbase::resetNumbering();
    ASMmxSquare pch1({1,1,1});
    pch1.raiseOrder(1,1);
    REQUIRE(pch1.generateFEMTopology());

    // The last basis of either pair is what the derivative of its vector
    // space lands in, taken patch by patch, so it is tied to nothing
    REQUIRE(pch1.isContinuousBasis(1));
    REQUIRE(pch1.isContinuousBasis(2));
    REQUIRE(!pch1.isContinuousBasis(3));

    // And either pair holds its vector DOFs as coefficients of a reference
    // basis, which change sign with the parameter they belong to
    pch1.setPiolaMapped(true);
    REQUIRE(pch1.dofsFollowParametrization(1));
    REQUIRE(pch1.dofsFollowParametrization(2));
    REQUIRE(!pch1.dofsFollowParametrization(3));
  }

  // A pair which is continuous in everything it has keeps every basis
  ASMmxBase::Type = ASMmxBase::REDUCED_CONT_RAISE_BASIS1;
  ASMbase::resetNumbering();
  ASMmxSquare pch2({1,1});
  REQUIRE(pch2.generateFEMTopology());

  pch2.setPiolaMapped(true);
  for (size_t b = 1; b <= 2; b++)
  {
    REQUIRE(pch2.isContinuousBasis(b));
    REQUIRE(!pch2.dofsFollowParametrization(b));
  }
}


TEST_CASE("TestASMs2Dmx.ContinuousOverBoundary")
{
  const ASMmxBase::MixedType param = GENERATE(
    ASMmxBase::DIV_COMPATIBLE,
    ASMmxBase::CURL_COMPATIBLE
  );

  const bool div = param == ASMmxBase::DIV_COMPATIBLE;
  SECTION(div ? "Div" : "Curl") {
    ASMmxBase::Type = param;
    ASMbase::resetNumbering();
    ASMmxSquare pch1({1,1,1});
    pch1.raiseOrder(1,1);
    REQUIRE(pch1.generateFEMTopology());

    // Edges 1 and 2 hold u constant and the component through them is on
    // the first basis, edges 3 and 4 hold v and theirs is the second. An
    // H(div) field carries the component through an edge over it and an
    // H(curl) field the one along it, so the two pairs answer the one way
    // and the other
    for (int edge = 1; edge <= 4; edge++)
    {
      const size_t through = (edge+1)/2;
      for (size_t b = 1; b <= 2; b++)
        REQUIRE(pch1.isContinuousOverBoundary(edge,b) ==
                (div ? b == through : b != through));

      // Neither carries the last basis over anything
      REQUIRE(!pch1.isContinuousOverBoundary(edge,3));
    }
  }

  // A pair continuous in everything it has is continuous over every edge
  ASMmxBase::Type = ASMmxBase::REDUCED_CONT_RAISE_BASIS1;
  ASMbase::resetNumbering();
  ASMmxSquare pch2({1,1});
  REQUIRE(pch2.generateFEMTopology());

  for (int edge = 1; edge <= 4; edge++)
    for (size_t b = 1; b <= 2; b++)
      REQUIRE(pch2.isContinuousOverBoundary(edge,b));
}


TEST_CASE("TestASMs2Dmx.SubgridCorners")
{
  // A curved biquadratic surface, which the refined basis does not represent
  const char* surface = "200 1 0 0\n2 0\n"
                        "4 3\n0 0 0 0.3 1 1 1\n"
                        "4 3\n0 0 0 0.6 1 1 1\n"
                        "0 0\n0.4 0.1\n1.1 -0.1\n1.5 0\n"
                        "0 0.5\n0.5 0.6\n1 0.4\n1.6 0.5\n"
                        "0.1 1.2\n0.4 1.1\n1.1 1.3\n1.5 1.2\n"
                        "0 2\n0.6 2.1\n1.2 1.9\n1.6 2\n";

  ASMmxBase::Type = ASMmxBase::SUBGRID;
  ASMs2Dmx pch(2,{2,1});
  std::istringstream is(surface);
  REQUIRE(pch.read(is));
  ASMbase::resetNumbering();
  REQUIRE(pch.generateFEMTopology());

  ElementCornersIntegrand itg(2);
  GlobalIntegral dummy;
  pch.generateThreadGroups(itg,true,false);
  REQUIRE(pch.integrate(itg,dummy,TimeDomain()));

  // The corners are on the geometry, at the borders of the refined elements
  const Go::SplineSurface* geo = pch.getBasis(ASM::GEOMETRY_BASIS);
  RealArray ku, kv;
  pch.getBasis(1)->basis_u().knotsSimple(ku);
  pch.getBasis(1)->basis_v().knotsSimple(kv);
  REQUIRE(itg.corners.size() == (ku.size()-1)*(kv.size()-1));
  int iel = 0;
  for (size_t j = 0; j+1 < kv.size(); j++)
    for (size_t i = 0; i+1 < ku.size(); i++)
    {
      INFO("element " << iel+1);
      const std::vector<Vec3>& XC = itg.corners[++iel];
      REQUIRE(XC.size() == 4);
      for (size_t c = 0; c < 4; c++)
      {
        Go::Point X;
        geo->point(X,ku[i+c%2],kv[j+c/2]);
        CHECK_THAT(XC[c].x, WithinAbs(X[0], 1e-14));
        CHECK_THAT(XC[c].y, WithinAbs(X[1], 1e-14));
      }
    }
}
