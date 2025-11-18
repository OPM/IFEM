//==============================================================================
//!
//! \file SplineFieldTests.h
//!
//! \date Mar 9 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for spline fields.
//!
//==============================================================================

#include "ASMenums.h"
#include "ASMmxBase.h"
#include "Field.h"
#include "Fields.h"
#include "ItgPoint.h"

#include "Catch2Support.h"
#include <catch2/generators/catch_generators_range.hpp>

#include <array>
#include <memory>
#include <vector>

namespace {

  class ItgPointLocInit : public ItgPoint
  {
  public:
    //! \brief Constructor initializing the domain parameters and element number.
    //! \details The local parameters are derived assuming a and b live in the 0..1 domain.
    ItgPointLocInit(double a, double b)
    {
      iGP = 0;
      idx = 0;
      u = a;
      v = b;
      xi = -1.0 + 2.0 * a;
      eta = -1.0 + 2.0 * b;
      iel = 1;
      w = zeta = 0.0;
    }

    //! \brief Constructor initializing the domain parameters and element number.
    //! \details The local parameters are derived assuming a and b live in the 0..1 domain.
    ItgPointLocInit(double a, double b, double c)
    {
      iGP = 0;
      idx = 0;
      u = a;
      v = b;
      w = c;
      xi = -1.0 + 2.0 * a;
      eta = -1.0 + 2.0 * b;
      zeta = -1.0 + 2.0 * c;
      iel = 1;
    }

  };

  const auto test_points2D = std::array{
    ItgPointLocInit{0.5, 0.5},
    ItgPointLocInit{1.0, 0.0},
    ItgPointLocInit{0.0, 1.0},
    ItgPointLocInit{1.0, 1.0},
  };

  const auto linear2D = std::vector{
    0.0, 1.0,
    1.0, 2.0
  }; // x + y

  const auto linear2DV = {
    0.0,  0.0,
    1.0,  1.0,
    1.0, -1.0,
    3.0, 1.0
  }; // {x+y+x*y, x-y+x*y}

  const auto linear2DVmx = std::vector{
    // x
    0.0, 0.5, 1.0,
    1.0, 2.0, 3.0,
    // y
    0.0,  1.0, -0.5,
    1.0, -1.0,  1.0,
    // p
    0.0, 1.0,
    1.0, 2.0
  }; // {x+y+x*y, x-y+x*y, x + y}

  const auto test_points3D = std::array{
    ItgPointLocInit{0.5, 0.5, 0.5},
    ItgPointLocInit{0.0, 0.0, 0.0},
    ItgPointLocInit{1.0, 0.0, 0.0},
    ItgPointLocInit{0.0, 1.0, 0.0},
    ItgPointLocInit{1.0, 1.0, 0.0},
    ItgPointLocInit{0.0, 0.0, 1.0},
    ItgPointLocInit{1.0, 0.0, 1.0},
    ItgPointLocInit{0.0, 1.0, 1.0},
    ItgPointLocInit{1.0, 1.0, 1.0},
  };

  const auto linear3D = std::vector{
    0.0, 1.0, 1.0, 2.0,
    1.0, 2.0, 2.0, 3.0
  }; // x + y + z

  const auto linear3DV = std::vector{
    0.0,  0.0,  0.0,
    1.0,  1.0,  1.0,
    1.0,  1.0, -1.0,
    2.0,  2.0,  0.0,
    1.0, -1.0,  1.0,
    2.0,  0.0,  2.0,
    2.0,  0.0,  0.0,
    4.0,  2.0,  2.0
  }; // {x+y+z+x*y*z, x+y-z+x*y*z, x-y+z+x*y*z}

  const auto linear3DVmx = std::vector{
    // x
    0.0, 0.5, 1.0,
    1.0, 1.5, 2.0,
    1.0, 1.5, 2.0,
    2.0, 3.0, 4.0,
    // y
    0.0, 1.0,
    0.5, 1.5,
    1.0, 2.0,
    -1.0, 0.0,
    -0.5, 1.0,
    0.0, 2.0,
    // z
     0.0, 1.0,
    -1.0, 0.0,
     0.5, 1.5,
    -0.5, 1.0,
     1.0, 2.0,
     0.0, 2.0,
    // p
    0.0, 1.0,
    1.0, 2.0,
    1.0, 2.0,
    2.0, 3.0
  }; // {x+y+z+x*y*z, x+y-z+x*y*z, x-y+z+x*y*z, x+y+z}
}


template<class Patch>
struct Field2DTests
{
  static void Value()
  {
    Patch patch(1);
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Field> flinear(Field::create(&patch, linear2D));
    REQUIRE(flinear);

    auto v = [](double x, double y) { return x + y; };

    const auto param = GENERATE(Catch::Generators::from_range(test_points2D));
    SECTION("(" + std::to_string(param.u) + "," + std::to_string(param.v) + ")")
    {
      REQUIRE_THAT(flinear->valueFE(param), WithinRel(v(param.u, param.v)));
    }
  }

  static void ValueQuad()
  {
    Patch patch(1);
    REQUIRE(patch.raiseOrder(1, 1));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Field> fquad(
      Field::create(&patch,
                    std::vector(Patch::quad.begin(), Patch::quad.end()))
    );
    REQUIRE(fquad);

    auto v = [](double x, double y) { return x*x + y*y + x*x*y; };

    const auto param = GENERATE(Catch::Generators::from_range(test_points2D));
    SECTION("(" + std::to_string(param.u) + "," + std::to_string(param.v) + ")")
    {
      REQUIRE_THAT(fquad->valueFE(param), WithinRel(v(param.u, param.v)));
    }
  }

  static void Grad()
  {
    Patch patch(1);
    REQUIRE(patch.raiseOrder(1,1));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Field> fscalar(
      Field::create(&patch,
                    std::vector(Patch::quad.begin(), Patch::quad.end()))
    );
    REQUIRE(fscalar);

    GradChecks(*fscalar);
  }

  static void GradSepGeom()
  {
    Patch patch(1);
    REQUIRE(patch.createProjectionBasis(true));
    REQUIRE(patch.raiseOrder(1, 1));
    REQUIRE(patch.createProjectionBasis(false));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Field> fscalar(
      Field::create(&patch,
                    std::vector(Patch::quad.begin(), Patch::quad.end()),
                    ASM::PROJECTION_BASIS)
    );
    REQUIRE(fscalar);

    GradChecks(*fscalar);
  }

  static void Hessian()
  {
    Patch patch(1);
    REQUIRE(patch.raiseOrder(1,1));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Field> fscalar(
      Field::create(&patch,
                    std::vector(Patch::quad.begin(), Patch::quad.end()))
    );
    REQUIRE(fscalar);

    HessianChecks(*fscalar);
  }

  static void HessianSepGeom()
  {
    Patch patch(1);
    REQUIRE(patch.createProjectionBasis(true));
    REQUIRE(patch.raiseOrder(1,1));
    REQUIRE(patch.createProjectionBasis(false));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Field> fscalar(
      Field::create(&patch,
                    std::vector(Patch::quad.begin(), Patch::quad.end()),
                    ASM::PROJECTION_BASIS)
    );
    REQUIRE(fscalar);

    HessianChecks(*fscalar);
  }

private:
  static void GradChecks(Field& fscalar)
  {
    auto dx = [](double x, double y) { return 2.0*x + 2.0*x*y; };
    auto dy = [](double x, double y) { return 2.0*y + x*x; };

    const auto param = GENERATE(Catch::Generators::from_range(test_points2D));
    SECTION("(" + std::to_string(param.u) + "," + std::to_string(param.v) + ")")
    {
      Vector v(2);
      REQUIRE(fscalar.gradFE(param, v));
      REQUIRE_THAT(v(1), WithinRel(dx(param.u, param.v)));
      REQUIRE_THAT(v(2), WithinRel(dy(param.u, param.v)));
    }
  }

  static void HessianChecks(Field& fscalar)
  {
    auto hxx = [](double x, double y) { return 2.0 + 2.0*y; };
    auto hxy = [](double x, double y) { return 2.0*x; };

    const auto param = GENERATE(Catch::Generators::from_range(test_points2D));
    SECTION("(" + std::to_string(param.u) + "," + std::to_string(param.v) + ")")
    {
      REQUIRE_THAT(fscalar.valueFE(param),
                    WithinRel(param.u*param.u + param.v*param.v +
                              param.u*param.u*param.v));
      Matrix H(2,2);
      REQUIRE(fscalar.hessianFE(param, H));
      REQUIRE_THAT(H(1,1), WithinRel(hxx(param.u, param.v)));
      REQUIRE_THAT(H(2,2), WithinRel(2.0));
      REQUIRE(H(1,2) == H(2,1));
      REQUIRE_THAT(H(1,2), WithinRel(hxy(param.u, param.v)));
    }
  }
};


template<class Patch>
struct Fields2DTests
{
  static void Value()
  {
    Patch patch;
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Fields> fvector(Fields::create(&patch,linear2DV));
    REQUIRE(fvector);
    std::unique_ptr<Field> fscalar(Field::create(&patch,linear2DV,1,2));
    REQUIRE(fscalar);

    auto vx = [](double x, double y) { return x + y + x*y; };
    auto vy = [](double x, double y) { return x - y + x*y; };

    const auto param = GENERATE(Catch::Generators::from_range(test_points2D));
    SECTION("(" + std::to_string(param.u) + "," + std::to_string(param.v) + ")")
    {
      Vector v(2);
      REQUIRE(fvector->valueFE(param,v));
      REQUIRE_THAT(v(1), WithinRel(vx(param.u, param.v)));
      REQUIRE_THAT(v(2), WithinRel(vy(param.u, param.v)));
      REQUIRE_THAT(fscalar->valueFE(param), WithinRel(vy(param.u, param.v)));
    }
  }

  static void ValueQuad()
  {
    Patch patch;
    REQUIRE(patch.raiseOrder(1,1));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Fields> fvector(
      Fields::create(&patch,
                     std::vector(Patch::quadv.begin(), Patch::quadv.end()))
    );
    REQUIRE(fvector);
    std::unique_ptr<Field> fscalar(
      Field::create(&patch,
                     std::vector(Patch::quadv.begin(), Patch::quadv.end()),
                     1,2)
    );
    REQUIRE(fscalar);

    auto vx = [](double x, double y) { return x*x + y*y + x*y*y; };
    auto vy = [](double x, double y) { return x*x + y*y + x*x*y; };

    const auto param = GENERATE(Catch::Generators::from_range(test_points2D));
    SECTION("(" + std::to_string(param.u) + "," + std::to_string(param.v) + ")") {
      Vector v(2);
      REQUIRE(fvector->valueFE(param,v));
      REQUIRE_THAT(v(1), WithinRel(vx(param.u, param.v)));
      REQUIRE_THAT(v(2), WithinRel(vy(param.u,param.v)));
      REQUIRE_THAT(fscalar->valueFE(param), WithinRel(vy(param.u, param.v)));
    }
  }

  static void Valuemx()
  {
    ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
    Patch patch({1,1,1});
    REQUIRE(patch.raiseOrder(1,1));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Fields> fvector(Fields::create(&patch,linear2DVmx,12));
    REQUIRE(fvector);
    std::unique_ptr<Field> fscalar(Field::create(&patch,linear2DVmx,3,1));
    REQUIRE(fscalar);
    const auto param = GENERATE(Catch::Generators::from_range(test_points2D));

    auto vx = [](double x, double y) { return x + y + x*y; };
    auto vy = [](double x, double y) { return x - y + x*y; };
    auto vp = [](double x, double y) { return x + y; };

    SECTION("(" + std::to_string(param.u) + "," + std::to_string(param.v) + ")")
    {
      Vector v(2);
      REQUIRE(fvector->valueFE(param,v));
      REQUIRE_THAT(v(1), WithinRel(vx(param.u, param.v)));
      REQUIRE_THAT(v(2), WithinRel(vy(param.u, param.v)));
      REQUIRE_THAT(fscalar->valueFE(param), WithinRel(vp(param.u, param.v)));
    }
  }

  static void Grad()
  {
    Patch patch;
    REQUIRE(patch.raiseOrder(1,1));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Fields> fvector(
      Fields::create(&patch,
                     std::vector(Patch::quadv.begin(), Patch::quadv.end()))
    );
    REQUIRE(fvector);

    GradChecks(*fvector);
  }

  static void Gradmx()
  {
    ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
    Patch patch({1,1,1});
    REQUIRE(patch.raiseOrder(2,2));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Fields> fvector(
      Fields::create(&patch,
                     std::vector(Patch::quadv.begin(), Patch::quadv.end()),
                     12)
    );
    REQUIRE(fvector);
    std::unique_ptr<Field> fscalar(
      Field::create(&patch,
                    std::vector(Patch::quadv.begin(), Patch::quadv.end()),
                    3, 1)
    );

    GradChecks(*fvector);

    auto dpdx = [](double x, double y) { return y + 2.0*x; };
    auto dpdy = [](double x, double y) { return x - 2.0*y; };

    const auto param = GENERATE(Catch::Generators::from_range(test_points2D));
    SECTION("Scalar (" + std::to_string(param.u) + "," + std::to_string(param.v) + ")")
    {
      Vector gradp(2);
      REQUIRE(fscalar->gradFE(param,gradp));
      REQUIRE_THAT(gradp(1), WithinRel(dpdx(param.u, param.v)));
      REQUIRE_THAT(gradp(2), WithinRel(dpdy(param.u, param.v)));
    }
  }

  static void GradSepGeom()
  {
    Patch patch;
    REQUIRE(patch.createProjectionBasis(true));
    REQUIRE(patch.raiseOrder(1,1));
    REQUIRE(patch.createProjectionBasis(false));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Fields> fvector(
      Fields::create(&patch,
                     std::vector(Patch::quadv.begin(), Patch::quadv.end()),
                     ASM::PROJECTION_BASIS)
    );
    REQUIRE(fvector);

    GradChecks(*fvector);
  }

  static void Hessian()
  {
    Patch patch;
    REQUIRE(patch.raiseOrder(1,1));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Fields> fvector(
      Fields::create(&patch,
                     std::vector(Patch::quadv.begin(), Patch::quadv.end()))
    );
    REQUIRE(fvector);

    HessianChecks(*fvector);
  }

  static void Hessianmx()
  {
    ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
    Patch patch({1,1,1});
    REQUIRE(patch.raiseOrder(2,2));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Fields> fvector(
      Fields::create(&patch,
                     std::vector(Patch::quadv.begin(), Patch::quadv.end()),
                     12)
    );
    REQUIRE(fvector);
    std::unique_ptr<Field> fscalar(
      Field::create(&patch,
                    std::vector(Patch::quadv.begin(), Patch::quadv.end()),
                    3, 1)
    );
    REQUIRE(fscalar);

    HessianChecks(*fvector);

    const auto param = GENERATE(Catch::Generators::from_range(test_points2D));
    SECTION("Scalar (" + std::to_string(param.u) + "," + std::to_string(param.v) + ")")
    {
      Matrix Hs(2,2);
      REQUIRE(fscalar->hessianFE(param, Hs));
      REQUIRE_THAT(Hs(1,1), WithinRel(2.0));
      REQUIRE(Hs(1,2) == Hs(2,1));
      REQUIRE_THAT(Hs(1,2), WithinRel(1.0));
      REQUIRE_THAT(Hs(2,2), WithinRel(-2.0));
    }
  }

  static void HessianSepGeom()
  {
    Patch patch;
    REQUIRE(patch.createProjectionBasis(true));
    REQUIRE(patch.raiseOrder(1,1));
    REQUIRE(patch.createProjectionBasis(false));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Fields> fvector(
      Fields::create(&patch,
                     std::vector(Patch::quadv.begin(), Patch::quadv.end()),
                     ASM::PROJECTION_BASIS)
    );
    REQUIRE(fvector);

    HessianChecks(*fvector);
  }

private:
  static void GradChecks(Fields& fvector)
  {
    auto xdx = [](double x, double y) { return 2.0*x + 2.0*x*y; };
    auto xdy = [](double x, double y) { return 2.0*y + x*x; };
    auto ydx = [](double x, double y) { return 2.0*x + y*y; };
    auto ydy = [](double x, double y) { return 2.0*y + 2.0*x*y; };

    const auto param = GENERATE(Catch::Generators::from_range(test_points2D));
    SECTION("(" + std::to_string(param.u) + "," + std::to_string(param.v) + ")")
    {
      Matrix gradu(2,2);
      REQUIRE(fvector.gradFE(param,gradu));
      REQUIRE_THAT(gradu(1,1), WithinRel(xdx(param.u, param.v)));
      REQUIRE_THAT(gradu(1,2), WithinRel(xdy(param.u, param.v)));
      REQUIRE_THAT(gradu(2,1), WithinRel(ydx(param.u, param.v)));
      REQUIRE_THAT(gradu(2,2), WithinRel(ydy(param.u, param.v)));
    }
  }

  static void HessianChecks(Fields& fvector)
  {
    auto vx = [](double x, double y) { return x*x + y*y + x*x*y; };
    auto vy = [](double x, double y) { return x*x + y*y + x*y*y; };

    auto xdxy = [](double x, double y) { return 2.0*x; };
    auto xdxx = [](double x, double y) { return 2.0*(y + 1.0); };

    auto ydxy = [](double x, double y) { return 2.0*y; };
    auto ydyy = [](double x, double y) { return 2.0*(x + 1.0); };

    const auto param = GENERATE(Catch::Generators::from_range(test_points2D));
    SECTION("(" + std::to_string(param.u) + "," + std::to_string(param.v) + ")")
    {
      Vector v(2);
      REQUIRE(fvector.valueFE(param, v));
      REQUIRE_THAT(v(1), WithinRel(vx(param.u, param.v)));
      REQUIRE_THAT(v(2), WithinRel(vy(param.u, param.v)));

      Matrix3D H(2,2,2);
      REQUIRE(fvector.hessianFE(param, H));
      REQUIRE_THAT(H(1,1,1), WithinRel(xdxx(param.u, param.v)));
      REQUIRE_THAT(H(1,2,2), WithinRel(2.0));
      REQUIRE(H(1,1,2) == H(1,2,1));
      REQUIRE_THAT(H(1,1,2), WithinRel(xdxy(param.u, param.v)));

      REQUIRE_THAT(H(2,1,1), WithinRel(2.0));
      REQUIRE_THAT(H(2,2,2), WithinRel(ydyy(param.u, param.v)));
      REQUIRE(H(2,1,2) == H(2,2,1));
      REQUIRE_THAT(H(2,1,2), WithinRel(ydxy(param.u, param.v)));
    }
  }
};


template<class Patch>
struct Field3DTests
{
  static void Value()
  {
    Patch patch(1);
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Field> fscalar(Field::create(&patch, linear3D));
    REQUIRE(fscalar);

    auto v = [](double x, double y, double z) { return x + y + z; };

    const auto param = GENERATE(Catch::Generators::from_range(test_points3D));
    SECTION("(" + std::to_string(param.u) + ","
                + std::to_string(param.v) + ","
                + std::to_string(param.w)  + ")")
    {
      REQUIRE_THAT(fscalar->valueFE(param),
                    WithinRel(v(param.u, param.v, param.w)));
    }
  }

  static void ValueQuad()
  {
    Patch patch(1);
    REQUIRE(patch.raiseOrder(1,1,1,false));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Field> fscalar(
      Field::create(&patch,
                    std::vector(Patch::quad.begin(), Patch::quad.end()))
    );
    REQUIRE(fscalar);

    auto v = [](double x, double y, double z)
    { return x*x + y*y + z*z + x*x*y*z + x*y*y*z + x*y*z*z; };

    const auto param = GENERATE(Catch::Generators::from_range(test_points3D));
    SECTION("(" + std::to_string(param.u) + ","
                + std::to_string(param.v) + ","
                + std::to_string(param.w)  + ")")
    {
      REQUIRE_THAT(fscalar->valueFE(param),
                   WithinRel(v(param.u, param.v, param.w)));
    }
  }

  static void Grad()
  {
    Patch patch(1);
    REQUIRE(patch.raiseOrder(1,1,1,false));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Field> fscalar(
      Field::create(&patch,
                    std::vector(Patch::quad.begin(), Patch::quad.end()))
    );
    REQUIRE(fscalar);

    GradChecks(*fscalar);
  }

  static void GradSepGeom()
  {
    Patch patch(1);
    REQUIRE(patch.createProjectionBasis(true));
    REQUIRE(patch.raiseOrder(1, 1, 1, false));
    REQUIRE(patch.createProjectionBasis(false));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Field> fscalar(
      Field::create(&patch,
                    std::vector(Patch::quad.begin(), Patch::quad.end()),
                    ASM::PROJECTION_BASIS)
    );
    REQUIRE(fscalar);

    GradChecks(*fscalar);
  }

  static void Hessian()
  {
    Patch patch(1);
    REQUIRE(patch.raiseOrder(1, 1, 1, false));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Field> fscalar(
      Field::create(&patch,
                    std::vector(Patch::quad.begin(), Patch::quad.end()))
    );
    REQUIRE(fscalar);

    HessianChecks(*fscalar);
  }

  static void HessianSepGeom()
  {
    Patch patch(1);
    REQUIRE(patch.createProjectionBasis(true));
    REQUIRE(patch.raiseOrder(1, 1, 1, false));
    REQUIRE(patch.createProjectionBasis(false));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Field> fscalar(
      Field::create(&patch,
                    std::vector(Patch::quad.begin(), Patch::quad.end()),
                    ASM::PROJECTION_BASIS)
    );
    REQUIRE(fscalar);

    HessianChecks(*fscalar);
  }

private:
  static void GradChecks(Field& fscalar)
  {
    auto dx = [](double x, double y, double z)
    { return 2.0*x + 2.0*x*y*z + y*y*z + y*z*z; };

    auto dy = [](double x, double y, double z)
    { return 2.0*y + x*x*z + 2.0*x*y*z + x*z*z; };

    auto dz = [](double x, double y, double z)
    { return 2.0*z + x*x*y + x*y*y + 2.0*x*y*z; };

    const auto param = GENERATE(Catch::Generators::from_range(test_points3D));
    SECTION("(" + std::to_string(param.u) + ","
                + std::to_string(param.v) + ","
                + std::to_string(param.w)  + ")")
    {
      Vector v(3);
      REQUIRE(fscalar.gradFE(param, v));
      REQUIRE_THAT(v(1), WithinRel(dx(param.u, param.v, param.w)));
      REQUIRE_THAT(v(2), WithinRel(dy(param.u, param.v, param.w)));
      REQUIRE_THAT(v(3), WithinRel(dz(param.u, param.v, param.w)));
    }
  }

  static void HessianChecks(Field& fscalar)
  {
    auto dxx = [](double x, double y, double z) { return 2.0*(y*z + 1.0); };
    auto dyy = [](double x, double y, double z) { return 2.0*(x*z + 1.0); };
    auto dzz = [](double x, double y, double z) { return 2.0*(x*y + 1.0); };
    auto dxy = [](double x, double y, double z) { return z*(2.0*(x + y) + z); };
    auto dxz = [](double x, double y, double z) { return y*(2.0*(x + z) + y); };
    auto dyz = [](double x, double y, double z) { return x*(2.0*(y + z) + x); };

    const auto param = GENERATE(Catch::Generators::from_range(test_points3D));
    SECTION("(" + std::to_string(param.u) + ","
                + std::to_string(param.v) + ","
                + std::to_string(param.w)  + ")")
    {
      Matrix H(3,3);
      REQUIRE(fscalar.hessianFE(param, H));
      REQUIRE_THAT(H(1,1), WithinRel(dxx(param.u, param.v, param.w)));
      REQUIRE_THAT(H(2,2), WithinRel(dyy(param.u, param.v, param.w)));
      REQUIRE_THAT(H(3,3), WithinRel(dzz(param.u, param.v, param.w)));
      REQUIRE(H(1,2) == H(2,1));
      REQUIRE(H(1,3) == H(3,1));
      REQUIRE(H(2,3) == H(3,2));
      REQUIRE_THAT(H(1,2), WithinRel(dxy(param.u, param.v, param.w)));
      REQUIRE_THAT(H(1,3), WithinRel(dxz(param.u, param.v, param.w)));
      REQUIRE_THAT(H(2,3), WithinRel(dyz(param.u, param.v, param.w)));
    }
  }
};


template<class Patch>
struct Fields3DTests
{
  static void Value()
  {
    Patch patch;
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Fields> fvector(Fields::create(&patch,linear3DV));
    REQUIRE(fvector);
    std::unique_ptr<Field> fscalar(Field::create(&patch,linear3DV,1,2));
    REQUIRE(fscalar);

    auto vx = [](double x, double y, double z) { return x + y + z + x*y*z; };
    auto vy = [](double x, double y, double z) { return x + y - z + x*y*z; };
    auto vz = [](double x, double y, double z) { return x - y + z + x*y*z; };

    const auto param = GENERATE(Catch::Generators::from_range(test_points3D));
    SECTION("(" + std::to_string(param.u) + ","
                + std::to_string(param.v) + ","
                + std::to_string(param.w)  + ")")
    {
      Vector v(3);
      REQUIRE(fvector->valueFE(param,v));
      REQUIRE_THAT(v(1), WithinRel(vx(param.u, param.v, param.w)));
      REQUIRE_THAT(v(2), WithinRel(vy(param.u, param.v, param.w)));
      REQUIRE_THAT(v(3), WithinRel(vz(param.u, param.v, param.w)));
      REQUIRE_THAT(fscalar->valueFE(param),
                   WithinRel(vy(param.u, param.v, param.w)));
    }
  }

  static void ValueQuad()
  {
    Patch patch;
    REQUIRE(patch.raiseOrder(1,1,1,false));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Fields> fvector(
      Fields::create(&patch,
                     std::vector(Patch::quadv.begin(), Patch::quadv.end()))
    );
    REQUIRE(fvector);
    std::unique_ptr<Field> fscalar(
      Field::create(&patch,
                    std::vector(Patch::quadv.begin(), Patch::quadv.end()),1,2)
    );
    REQUIRE(fscalar);

    auto vx = [](double x, double y, double z)
    { return x*x + y*y + z*z + x*x*y*z; };
    auto vy = [](double x, double y, double z)
    { return x*x + y*y + z*z + x*y*y*z; };
    auto vz = [](double x, double y, double z)
    { return x*x + y*y + z*z + x*y*z*z; };

    const auto param = GENERATE(Catch::Generators::from_range(test_points3D));
    SECTION("Quad (" + std::to_string(param.u) + ","
                     + std::to_string(param.v) + ","
                     + std::to_string(param.w)  + ")")
    {
      Vector v(3);
      REQUIRE(fvector->valueFE(param,v));
      REQUIRE_THAT(v(1), WithinRel(vx(param.u, param.v, param.w)));
      REQUIRE_THAT(v(2), WithinRel(vy(param.u, param.v, param.w)));
      REQUIRE_THAT(v(3), WithinRel(vz(param.u, param.v, param.w)));
      REQUIRE_THAT(fscalar->valueFE(param),
                   WithinRel(vy(param.u, param.v, param.w)));
    }
  }

  static void Valuemx()
  {
    ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
    Patch patch({1,1,1,1});
    REQUIRE(patch.raiseOrder(1,1,1,false));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Fields> fvector(Fields::create(&patch, linear3DVmx, 123));
    REQUIRE(fvector);
    std::unique_ptr<Field> fscalar(Field::create(&patch, linear3DVmx, 4, 1));
    REQUIRE(fscalar);

    auto vx = [](double x, double y, double z) { return x + y + z + x*y*z; };
    auto vy = [](double x, double y, double z) { return x + y - z + x*y*z; };
    auto vz = [](double x, double y, double z) { return x - y + z + x*y*z; };
    auto vp = [](double x, double y, double z) { return x + y + z; };

    const auto param = GENERATE(Catch::Generators::from_range(test_points3D));
    SECTION("(" + std::to_string(param.u) + ","
                + std::to_string(param.v) + ","
                + std::to_string(param.w)  + ")")
    {
      Vector v(3);
      REQUIRE(fvector->valueFE(param, v));
      REQUIRE_THAT(v(1), WithinRel(vx(param.u, param.v, param.w)));
      REQUIRE_THAT(v(2), WithinRel(vy(param.u, param.v, param.w)));
      REQUIRE_THAT(v(3), WithinRel(vz(param.u, param.v, param.w)));
      REQUIRE_THAT(fscalar->valueFE(param),
                   WithinRel(vp(param.u, param.v, param.w)));
    }
  }

  static void Grad()
  {
    Patch patch;
    REQUIRE(patch.raiseOrder(1,1,1,false));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Fields> fvector(
      Fields::create(&patch,
                     std::vector(Patch::quadv.begin(), Patch::quadv.end()))
    );
    REQUIRE(fvector);

    GradChecks(*fvector);
  }

  static void Gradmx()
  {
    ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
    Patch patch({1,1,1,1});
    REQUIRE(patch.raiseOrder(2,2,2,false));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Fields> fvector(
      Fields::create(&patch,
                     std::vector(Patch::quadv.begin(), Patch::quadv.end()),
                     123)
    );
    REQUIRE(fvector);
    std::unique_ptr<Field> fscalar(
      Field::create(&patch,
                    std::vector(Patch::quadv.begin(), Patch::quadv.end()),4,1)
    );
    REQUIRE(fscalar);

    GradChecks(*fvector);

    auto dpdx =  [](double x, double y, double z) { return y*z + 2.0*x; };
    auto dpdy =  [](double x, double y, double z) { return x*z - 2.0*y; };
    auto dpdz =  [](double x, double y, double z) { return x*y + 2.0*z; };

    const auto param = GENERATE(Catch::Generators::from_range(test_points3D));
    SECTION("Scalar (" + std::to_string(param.u) + ","
                       + std::to_string(param.v) + ","
                       + std::to_string(param.w)  + ")")
    {
      Vector gradp(3);
      REQUIRE(fscalar->gradFE(param,gradp));
      REQUIRE_THAT(gradp(1), WithinRel(dpdx(param.u, param.v, param.w)));
      REQUIRE_THAT(gradp(2), WithinRel(dpdy(param.u, param.v, param.w)));
      REQUIRE_THAT(gradp(3), WithinRel(dpdz(param.u, param.v, param.w)));
    }
  }

  static void GradSepGeom()
  {
    Patch patch;
    REQUIRE(patch.createProjectionBasis(true));
    REQUIRE(patch.raiseOrder(1,1,1,false));
    REQUIRE(patch.createProjectionBasis(false));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Fields> fvector(
      Fields::create(&patch,
                     std::vector(Patch::quadv.begin(), Patch::quadv.end()),
                     ASM::PROJECTION_BASIS)
    );
    REQUIRE(fvector);

    GradChecks(*fvector);
  }

  static void Hessian()
  {
    Patch patch;
    REQUIRE(patch.raiseOrder(1,1,1,false));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Fields> fvector(
      Fields::create(&patch,
                     std::vector(Patch::quadv.begin(), Patch::quadv.end()))
    );
    REQUIRE(fvector);

    HessianChecks(*fvector);
  }

  static void Hessianmx()
  {
    ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
    Patch patch({1,1,1,1});
    REQUIRE(patch.raiseOrder(2,2,2,false));
    REQUIRE(patch.generateFEMTopology());

    std::unique_ptr<Fields> fvector(
      Fields::create(&patch,
                     std::vector(Patch::quadv.begin(), Patch::quadv.end()),
                     123)
    );
    REQUIRE(fvector);
    std::unique_ptr<Field> fscalar(
      Field::create(&patch,
                    std::vector(Patch::quadv.begin(), Patch::quadv.end()),
                    4, 1)
    );
    REQUIRE(fscalar);

    HessianChecks(*fvector);

    const auto param = GENERATE(Catch::Generators::from_range(test_points3D));
    SECTION("Scalar (" + std::to_string(param.u) + ","
                + std::to_string(param.v) + ","
                + std::to_string(param.w)  + ")")
    {
      Matrix Hs(3,3);
      REQUIRE(fscalar->hessianFE(param, Hs));
      REQUIRE_THAT(Hs(1,1), WithinRel(2.0));
      REQUIRE_THAT(Hs(2,2), WithinRel(-2.0));
      REQUIRE_THAT(Hs(3,3), WithinRel(2.0));
      REQUIRE(Hs(1,2) == Hs(2,1));
      REQUIRE(Hs(1,3) == Hs(3,1));
      REQUIRE(Hs(2,3) == Hs(3,2));
      REQUIRE_THAT(Hs(1,2), WithinAbs(param.w, 1e-12));
      REQUIRE_THAT(Hs(1,3), WithinAbs(param.v, 1e-12));
      REQUIRE_THAT(Hs(2,3), WithinAbs(param.u, 1e-12));
    }
  }

  static void HessianSepGeom()
  {
    Patch patch;
    REQUIRE(patch.createProjectionBasis(true));
    REQUIRE(patch.raiseOrder(1,1,1,false));
    REQUIRE(patch.createProjectionBasis(false));
    REQUIRE(patch.generateFEMTopology());
    std::unique_ptr<Fields> fvector(
      Fields::create(&patch,
                     std::vector(Patch::quadv.begin(), Patch::quadv.end()),
                     ASM::PROJECTION_BASIS)
    );
    REQUIRE(fvector);

    HessianChecks(*fvector);
  }

private:
  static void GradChecks(Fields& fvector)
  {
    auto xdx = [](double x, double y, double z) { return 2.0*(x + x*y*z); };
    auto xdy = [](double x, double y, double z) { return 2.0*y + x*x*z; };
    auto xdz = [](double x, double y, double z) { return 2.0*z + x*x*y; };
    auto ydx = [](double x, double y, double z) { return 2.0*x + y*y*z; };
    auto ydy = [](double x, double y, double z) { return 2.0*(y + x*y*z); };
    auto ydz = [](double x, double y, double z) { return 2.0*z + x*y*y; };
    auto zdx = [](double x, double y, double z) { return 2.0*x + y*z*z; };
    auto zdy = [](double x, double y, double z) { return 2.0*y + x*z*z; };
    auto zdz = [](double x, double y, double z) { return 2.0*(z + x*y*z); };

    const auto param = GENERATE(Catch::Generators::from_range(test_points3D));
    SECTION("(" + std::to_string(param.u) + ","
                + std::to_string(param.v) + ","
                + std::to_string(param.w)  + ")")
    {
      Matrix gradu(3,3);
      REQUIRE(fvector.gradFE(param,gradu));
      REQUIRE_THAT(gradu(1,1), WithinRel(xdx(param.u, param.v, param.w)));
      REQUIRE_THAT(gradu(1,2), WithinRel(xdy(param.u, param.v, param.w)));
      REQUIRE_THAT(gradu(1,3), WithinRel(xdz(param.u, param.v, param.w)));
      REQUIRE_THAT(gradu(2,1), WithinRel(ydx(param.u, param.v, param.w)));
      REQUIRE_THAT(gradu(2,2), WithinRel(ydy(param.u, param.v, param.w)));
      REQUIRE_THAT(gradu(2,3), WithinRel(ydz(param.u, param.v, param.w)));
      REQUIRE_THAT(gradu(3,1), WithinRel(zdx(param.u, param.v, param.w)));
      REQUIRE_THAT(gradu(3,2), WithinRel(zdy(param.u, param.v, param.w)));
      REQUIRE_THAT(gradu(3,3), WithinRel(zdz(param.u, param.v, param.w)));
    }
  }

  static void HessianChecks(Fields& fvector)
  {
    auto xdxx = [](double x, double y, double z) { return 2.0*(y*z + 1.0); };
    auto xdxy = [](double x, double y, double z) { return 2.0*x*z; };
    auto xdxz = [](double x, double y, double z) { return 2.0*x*y; };
    auto xdyz = [](double x, double y, double z) { return x*x; };

    auto ydyy = [](double x, double y, double z) { return 2.0*(x*z + 1.0); };
    auto ydxy = [](double x, double y, double z) { return 2.0*y*z; };
    auto ydxz = [](double x, double y, double z) { return y*y; };
    auto ydyz = [](double x, double y, double z) { return 2.0*x*y; };

    auto zdzz = [](double x, double y, double z) { return 2.0*(x*y + 1.0); };
    auto zdxy = [](double x, double y, double z) { return z*z; };
    auto zdxz = [](double x, double y, double z) { return 2.0*y*z; };
    auto zdyz = [](double x, double y, double z) { return 2.0*x*z; };

    const auto param = GENERATE(Catch::Generators::from_range(test_points3D));
    SECTION("(" + std::to_string(param.u) + ","
                + std::to_string(param.v) + ","
                + std::to_string(param.w)  + ")")
    {
      Matrix3D H(3,3,3);
      REQUIRE(fvector.hessianFE(param,H));
      REQUIRE_THAT(H(1,1,1), WithinRel(xdxx(param.u, param.v, param.w)));
      REQUIRE_THAT(H(1,2,2), WithinRel(2.0));
      REQUIRE_THAT(H(1,3,3), WithinRel(2.0));
      REQUIRE(H(1,1,2) == H(1,2,1));
      REQUIRE(H(1,1,3) == H(1,3,1));
      REQUIRE(H(1,2,3) == H(1,3,2));
      REQUIRE_THAT(H(1,1,2), WithinRel(xdxy(param.u, param.v, param.w)));
      REQUIRE_THAT(H(1,1,3), WithinRel(xdxz(param.u, param.v, param.w)));
      REQUIRE_THAT(H(1,2,3), WithinRel(xdyz(param.u, param.v, param.w)));

      REQUIRE_THAT(H(2,1,1), WithinRel(2.0));
      REQUIRE_THAT(H(2,2,2), WithinRel(ydyy(param.u, param.v, param.w)));
      REQUIRE_THAT(H(2,3,3), WithinRel(2.0));
      REQUIRE(H(2,1,2) == H(2,2,1));
      REQUIRE(H(2,1,3) == H(2,3,1));
      REQUIRE(H(2,2,3) == H(2,3,2));
      REQUIRE_THAT(H(2,1,2), WithinRel(ydxy(param.u, param.v, param.w)));
      REQUIRE_THAT(H(2,1,3), WithinRel(ydxz(param.u, param.v, param.w)));
      REQUIRE_THAT(H(2,2,3), WithinRel(ydyz(param.u, param.v, param.w)));

      REQUIRE_THAT(H(3,1,1), WithinRel(2.0));
      REQUIRE_THAT(H(3,2,2), WithinRel(2.0));
      REQUIRE_THAT(H(3,3,3), WithinRel(zdzz(param.u, param.v, param.w)));
      REQUIRE(H(3,1,2) == H(3,2,1));
      REQUIRE(H(3,1,3) == H(3,3,1));
      REQUIRE(H(3,2,3) == H(3,3,2));
      REQUIRE_THAT(H(3,1,2), WithinRel(zdxy(param.u, param.v, param.w)));
      REQUIRE_THAT(H(3,1,3), WithinRel(zdxz(param.u, param.v, param.w)));
      REQUIRE_THAT(H(3,2,3), WithinRel(zdyz(param.u, param.v, param.w)));
    }
  }
};
