//==============================================================================
//!
//! \file ASMCube.h
//!
//! \date Oct 2 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Default tri-unit cube patch for unit-testing.
//!
//==============================================================================

#ifndef _ASM_CUBE_H
#define _ASM_CUBE_H

#include "ASMs3Dmx.h"
#include "ASMs3DLag.h"

#include <array>
#include <sstream>

template<class Base>
class ASMCubeBase : public Base
{
public:
  // formatting matches write routine, do not change
  static constexpr const char* cube = "700 1 0 0\n"
                                      "3 0\n"
                                      "2 2\n"
                                      "0 0 1 1\n"
                                      "2 2\n"
                                      "0 0 1 1\n"
                                      "2 2\n"
                                      "0 0 1 1\n"
                                      "0 0 0\n"
                                      "1 0 0\n"
                                      "0 1 0\n"
                                      "1 1 0\n"
                                      "0 0 1\n"
                                      "1 0 1\n"
                                      "0 1 1\n"
                                      "1 1 1\n\n";

  static constexpr auto quad = std::array{
    0.0,  0.0,  1.0,
    0.0,  0.0,  1.0,
    1.0,  1.0,  2.0,
    0.0,  0.0,  1.0,
    0.0,  0.0, 1.25,
    1.0, 1.25,  3.0,
    1.0,  1.0,  2.0,
    1.0, 1.25,  3.0,
    2.0,  3.0,  6.0,
  }; // x^2 + y^2 + z^2 + x^2yz + xy^2z + xyz^2

  static constexpr auto quadv = std::array{
     0.0,  0.0, 0.0,
     0.0,  0.0, 0.0,
     1.0,  1.0, 1.0,
     0.0,  0.0, 0.0,
     0.0,  0.0, 0.0,
     1.0,  1.0, 1.0,
     1.0,  1.0, 1.0,
     1.0,  1.0, 1.0,
     2.0,  2.0, 2.0,
     0.0,  0.0, 0.0,
     0.0,  0.0, 0.0,
     1.0,  1.0, 1.0,
     0.0,  0.0, 0.0,
     0.0,  0.0, 0.0,
    1.25,  1.0, 1.0,
     1.0,  1.0, 1.0,
     1.0, 1.25, 1.0,
     2.5,  2.5, 2.0,
     1.0,  1.0, 1.0,
     1.0,  1.0, 1.0,
     2.0,  2.0, 2.0,
     1.0,  1.0, 1.0,
     1.0,  1.0, 1.25,
     2.5,  2.0, 2.5,
     2.0,  2.0, 2.0,
     2.0,  2.5, 2.5,
     4.0,  4.0, 4.0,
  }; // {x^2 + y^2 + z^2 + x^2yz,
     //  x^2 + y^2 + z^2 + xy^2z,
     //  x^2 + y^2 + z^2 + xyz^2}

  explicit ASMCubeBase(unsigned char n_f = 3) : Base(n_f)
  {
    std::stringstream geo(cube);
    this->read(geo);
  }

  void shiftElemNumbers(int shift)
  {
    for (int& e : this->myMLGE)
      e += (e == -1 ? 0 : shift);
  }
};

using ASMCube = ASMCubeBase<ASMs3D>;


class ASMmxCube : public ASMs3Dmx
{
public:
  static constexpr auto quadv = std::array{
    // x
    0.0, 0.0, 1.0 / 3.0,   1.0,
    0.0, 0.0, 1.0 / 3.0,   1.0,
    1.0, 1.0, 4.0 / 3.0,   2.0,
    0.0, 0.0, 1.0 / 3.0,   1.0,
    0.0, 0.0, 5.0 / 12.0, 1.25,
    1.0, 1.0, 1.5,         2.5,
    1.0, 1.0, 4.0 / 3.0,   2.0,
    1.0, 1.0, 1.5,         2.5,
    2.0, 2.0, 8.0 / 3.0,   4.0,
    // y
    0.0,       0.0,        1.0,
    0.0,       0.0,        1.0,
    1.0 / 3.0, 1.0 / 3.0,  4.0 / 3.0,
    1.0,       1.0,        2.0,
    0.0,       0.0,        1.0,
    0.0,       0.0,        1.0,
    1.0 / 3.0, 5.0 / 12.0, 1.5,
    1.0,       1.25,       2.5,
    1.0,       1.0,        2.0,
    1.0,       1.0,        2.0,
    4.0 / 3.0, 1.5,        8.0 / 3.0,
    2.0,       2.5,        4.0,
    // z
    0.0,       0.0,        1.0,
    0.0,       0.0,        1.0,
    1.0,       1.0,        2.0,
    0.0,       0.0,        1.0,
    0.0,       0.0,        1.0,
    1.0,       1.0,        2.0,
    1.0 / 3.0, 1.0 / 3.0,  4.0 / 3.0,
    1.0 / 3.0, 5.0 / 12.0, 1.5,
    4.0 / 3.0, 1.5,        8.0 / 3.0,
    1.0,       1.0,        2.0,
    1.0,       1.25,       2.5,
    2.0,       2.5,        4.0,
    // p
     0.0,   0.0,   1.0,
     0.0,   0.0,   1.0,
    -1.0,  -1.0,   0.0,
     0.0,   0.0,   1.0,
     0.0,   0.125, 1.25,
    -1.0,  -0.75,  0.5,
    1.0,    1.0,   2.0,
    1.0,    1.25,  2.5,
    0.0,    0.5,   2.0,
  }; // {x^2 + y^2 + z^2 + x^2yz,
     //  x^2 + y^2 +xy^2z,
     //  x^2 + y^2 + xyz^2,
     //  x^2 - y^2 + z^2 + xyz}

  explicit ASMmxCube(const CharVec& n_f) : ASMs3Dmx(n_f)
  {
    std::stringstream geo(ASMCube::cube);
    this->read(geo);
  }
  virtual ~ASMmxCube() {}
};


namespace {
  constexpr double quad1(double x, double y, double z)
  { return x*x + y*y + z*z + x*x*y*z + x*y*y*z + x*y*z*z; }
  constexpr double quadv1(double x, double y, double z)
  { return x*x + y*y + z*z + x*x*y*z; }
  constexpr double quadv2(double x, double y, double z)
  { return x*x + y*y  + z*z + x*y*y*z; }
  constexpr double quadv3(double x, double y, double z)
  { return x*x + y*y  + z*z + x*y*z*z; }
}


class ASMCubeLag : public ASMCubeBase<ASMs3DLag>
{
public:
    static constexpr auto quad = std::array{
      quad1(0.0, 0.0, 0.0), quad1(0.5, 0.0, 0.0), quad1(1.0, 0.0, 0.0),
      quad1(0.0, 0.5, 0.0), quad1(0.5, 0.5, 0.0), quad1(1.0, 0.5, 0.0),
      quad1(0.0, 1.0, 0.0), quad1(0.5, 1.0, 0.0), quad1(1.0, 1.0, 0.0),

      quad1(0.0, 0.0, 0.5), quad1(0.5, 0.0, 0.5), quad1(1.0, 0.0, 0.5),
      quad1(0.0, 0.5, 0.5), quad1(0.5, 0.5, 0.5), quad1(1.0, 0.5, 0.5),
      quad1(0.0, 1.0, 0.5), quad1(0.5, 1.0, 0.5), quad1(1.0, 1.0, 0.5),

      quad1(0.0, 0.0, 1.0), quad1(0.5, 0.0, 1.0), quad1(1.0, 0.0, 1.0),
      quad1(0.0, 0.5, 1.0), quad1(0.5, 0.5, 1.0), quad1(1.0, 0.5, 1.0),
      quad1(0.0, 1.0, 1.0), quad1(0.5, 1.0, 1.0), quad1(1.0, 1.0, 1.0),
  };

  static constexpr auto quadv = std::array{
    quadv1(0.0, 0.0, 0.0), quadv2(0.0, 0.0, 0.0), quadv3(0.0, 0.0, 0.0),
    quadv1(0.5, 0.0, 0.0), quadv2(0.5, 0.0, 0.0), quadv3(0.5, 0.0, 0.0),
    quadv1(1.0, 0.0, 0.0), quadv2(1.0, 0.0, 0.0), quadv3(1.0, 0.0, 0.0),
    quadv1(0.0, 0.5, 0.0), quadv2(0.0, 0.5, 0.0), quadv3(0.0, 0.5, 0.0),
    quadv1(0.5, 0.5, 0.0), quadv2(0.5, 0.5, 0.0), quadv3(0.5, 0.5, 0.0),
    quadv1(1.0, 0.5, 0.0), quadv2(1.0, 0.5, 0.0), quadv3(1.0, 0.5, 0.0),
    quadv1(0.0, 1.0, 0.0), quadv2(0.0, 1.0, 0.0), quadv3(0.0, 1.0, 0.0),
    quadv1(0.5, 1.0, 0.0), quadv2(0.5, 1.0, 0.0), quadv3(0.5, 1.0, 0.0),
    quadv1(1.0, 1.0, 0.0), quadv2(1.0, 1.0, 0.0), quadv3(1.0, 1.0, 0.0),

    quadv1(0.0, 0.0, 0.5), quadv2(0.0, 0.0, 0.5), quadv3(0.0, 0.0, 0.5),
    quadv1(0.5, 0.0, 0.5), quadv2(0.5, 0.0, 0.5), quadv3(0.5, 0.0, 0.5),
    quadv1(1.0, 0.0, 0.5), quadv2(1.0, 0.0, 0.5), quadv3(1.0, 0.0, 0.5),
    quadv1(0.0, 0.5, 0.5), quadv2(0.0, 0.5, 0.5), quadv3(0.0, 0.5, 0.5),
    quadv1(0.5, 0.5, 0.5), quadv2(0.5, 0.5, 0.5), quadv3(0.5, 0.5, 0.5),
    quadv1(1.0, 0.5, 0.5), quadv2(1.0, 0.5, 0.5), quadv3(1.0, 0.5, 0.5),
    quadv1(0.0, 1.0, 0.5), quadv2(0.0, 1.0, 0.5), quadv3(0.0, 1.0, 0.5),
    quadv1(0.5, 1.0, 0.5), quadv2(0.5, 1.0, 0.5), quadv3(0.5, 1.0, 0.5),
    quadv1(1.0, 1.0, 0.5), quadv2(1.0, 1.0, 0.5), quadv3(1.0, 1.0, 0.5),

    quadv1(0.0, 0.0, 1.0), quadv2(0.0, 0.0, 1.0), quadv3(0.0, 0.0, 1.0),
    quadv1(0.5, 0.0, 1.0), quadv2(0.5, 0.0, 1.0), quadv3(0.5, 0.0, 1.0),
    quadv1(1.0, 0.0, 1.0), quadv2(1.0, 0.0, 1.0), quadv3(1.0, 0.0, 1.0),
    quadv1(0.0, 0.5, 1.0), quadv2(0.0, 0.5, 1.0), quadv3(0.0, 0.5, 1.0),
    quadv1(0.5, 0.5, 1.0), quadv2(0.5, 0.5, 1.0), quadv3(0.5, 0.5, 1.0),
    quadv1(1.0, 0.5, 1.0), quadv2(1.0, 0.5, 1.0), quadv3(1.0, 0.5, 1.0),
    quadv1(0.0, 1.0, 1.0), quadv2(0.0, 1.0, 1.0), quadv3(0.0, 1.0, 1.0),
    quadv1(0.5, 1.0, 1.0), quadv2(0.5, 1.0, 1.0), quadv3(0.5, 1.0, 1.0),
    quadv1(1.0, 1.0, 1.0), quadv2(1.0, 1.0, 1.0), quadv3(1.0, 1.0, 1.0),
  };

  explicit ASMCubeLag(unsigned char n_f = 3) : ASMCubeBase<ASMs3DLag>(n_f)
  {
  }
};

#endif
