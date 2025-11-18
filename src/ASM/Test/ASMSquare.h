//==============================================================================
//!
//! \file ASMSquare.h
//!
//! \date Oct 2 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Default bi-unit square patch for unit-testing.
//!
//==============================================================================

#ifndef ASM_SQUARE_H
#define ASM_SQUARE_H

#include "ASMs2Dmx.h"
#include "ASMs2DLag.h"

#include <array>
#include <sstream>

template<class Base>
class ASMSquareBase : public Base
{
public:
  // formatting matches write routine, do not change
  static constexpr const char* square = "200 1 0 0\n"
                                        "2 0\n"
                                        "2 2\n"
                                        "0 0 1 1\n"
                                        "2 2\n"
                                        "0 0 1 1\n"
                                        "0 0\n"
                                        "1 0\n"
                                        "0 1\n"
                                        "1 1\n\n";

  static constexpr auto quad = std::array{
    0.0, 0.0, 1.0,
    0.0, 0.0, 1.5,
    1.0, 1.0, 3.0
  }; // x^2 + y^2 + x^2*y

  static constexpr auto quadv = std::array{
    0.0, 0.0,
    0.0, 0.0,
    1.0, 1.0,
    0.0, 0.0,
    0.0, 0.0,
    1.5, 1.0,
    1.0, 1.0,
    1.0, 1.5,
    3.0, 3.0,
  }; // {x^2 + y^2 + x^2y, x^2 + y^2 + xy^2}

  explicit ASMSquareBase(unsigned char n_f = 2) : Base(2,n_f)
  {
    std::stringstream geo(square);
    this->read(geo);
  }

  void shiftElemNumbers(int shift)
  {
    for (int& e : this->myMLGE)
      e += (e == -1 ? 0 : shift);
  }
};

using ASMSquare = ASMSquareBase<ASMs2D>;


class ASMmxSquare : public ASMs2Dmx
{
public:
  static constexpr auto quadv = std::array{
    // x
    0.0, 0.0, 1.0 / 3.0, 1.0,
    0.0, 0.0, 0.5,       1.5,
    1.0, 1.0, 5.0 / 3.0, 3.0,
    // y
    0.0,       0.0, 1.0,
    0.0,       0.0, 1.0,
    1.0 / 3.0, 0.5, 5.0 / 3.0,
    1.0,       1.5, 3.0,
    // p
     0.0,  0.0, 1.0,
     0.0, 0.25, 1.5,
    -1.0, -0.5, 1.0,
  }; // {x^2 + y^2 + x^2y, x^2 + y^2 + xy^2, x^2 - y^2 + x*y}

  explicit ASMmxSquare(const CharVec& n_f) : ASMs2Dmx(2,n_f)
  {
    std::stringstream geo(ASMSquare::square);
    this->read(geo);
  }

  virtual ~ASMmxSquare() {}
};


namespace {
  constexpr double quad1(double x, double y)
  { return x*x + y*y + x*x*y; }
  constexpr double quadv1(double x, double y)
  { return x*x + y*y + x*x*y; }
  constexpr double quadv2(double x, double y)
  { return x*x + y*y + x*y*y; }
}

class ASMSquareLag : public ASMSquareBase<ASMs2DLag>
{
public:
  static constexpr auto quad = std::array{
    quad1(0.0, 0.0), quad1(0.5, 0.0), quad1(1.0, 0.0),
    quad1(0.0, 0.5), quad1(0.5, 0.5), quad1(1.0, 0.5),
    quad1(0.0, 1.0), quad1(0.5, 1.0), quad1(1.0, 1.0),
  };

  static constexpr auto quadv = std::array{
    quadv1(0.0, 0.0), quadv2(0.0, 0.0),
    quadv1(0.5, 0.0), quadv2(0.5, 0.0),
    quadv1(1.0, 0.0), quadv2(1.0, 0.0),

    quadv1(0.0, 0.5), quadv2(0.0, 0.5),
    quadv1(0.5, 0.5), quadv2(0.5, 0.5),
    quadv1(1.0, 0.5), quadv2(1.0, 0.5),

    quadv1(0.0, 1.0), quadv2(0.0, 1.0),
    quadv1(0.5, 1.0), quadv2(0.5, 1.0),
    quadv1(1.0, 1.0), quadv2(1.0, 1.0),
  };

  explicit ASMSquareLag(unsigned char n_f = 2)
    : ASMSquareBase<ASMs2DLag>(n_f)
  {
  }
};

#endif
