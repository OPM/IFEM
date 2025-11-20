//==============================================================================
//!
//! \file ASMuCube.h
//!
//! \date Oct 2 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Default tri-unit cube LR patch for unit-testing.
//!
//==============================================================================

#ifndef _ASM_UCUBE_H
#define _ASM_UCUBE_H

#include "LR/ASMu3Dmx.h"
#include "Test/ASMCube.h"
#include <sstream>

class ASMuCube : public ASMu3D
{
public:
  // formatting matches write routine, do not change
  static constexpr const char* cube =
    "# LRSPLINE VOLUME\n"
    "#\tp1\tp2\tp3\tNbasis\tNline\tNel\tdim\trat\n"
    "\t2\t2\t2\t8\t6\t1\t3\t0\n"
    "# Basis functions:\n"
    "0: [0 0 1 ] x [0 0 1 ] x [0 0 1 ] 0 0 0 (1)\n"
    "1: [0 1 1 ] x [0 0 1 ] x [0 0 1 ] 1 0 0 (1)\n"
    "2: [0 0 1 ] x [0 1 1 ] x [0 0 1 ] 0 1 0 (1)\n"
    "3: [0 1 1 ] x [0 1 1 ] x [0 0 1 ] 1 1 0 (1)\n"
    "4: [0 0 1 ] x [0 0 1 ] x [0 1 1 ] 0 0 1 (1)\n"
    "5: [0 1 1 ] x [0 0 1 ] x [0 1 1 ] 1 0 1 (1)\n"
    "6: [0 0 1 ] x [0 1 1 ] x [0 1 1 ] 0 1 1 (1)\n"
    "7: [0 1 1 ] x [0 1 1 ] x [0 1 1 ] 1 1 1 (1)\n"
    "# Mesh rectangles:\n"
    "[0, 0] x [0, 1] x [0, 1] (2)\n"
    "[1, 1] x [0, 1] x [0, 1] (2)\n"
    "[0, 1] x [0, 0] x [0, 1] (2)\n"
    "[0, 1] x [1, 1] x [0, 1] (2)\n"
    "[0, 1] x [0, 1] x [0, 0] (2)\n"
    "[0, 1] x [0, 1] x [1, 1] (2)\n"
    "# Elements:\n"
    "0 [3] : (0, 0, 0) x (1, 1, 1)    {0, 1, 2, 3, 4, 5, 6, 7}\n";

  static constexpr auto quad = ASMCube::quad;
  static constexpr auto quadv = ASMCube::quadv;

  explicit ASMuCube(unsigned char n_f = 3,
                    double xshift = 0.0,
                    double yshift = 0.0,
                    double zshift = 0.0)
    : ASMu3D(n_f)
  {
    auto&& addLine = [](double x, double y, double z)
    {
      return std::to_string(x) + " " +
             std::to_string(y) + " " +
             std::to_string(z) + '\n';
    };
    std::stringstream geo("700 1 0 0\n3 0\n"
                          "2 2\n0 0 1 1\n"
                          "2 2\n0 0 1 1\n"
                          "2 2\n0 0 1 1\n" +
                          addLine(xshift,       yshift,       zshift) +
                          addLine(xshift + 1.0, yshift,       zshift) +
                          addLine(xshift,       yshift + 1.0, zshift) +
                          addLine(xshift + 1.0, yshift + 1.0, zshift) +
                          addLine(xshift,       yshift,       zshift + 1.0) +
                          addLine(xshift + 1.0, yshift,       zshift + 1.0) +
                          addLine(xshift,       yshift + 1.0, zshift + 1.0) +
                          addLine(xshift + 1.0, yshift + 1.0, zshift + 1.0));
    this->read(geo);
  }

  void shiftElemNumbers(int shift)
  {
    for (int& e : myMLGE)
      e += (e == -1 ? 0 : shift);
  }

  std::array<int,4> getFaceCorners(int dir)
  {
    DirichletFace df(this->getBasis(),dir);
    std::array<int,4> corners;
    std::copy(df.corners, df.corners + 4, corners.begin());
    return corners;
  }
};


class ASMmxuCube : public ASMu3Dmx
{
public:
  static constexpr auto quadv = ASMmxCube::quadv;
  explicit ASMmxuCube(const CharVec& n_f) : ASMu3Dmx(n_f)
  {
    std::stringstream geo(ASMCube::cube);
    this->read(geo);
  }
};

#endif
