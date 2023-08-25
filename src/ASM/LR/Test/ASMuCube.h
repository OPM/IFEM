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

#include "ASMu3Dmx.h"
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

  explicit ASMuCube(unsigned char n_f = 3) : ASMu3D(n_f)
  {
    std::stringstream geo("700 1 0 0\n3 0\n"
                          "2 2\n0 0 1 1\n"
                          "2 2\n0 0 1 1\n"
                          "2 2\n0 0 1 1\n"
                          "0 0 0\n1 0 0\n0 1 0\n1 1 0\n"
                          "0 0 1\n1 0 1\n0 1 1\n1 1 1\n");
    this->read(geo);
  }
  virtual ~ASMuCube() {}
};


class ASMmxuCube : public ASMu3Dmx
{
public:
  explicit ASMmxuCube(const CharVec& n_f) : ASMu3Dmx(n_f)
  {
    std::stringstream geo("700 1 0 0\n3 0\n"
                          "2 2\n0 0 1 1\n"
                          "2 2\n0 0 1 1\n"
                          "2 2\n0 0 1 1\n"
                          "0 0 0\n1 0 0\n0 1 0\n1 1 0\n"
                          "0 0 1\n1 0 1\n0 1 1\n1 1 1\n");
    this->read(geo);
  }
  virtual ~ASMmxuCube() {}
};

#endif
