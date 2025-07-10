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
#include <sstream>

class ASMCube : public ASMs3D
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

  explicit ASMCube(unsigned char n_f = 3) : ASMs3D(n_f)
  {
    std::stringstream geo(cube);
    this->read(geo);
  }

  void shiftElemNumbers(int shift)
  {
    for (int& e : myMLGE)
      e += (e == -1 ? 0 : shift);
  }

  virtual ~ASMCube() {}
};


class ASMmxCube : public ASMs3Dmx
{
public:
  explicit ASMmxCube(const CharVec& n_f) : ASMs3Dmx(n_f)
  {
    std::stringstream geo(ASMCube::cube);
    this->read(geo);
  }
  virtual ~ASMmxCube() {}
};

#endif
