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

static const char* cube = "700 1 0 0 3 0\n"
                          "2 2 0 0 1 1\n"
                          "2 2 0 0 1 1\n"
                          "2 2 0 0 1 1\n"
                          "0 0 0  1 0 0  0 1 0  1 1 0\n"
                          "0 0 1  1 0 1  0 1 1  1 1 1\n";


class ASMCube : public ASMs3D
{
public:
  explicit ASMCube(unsigned char n_f = 3) : ASMs3D(n_f)
  {
    std::stringstream geo(cube);
    this->read(geo);
  }
  virtual ~ASMCube() {}
};


class ASMmxCube : public ASMs3Dmx
{
public:
  explicit ASMmxCube(const CharVec& n_f) : ASMs3Dmx(n_f)
  {
    std::stringstream geo(cube);
    this->read(geo);
  }
  virtual ~ASMmxCube() {}
};

#endif
