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
#include <sstream>

static const char* square = "200 1 0 0 2 0\n"
                            "2 2 0 0 1 1\n"
                            "2 2 0 0 1 1\n"
                            "0 0 1 0 0 1 1 1\n";

class ASMSquare : public ASMs2D
{
public:
  explicit ASMSquare(unsigned char n_f = 2) : ASMs2D(2,n_f)
  {
    std::stringstream geo(square);
    this->read(geo);
  }
  virtual ~ASMSquare() {}
};


class ASMmxSquare : public ASMs2Dmx
{
public:
  explicit ASMmxSquare(const CharVec& n_f) : ASMs2Dmx(2,n_f)
  {
    std::stringstream geo(square);
    this->read(geo);
    this->generateFEMTopology();
  }
  virtual ~ASMmxSquare() {}
};

#endif
