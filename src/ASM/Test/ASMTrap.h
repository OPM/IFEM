//==============================================================================
//!
//! \file ASMTrap.h
//!
//! \date Nov 11 2025
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief 2D trapezoid patch for unit-testing.
//!
//==============================================================================

#ifndef ASM_TRAP_H
#define ASM_TRAP_H

#include "ASMbase.h"
#include "Vec3.h"

#include <sstream>
#include <vector>


template<class Base>
class ASMTrap : public Base
{
public:
  ASMTrap()
  {
    constexpr const char* trap =
      "200 1 0 0\n"
      "2 0\n"
      "2 2\n"
      "0 0 1 1\n"
      "2 2\n"
      "0 0 1 1\n"
      "0.0 0.0\n"
      "6.0 0.0\n"
      "1.0 3.0\n"
      "4.0 3.0\n";

    ASMbase::resetNumbering();
    std::stringstream str;
    str << trap;
    this->read(str);
  }

};

#endif
