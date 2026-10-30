//==============================================================================
//!
//! \file ASMTrap3.h
//!
//! \date Nov 11 2025
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief 3D trapezoidal patch for unit-testing.
//!
//==============================================================================

#ifndef ASM_TRAP3_H
#define ASM_TRAP3_H

#include "ASMbase.h"
#include "Vec3.h"

#include <sstream>
#include <vector>


template<class Base>
class ASMTrap3 : public Base
{
public:
  ASMTrap3()
  {
    constexpr const char* trap =
      "700 1 0 0\n"
      "3 0\n"
      "2 2\n"
      "0 0 1 1\n"
      "2 2\n"
      "0 0 1 1\n"
      "2 2\n"
      "0 0 1 1\n"
      "0 0 0\n"
      "6 0 0\n"
      "1 3 0\n"
      "4 3 0\n"
      "1 0 1\n"
      "5 0 1\n"
      "2 3 1\n"
      "4 3 1\n\n";

    ASMbase::resetNumbering();
    std::stringstream str;
    str << trap;
    this->read(str);
  }

};

#endif
