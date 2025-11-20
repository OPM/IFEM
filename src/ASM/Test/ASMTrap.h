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
#include <tuple>
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

  std::tuple<double, std::vector<Vec3>, RealArray>
  getElementCorners(int iel) const
  {
     RealArray prm;
     std::vector<Vec3> XC;

     int p1, p2, n1, n2, dummy;
     this->getOrder(p1,p2,dummy);
     this->getSize(n1,n2,dummy);
     const int nel1 = n1 - p1 + 1;

     --iel;
     const int i1 = p1 + iel % nel1;
     const int i2 = p2 + iel / nel1;

     return {this->Base::getElementCorners(i1-1, i2-1, XC, &prm), XC, prm};
  }
};

#endif
