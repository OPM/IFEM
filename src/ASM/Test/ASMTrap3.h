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
#include <tuple>
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

  std::tuple<double, std::vector<Vec3>, RealArray>
  getElementCorners(int iel) const
  {
     RealArray prm;
     std::vector<Vec3> XC;

     int n1, n2, n3;
     int p1,p2,p3;
     this->getOrder(p1,p2,p3);
     this->getSize(n1,n2,n3);
     const int nel1 = n1 - p1 + 1;
     const int nel2 = n2 - p2 + 1;

     --iel;
     int i1 = p1 + iel % nel1;
     int i2 = p2 + (iel / nel1) % nel2;
     int i3 = p3 + iel / (nel1*nel2);
     return {this->Base::getElementCorners(i1-1, i2-1, i3-1, XC, &prm), XC, prm};
  }
};

#endif
