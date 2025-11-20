//==============================================================================
//!
//! \file ASMuTrap3.h
//!
//! \date Nov 11 2025
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief 3D trapezoid patch for LR unit-testing.
//!
//==============================================================================

#ifndef ASM_UTRAP3_H
#define ASM_UTRAP3_H

#include "Test/ASMTrap3.h"

#include <tuple>

template<class Base>
class ASMuTrap3 : public ASMTrap3<Base>
{
public:
  std::tuple<double, std::vector<Vec3>, RealArray>
  getElementCorners(int iel) const
  {
     RealArray prm;
     std::vector<Vec3> XC;
     return {this->Base::getElementCorners(iel, XC, &prm), XC, prm};
  }
};

#endif
