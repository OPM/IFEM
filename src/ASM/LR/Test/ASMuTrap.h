//==============================================================================
//!
//! \file ASMuTrap.h
//!
//! \date Nov 11 2025
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief 2D trapezoid patch for LR unit-testing.
//!
//==============================================================================

#ifndef ASM_UTRAP_H
#define ASM_UTRAP_H

#include "Test/ASMTrap.h"

#include <tuple>

template<class Base>
class ASMuTrap : public ASMTrap<Base>
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
