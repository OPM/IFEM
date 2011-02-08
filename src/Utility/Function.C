// $Id: Function.C,v 1.6 2010-10-14 19:10:55 kmo Exp $
//==============================================================================
//!
//! \file Function.C
//!
//! \date Jul 1 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief General functions with arbitrary argument and value type.
//!
//==============================================================================

#include "Function.h"
#include "Tensor.h"
#include "Vec3.h"
#include "Vec3Oper.h"


bool PressureField::isNormalPressure () const
{
  return pdir < 1 || pdir > 3;
}


Vec3 PressureField::evaluate (const Vec3& x, const Vec3& n) const
{
  const RealFunc& p = *pressure;

  if (pdir < 1) // normal pressure
    return p(x) * n;

  Vec3 t; // pressure acting in global pdir direction
  t[(pdir-1)%3] = p(x);

  if (pdir > 3) // normal pressure in global pdir direction
    t = (t*n) * n;

  return t;
}


Vec3 TractionField::evaluate (const Vec3& x, const Vec3& n) const
{
  return stress(x) * n;
}
