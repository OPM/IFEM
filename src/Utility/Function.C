// $Id$
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
#include "TensorFunction.h"
#include "Vec3Oper.h"


Vec3 PressureField::evaluate (const Vec3& x, const Vec3& n) const
{
  if (!pressure) // zero pressure field
    return Vec3();

  const RealFunc& p = *pressure;

  if (pdir < 1) // normal pressure
    return p(x) * n;

  Vec3 t; // pressure acting in global pdir direction
  t[(pdir-1)%3] = p(x);

  if (pdir > 3) // normal pressure in global pdir direction
    t = (t*n) * n;

  return t;
}


TractionField::TractionField (const STensorFunc& field)
{
  sigma = &field;
  sigmaN = nullptr;
}


TractionField::TractionField (const TensorFunc& field)
{
  sigma = nullptr;
  sigmaN = &field;
}


Vec3 TractionField::evaluate (const Vec3& x, const Vec3& n) const
{
  if (sigma) // symmetric tensor field
    return (*sigma)(x) * n;
  else if (sigmaN) // non-symmetric tensor field
    return (*sigmaN)(x) * n;
  else // zero tensor field
    return Vec3();
}


bool TractionField::isZero () const
{
  if (sigma)
    return sigma->isZero();
  else if (sigmaN)
    return sigmaN->isZero();
  else
    return true;
}
