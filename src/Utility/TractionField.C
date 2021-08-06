// $Id$
//==============================================================================
//!
//! \file TractionField.C
//!
//! \date Jul 1 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Function interfaces for representation of explicit traction fields.
//!
//==============================================================================

#include "TractionField.h"
#include "TensorFunction.h"
#include "Vec3Oper.h"


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


Vec3 PressureField::evaluate (const Vec3& x, const Vec3& n) const
{
  if (!pressure) // zero pressure field
    return Vec3();

  const RealFunc& p = *pressure;

  if (pdir < 1 && !pdfn) // normal pressure
    return p(x) * n;

  Vec3 t;
  if (pdfn) // pressure direction specified as a function
  {
    t = (*pdfn)(x);
    t.normalize();
    t *= p(x);
  }
  else if (pdir > 0) // pressure acting in global pdir direction
    t[(pdir-1)%3] = p(x);

  if (pdir > 3) // normal pressure in global pdir direction
    t = (t*n) * n;

  return t;
}


ForceDirField::~ForceDirField ()
{
  delete force;
  delete fdir;
  delete shape;
}


Vec3 ForceDirField::evaluate (const Vec3& x, const Vec3&) const
{
  if (!force || !shape) // zero force
    return Vec3();

  const ScalarFunc& F = *force;
  const RealFunc& Shp = *shape;

  const Vec4* Xt = dynamic_cast<const Vec4*>(&x);
  Real time = Xt ? Xt->t : Real(0);

  Vec3   Fdir;
  Tensor Trot(3,true);
  if (fdir)
  {
    Fdir = (*fdir)(time);
    if (dirVec) // fdir gives the force direction
      Trot = Tensor(Tlg[2],Fdir,true);
    else // fdir gives rotation angles w.r.t. the local axes
    {
      Trot = Tensor(Fdir.x,Fdir.y,Fdir.z);
      Fdir = Trot[1];
    }
  }

  // Calculate updated local coordinates of evaluation point
  Vec4 Xloc(Trot * ((x-X0) * Tlg), time);
  // Evaluate traction magnitude at current point
  Real Fa = F(time)*Shp(Xloc);
  // Evaluate traction vector in local axes
  Vec3 Fv = Fa*Fdir;
  // Transform to global coordinates
  return Tlg*Fv;
}
