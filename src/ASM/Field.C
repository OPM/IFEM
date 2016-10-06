// $Id$
//==============================================================================
//!
//! \file Field.C
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Base class for scalar fields.
//!
//==============================================================================

#include "SplineField2D.h"
#include "SplineField3D.h"
#include "LagrangeField2D.h"
#include "LagrangeField3D.h"
#include "ASMs2DLag.h"
#include "ASMs3DLag.h"


Field* Field::create (const ASMbase* pch, const RealArray& v,
                      char basis, char cmp, const char* name)
{
  const ASMs2DLag* pl2 = dynamic_cast<const ASMs2DLag*>(pch);
  if (pl2) return new LagrangeField2D(pl2,v,basis,cmp,name);

  const ASMs2D* ps2 = dynamic_cast<const ASMs2D*>(pch);
  if (ps2) return new SplineField2D(ps2,v,basis,cmp,name);

  const ASMs3DLag* pl3 = dynamic_cast<const ASMs3DLag*>(pch);
  if (pl3) return new LagrangeField3D(pl3,v,basis,cmp,name);

  const ASMs3D* ps3 = dynamic_cast<const ASMs3D*>(pch);
  if (ps3) return new SplineField3D(ps3,v,basis,cmp,name);

  return nullptr;
}
