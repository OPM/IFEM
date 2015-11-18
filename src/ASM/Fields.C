// $Id$
//==============================================================================
//!
//! \file Fields.C
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Base class for vector fields.
//!
//==============================================================================

#include "SplineFields2D.h"
#include "SplineFields3D.h"
#include "LagrangeFields2D.h"
#include "LagrangeFields3D.h"
#include "ASMs2DLag.h"
#include "ASMs3DLag.h"


Fields* Fields::create (const ASMbase* pch, const RealArray& v,
			char basis, const char* name)
{
  const ASMs2DLag* pl2 = dynamic_cast<const ASMs2DLag*>(pch);
  if (pl2) return new LagrangeFields2D(pl2,v,basis,name);

  const ASMs2D* ps2 = dynamic_cast<const ASMs2D*>(pch);
  if (ps2) return new SplineFields2D(ps2,v,basis,name);

  const ASMs3DLag* pl3 = dynamic_cast<const ASMs3DLag*>(pch);
  if (pl3) return new LagrangeFields3D(pl3,v,basis,name);

  const ASMs3D* ps3 = dynamic_cast<const ASMs3D*>(pch);
  if (ps3) return new SplineFields3D(ps3,v,basis,name);

  return nullptr;
}
