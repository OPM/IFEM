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
#include "SplineFields2Dmx.h"
#include "SplineFields3D.h"
#include "SplineFields3Dmx.h"
#include "LagrangeFields2D.h"
#include "LagrangeFields3D.h"
#include "ASMs2DLag.h"
#include "ASMs3DLag.h"
#include "ASMs2Dmx.h"
#include "ASMs3Dmx.h"
#ifdef HAS_LRSPLINE
#include "LR/LRSplineFields2D.h"
#include "LR/LRSplineFields2Dmx.h"
#include "LR/LRSplineFields3D.h"
#include "LR/LRSplineFields3Dmx.h"
#include "LR/ASMu2D.h"
#include "LR/ASMu2Dmx.h"
#include "LR/ASMu3D.h"
#include "LR/ASMu3Dmx.h"
#endif


Fields* Fields::create (const ASMbase* pch, const RealArray& v,
			char basis, int nf, const char* name)
{
#ifdef HAS_LRSPLINE
  const ASMu2Dmx* pu2mx = dynamic_cast<const ASMu2Dmx*>(pch);
  if (basis > 10 && pu2mx) return new LRSplineFields2Dmx(pu2mx,v,basis,name);

  const ASMu2D* pu2 = dynamic_cast<const ASMu2D*>(pch);
  if (pu2) return new LRSplineFields2D(pu2,v,basis,nf,name);

  const ASMu3Dmx* pu3mx = dynamic_cast<const ASMu3Dmx*>(pch);
  if (basis > 10 && pu3mx) return new LRSplineFields3Dmx(pu3mx,v,basis,name);

  const ASMu3D* pu3 = dynamic_cast<const ASMu3D*>(pch);
  if (pu3) return new LRSplineFields3D(pu3,v,basis,nf,name);
#endif
  const ASMs2DLag* pl2 = dynamic_cast<const ASMs2DLag*>(pch);
  if (pl2) return new LagrangeFields2D(pl2,v,basis,name);

  const ASMs2Dmx* ps2mx = dynamic_cast<const ASMs2Dmx*>(pch);
  if (basis > 10 && ps2mx) return new SplineFields2Dmx(ps2mx,v,basis,name);

  const ASMs2D* ps2 = dynamic_cast<const ASMs2D*>(pch);
  if (ps2) return new SplineFields2D(ps2,v,basis,nf,name);

  const ASMs3DLag* pl3 = dynamic_cast<const ASMs3DLag*>(pch);
  if (pl3) return new LagrangeFields3D(pl3,v,basis,name);

  const ASMs3Dmx* ps3mx = dynamic_cast<const ASMs3Dmx*>(pch);
  if (basis > 10 && ps3mx) return new SplineFields3Dmx(ps3mx,v,basis,name);

  const ASMs3D* ps3 = dynamic_cast<const ASMs3D*>(pch);
  if (ps3) return new SplineFields3D(ps3,v,basis,nf,name);

  return nullptr;
}
