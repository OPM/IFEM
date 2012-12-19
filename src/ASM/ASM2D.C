// $Id$
//==============================================================================
//!
//! \file ASM2D.C
//!
//! \date Oct 23 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Abstract interface for 2D patches.
//!
//==============================================================================

#include "ASM2D.h"
#include "ASMs2DIB.h"
#include "ASMs2DC1.h"
#include "ASMs2Dmx.h"
#include "ASMs2DmxLag.h"
#include "ASMs2DSpec.h"
#ifdef HAS_LRSPLINE
#include "LR/ASMu2D.h"
#endif


ASMbase* ASM2D::create (ASM::Discretization discretization,
			const unsigned char* nf, bool mixedFEM)
{
  switch (discretization) {
  case ASM::SplineC1:
    return new ASMs2DC1(2,nf[0]);

  case ASM::Lagrange:
    if (nf[1] > 0 || mixedFEM)
      return new ASMs2DmxLag(2,nf[0],nf[1]);
    else
      return new ASMs2DLag(2,nf[0]);

  case ASM::Spectral:
    return new ASMs2DSpec(2,nf[0]);

#ifdef HAS_LRSPLINE
  case ASM::LRSpline:
    return new ASMu2D(2,nf[0]);
#endif

  default:
    if (nf[1] == 'I') // hack for immersed boundary approach
      return new ASMs2DIB(2,nf[0],nf[2]);
    else if (nf[1] > 0 || mixedFEM)
      return new ASMs2Dmx(2,nf[0],nf[1]);
    else
      return new ASMs2D(2,nf[0]);
  }
}


#define TRY_CLONE1(classType,n) {					\
    const classType* p = dynamic_cast<const classType*>(this);		\
    if (p) return n ? new classType(*p,n[0]) : new classType(*p);	\
  }
#define TRY_CLONE2(classType,n) {					\
    const classType* p = dynamic_cast<const classType*>(this);		\
    if (p) return n ? new classType(*p,n[0],n[1]) : new classType(*p);	\
  }

ASMbase* ASM2D::clone (unsigned char* nf) const
{
  TRY_CLONE2(ASMs2DmxLag,nf)
  TRY_CLONE2(ASMs2Dmx,nf)
  TRY_CLONE1(ASMs2DSpec,nf)
  TRY_CLONE1(ASMs2DLag,nf)
  TRY_CLONE1(ASMs2DC1,nf)
  TRY_CLONE1(ASMs2DIB,nf)
  TRY_CLONE1(ASMs2D,nf)
#ifdef HAS_LRSPLINE
  TRY_CLONE1(ASMu2D,nf)
#endif

  std::cerr <<" *** ASM2D::clone: Failure, probably not a 2D patch"<< std::endl;
  return 0;
}

#undef TRY_CLONE1
#undef TRY_CLONE2
