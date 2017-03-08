// $Id$
//==============================================================================
//!
//! \file ASM3D.C
//!
//! \date January 2013
//!
//! \author Kjetil A. Johannessen / NTNU
//!
//! \brief Abstract interface for 3D patches.
//!
//==============================================================================

#include "ASM3D.h"
#include "ASMs3DSpec.h"
#include "ASMs3Dmx.h"
#include "ASMs3DmxLag.h"
#ifdef HAS_LRSPLINE
#include "LR/ASMu3D.h"
#include "LR/ASMu3Dmx.h"
#endif
#include "Vec3Oper.h"


ASMbase* ASM3D::create (ASM::Discretization discretization, unsigned char nf)
{
  return ASM3D::create(discretization,CharVec(1,nf));
}


ASMbase* ASM3D::create (ASM::Discretization discretization,
                        const CharVec& nf, bool mixedFEM)
{
  switch (discretization) {

  case ASM::Lagrange:
    if (nf.size() > 1 || mixedFEM)
      return new ASMs3DmxLag(nf);
    else
      return new ASMs3DLag(nf[0]);

  case ASM::Spectral:
    return new ASMs3DSpec(nf[0]);

#ifdef HAS_LRSPLINE
  case ASM::LRSpline:
    if (nf.size() > 1 || mixedFEM)
      return new ASMu3Dmx(nf);
    else
      return new ASMu3D(nf[0]);
#endif

  default:
    if (nf.size() > 1 || mixedFEM)
      return new ASMs3Dmx(nf);
    else
      return new ASMs3D(nf[0]);
  }
}


#define TRY_CLONE1(classType,n) {					\
    const classType* p = dynamic_cast<const classType*>(this);		\
    if (p) return n.empty() ? new classType(*p) : new classType(*p,n[0]);\
  }
#define TRY_CLONE2(classType,n) {					\
    const classType* p = dynamic_cast<const classType*>(this);		\
    if (p) return n.empty() ? new classType(*p) : new classType(*p,n);	\
  }

ASMbase* ASM3D::clone (const CharVec& nf) const
{
  TRY_CLONE2(ASMs3DmxLag,nf)
  TRY_CLONE2(ASMs3Dmx,nf)
  TRY_CLONE1(ASMs3DSpec,nf)
  TRY_CLONE1(ASMs3DLag,nf)
  TRY_CLONE1(ASMs3D,nf)
#ifdef HAS_LRSPLINE
  TRY_CLONE1(ASMu3D,nf)
#endif

  std::cerr <<" *** ASM3D::clone: Failure, probably not a 3D patch"<< std::endl;
  return 0;
}

#undef TRY_CLONE1
#undef TRY_CLONE2


double ASM3D::getElementSize (const std::vector<Vec3>& XC)
{
  // Find the longest diagonal (diameter av minste omskrivende kule)
  double siz = (XC[7] - XC[0]).length();
  for (int c = 1; c < 4; c++)
    siz = std::max(siz,(XC[7-c]-XC[c]).length());

  return siz / ASMbase::modelSize; // Scale down by model size
}
