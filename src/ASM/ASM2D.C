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
#include "ASMs2DTri.h"
#include "ASMu2DLag.h"
#include "ASMs2Dmx.h"
#include "ASMs2DmxLag.h"
#include "ASMs2DSpec.h"
#ifdef HAS_LRSPLINE
#include "LR/ASMu2DIB.h"
#include "LR/ASMu2Dmx.h"
#include "ASMu2Dnurbs.h"
#endif
#include "Vec3Oper.h"


ASMbase* ASM2D::create (ASM::Discretization discretization, unsigned char nf)
{
  return ASM2D::create(discretization,2,CharVec(1,nf),false);
}


ASMbase* ASM2D::create (ASM::Discretization discretization,
                        unsigned char nd, const CharVec& nf, bool mixedFEM)
{
  switch (discretization) {
  case ASM::SplineC1:
    return new ASMs2DC1(nd,nf.front());

  case ASM::Lagrange:
    if (nf.size() > 1 && nf[1] > 64) // hack for mesh input from file
      return new ASMu2DLag(nd,nf[0],nf[1]);
    else if (nf.size() > 1 || mixedFEM)
      return new ASMs2DmxLag(nd,nf);
    else
      return new ASMs2DLag(nd,nf.front());

  case ASM::Triangle:
    return new ASMs2DTri(nd,nf.front());

  case ASM::Spectral:
    return new ASMs2DSpec(nd,nf.front());

#ifdef HAS_LRSPLINE
  case ASM::LRSpline:
    if (nf.size() > 1 && nf[1] == 'I') // hack for immersed boundary approach
      return new ASMu2DIB(nd,nf[0],nf[2]);
    else if (nf.size() > 1 || mixedFEM)
      return new ASMu2Dmx(nd,nf);
    else
      return new ASMu2D(nd,nf.front());
  case ASM::LRNurbs:
    return new ASMu2Dnurbs(nd,nf.front());
#endif

  default:
    if (nf.size() > 1 && nf[1] == 'I') // hack for immersed boundary approach
      return new ASMs2DIB(nd,nf[0],nf[2]);
    else if (nf.size() > 1 || mixedFEM)
      return new ASMs2Dmx(nd,nf);
    else
      return new ASMs2D(nd,nf.front());
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

ASMbase* ASM2D::clone (const CharVec& nf) const
{
  TRY_CLONE2(ASMs2DmxLag,nf)
  TRY_CLONE2(ASMs2Dmx,nf)
  TRY_CLONE1(ASMs2DSpec,nf)
  TRY_CLONE1(ASMs2DTri,nf)
  TRY_CLONE1(ASMu2DLag,nf)
  TRY_CLONE1(ASMs2DLag,nf)
  TRY_CLONE1(ASMs2DC1,nf)
  TRY_CLONE1(ASMs2DIB,nf)
  TRY_CLONE1(ASMs2D,nf)
#ifdef HAS_LRSPLINE
  TRY_CLONE2(ASMu2Dmx,nf)
  TRY_CLONE1(ASMu2DIB,nf)
  TRY_CLONE1(ASMu2D,nf)
#endif

  std::cerr <<" *** ASM2D::clone: Failure, probably not a 2D patch"<< std::endl;
  return 0;
}

#undef TRY_CLONE1
#undef TRY_CLONE2


double ASM2D::getElementSize (const std::vector<Vec3>& XC)
{
  // Find longest diagonal (diameter av minste omskrivende sirkel)
  double siz = std::max((XC[3]-XC[0]).length(),(XC[2]-XC[1]).length());
  return siz / ASMbase::modelSize; // Scale down by model dimension
}
