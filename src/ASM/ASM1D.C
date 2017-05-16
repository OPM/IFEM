// $Id$
//==============================================================================
//!
//! \file ASM1D.C
//!
//! \date May 15 2014
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Abstract interface for 1D patches.
//!
//==============================================================================

#include "ASM1D.h"
#include "ASMs1DC1.h"
#include "ASMs1DLag.h"
#include "ASMs1DSpec.h"
#include "Vec3Oper.h"


ASMbase* ASM1D::create (ASM::Discretization discretization, unsigned char nf)
{
  return create(discretization,1,nf);
}


ASMbase* ASM1D::create (ASM::Discretization discretization,
                        unsigned char nd, unsigned char nf)
{
  switch (discretization) {
  case ASM::SplineC1:
    return new ASMs1DC1(nd,nf);

  case ASM::Lagrange:
    return new ASMs1DLag(nd,nf);

  case ASM::Spectral:
    return new ASMs1DSpec(nd,nf);

  default:
    return new ASMs1D(nd,nf);
  }
}


#define TRY_CLONE1(classType,n) {				\
    const classType* p = dynamic_cast<const classType*>(this);	\
    if (p) return n ? new classType(*p,*n) : new classType(*p);	\
  }

ASMbase* ASM1D::clone (unsigned char* nf) const
{
  TRY_CLONE1(ASMs1DSpec,nf)
  TRY_CLONE1(ASMs1DLag,nf)
  TRY_CLONE1(ASMs1DC1,nf)
  TRY_CLONE1(ASMs1D,nf)

  std::cerr <<" *** ASM1D::clone: Failure, probably not a 1D patch"<< std::endl;
  return 0;
}

#undef TRY_CLONE1


double ASM1D::getElementSize (const std::vector<Vec3>& XC)
{
  // Find the element length
  return (XC.back() - XC.front()).length() / ASMbase::modelSize;
}
