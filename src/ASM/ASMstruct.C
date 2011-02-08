// $Id: ASMstruct.C,v 1.4 2010-03-24 12:27:07 rho Exp $
//==============================================================================
//!
//! \file ASMstruct.C
//!
//! \date Sep 20 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for structured spline-based FE assembly drivers.
//!
//==============================================================================

#include "ASMstruct.h"
#include "GoTools/geometry/GeomObject.h"


int ASMstruct::gEl = 0;
int ASMstruct::gNod = 0;


ASMstruct::ASMstruct (unsigned char n_p, unsigned char n_s, unsigned char n_f)
  : ASMbase(n_p,n_s,n_f)
{
  geo = 0;
}


ASMstruct::~ASMstruct ()
{
  if (geo) delete geo;
}
