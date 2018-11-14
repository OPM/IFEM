// $Id$
//==============================================================================
//!
//! \file ASMunstruct.C
//!
//! \date December 2010
//!
//! \author Kjetil Andre Johannessen / SINTEF
//!
//! \brief Base class for unstructured spline-based FE assembly drivers.
//!
//==============================================================================

#include "ASMunstruct.h"


int ASMunstruct::gEl = 0;
int ASMunstruct::gNod = 0;


ASMunstruct::ASMunstruct (unsigned char n_p, unsigned char n_s,
                          unsigned char n_f)
  : ASMbase(n_p,n_s,n_f)
{
}


ASMunstruct::ASMunstruct (const ASMunstruct& patch, unsigned char n_f)
  : ASMbase(patch,n_f)
{
  nGauss = patch.nGauss;
}
