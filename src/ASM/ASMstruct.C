// $Id$
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


ASMstruct::ASMstruct (const ASMstruct& patch, unsigned char n_f)
  : ASMbase(patch,n_f)
{
  nGauss = patch.nGauss;
  geo = patch.geo;
}


ASMstruct::~ASMstruct ()
{
  if (geo && !shareFE) delete geo;
}


void ASMstruct::getNoIntPoints (size_t& nPt)
{
  size_t nGp = 1;
  for (unsigned char d = 0; d < ndim; d++)
    nGp *= nGauss;

  firstIp = nPt;

  nPt += this->getNoElms(true)*nGp; // Note: Includes also the 0-span elements
}


void ASMstruct::getNoBouPoints (size_t& nPt, int ldim, int lindx)
{
  firstBp = nPt;
  if (ldim == 0) nPt ++; // This is correct for 1D only, 2D and 3D must take
  // care of how many element do we have on the boundary defined by lindx
}
