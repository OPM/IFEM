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
  firstIp = 0;
  nGauss = 0;
  geo = 0;
}


ASMstruct::ASMstruct (const ASMstruct& patch, unsigned char n_f)
  : ASMbase(patch,n_f)
{
  firstIp = patch.nGauss;
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


void ASMstruct::getNoBouPoints (size_t& nPt, char ldim, char lindx)
{
  size_t nGp = 1;
  for (char d = 0; d < ldim; d++)
    nGp *= nGauss;

  firstBp[lindx] = nPt;

  nPt += this->getNoBoundaryElms(lindx,ldim)*nGp; // Includes 0-span elements
}
