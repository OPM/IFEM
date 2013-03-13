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
std::map<int,int> ASMstruct::xNode;


ASMstruct::ASMstruct (unsigned char n_p, unsigned char n_s, unsigned char n_f)
  : ASMbase(n_p,n_s,n_f)
{
  geo = NULL;
}


ASMstruct::ASMstruct (const ASMstruct& patch, unsigned char n_f)
  : ASMbase(patch,n_f)
{
  geo = patch.geo;
}


ASMstruct::~ASMstruct ()
{
  if (geo && !shareFE) delete geo;
}


void ASMstruct::resetNumbering (int nnod)
{
  gEl = 0;
  gNod = nnod;
  xNode.clear();
}
