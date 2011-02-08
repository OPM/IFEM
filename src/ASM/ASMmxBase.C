// $Id: ASMmxBase.C,v 1.1 2010-12-29 18:41:38 kmo Exp $
//==============================================================================
//!
//! \file ASMmxBase.C
//!
//! \date Dec 28 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for spline-based mixed finite element assembly drivers.
//!
//==============================================================================

#include "ASMmxBase.h"


ASMmxBase::ASMmxBase (unsigned char n_f1, unsigned char n_f2, bool geo1)
{
  nf1 = n_f1;
  nf2 = n_f2;

  geoUsesBasis1 = geo1;
}


void ASMmxBase::init (const std::vector<int>& MLGN, const int* sysMadof)
{
  MADOF.resize(MLGN.size());
  for (size_t i = 0; i < MADOF.size(); i++)
    MADOF[i] = sysMadof[MLGN[i]-1]-1;
}


void ASMmxBase::extrNodeVec (const Vector& globRes, Vector& nodeVec) const
{
  nodeVec.resize(nf1*nb1 + nf2*nb2);

  size_t i, j;
  int idof, ldof = 0;
  for (i = 0; i < nb1; i++)
  {
    idof = MADOF[i];
    for (j = 0; j < nf1; j++, ldof++)
      nodeVec[ldof] = globRes[idof++];
  }

  for (i = nb1; i < nb1+nb2; i++)
  {
    idof = MADOF[i];
    for (j = 0; j < nf2; j++, ldof++)
      nodeVec[ldof] = globRes[idof++];
  }
}
