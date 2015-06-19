// $Id$
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


char ASMmxBase::geoBasis            = 2;
ASMmxBase::MixedType ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS1;


ASMmxBase::ASMmxBase (unsigned char n_f1, unsigned char n_f2)
{
  nf1 = n_f1;
  nf2 = n_f2;
}


void ASMmxBase::initMx (const std::vector<int>& MLGN, const int* sysMadof)
{
  MADOF.resize(MLGN.size());
  for (size_t i = 0; i < MADOF.size(); i++)
    MADOF[i] = sysMadof[MLGN[i]-1]-1;
}


void ASMmxBase::extractNodeVecMx (const Vector& globRes, Vector& nodeVec,
				  int basis) const
{
  switch (basis) {
  case  1: nodeVec.resize(nf1*nb1); break;
  case  2: nodeVec.resize(nf2*nb2); break;
  default: nodeVec.resize(nf1*nb1 + nf2*nb2);
  }

  size_t i, j;
  int idof, ldof = 0;
  if (basis != 2)
    for (i = 0; i < nb1; i++)
    {
      idof = MADOF[i];
      for (j = 0; j < nf1; j++, ldof++)
	nodeVec[ldof] = globRes[idof++];
    }

  if (basis != 1)
    for (i = nb1; i < nb1+nb2; i++)
    {
      idof = MADOF[i];
      for (j = 0; j < nf2; j++, ldof++)
	nodeVec[ldof] = globRes[idof++];
    }
}


void ASMmxBase::injectNodeVecMx (Vector& globRes, const Vector& nodeVec,
                                 int basis) const
{
  size_t i, j;
  int idof, ldof = 0;
  if (basis != 2)
    for (i = 0; i < nb1; i++)
    {
      idof = MADOF[i];
      for (j = 0; j < nf1; j++, ldof++)
	globRes[idof++] = nodeVec[ldof];
    }

  if (basis != 1)
    for (i = nb1; i < nb1+nb2; i++)
    {
      idof = MADOF[i];
      for (j = 0; j < nf2; j++, ldof++)
	globRes[idof++] = nodeVec[ldof];
    }
}


bool ASMmxBase::getSolutionMx (Matrix& sField, const Vector& locSol,
			       const std::vector<int>& nodes) const
{
  if (nodes.empty()) return true;

  int low, high, nvar;
  if ((size_t)nodes.front() <= nb1)
  {
    nvar = nf1;
    low  = 1;
    high = nb1;
  }
  else
  {
    nvar = nf2;
    low  = nb1+1;
    high = nb1+nb2;
  }

  sField.resize(nvar,nodes.size());
  for (size_t i = 0; i < nodes.size(); i++)
    if (nodes[i] < low || nodes[i] > high)
    {
      std::cerr <<" *** ASMmxBase::getSolutionMx: Node #"<< nodes[i]
		<<" is out of range ["<< low <<","<< high <<"]."<< std::endl;
      return false;
    }
    else
    {
      int idof = nvar*(nodes[i]-1);
      if (low > 1) idof += nf1*nb1;
      for (int j = 0; j < nvar; j++)
	sField(j+1,i+1) = locSol[idof++];
    }

  return true;
}
