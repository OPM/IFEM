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


ASMmxBase::ASMmxBase (const std::vector<unsigned char>& n_f)
{
  nfx = n_f;
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
  if (basis > (int)nfx.size())
    basis = 0;

  size_t len=0;
  if (basis == 0)
    for (size_t i=0;i<nfx.size();++i)
      len += nfx[i]*nb[i];
  else
    len = nfx[basis-1]*nb[basis-1];

  nodeVec.resize(len);

  size_t i, j;
  int idof, ldof = 0;
  size_t k=basis==0?1:basis;
  size_t ofs = std::accumulate(nb.begin(), nb.begin()+k-1, 0);
  for (; k < (basis==0?nfx.size()+1:(size_t)basis+1); ++k) {
    for (i = ofs; i < nb[k-1]+ofs; i++)
    {
      idof = MADOF[i];
      for (j = 0; j < nfx[k-1]; j++, ldof++)
	nodeVec[ldof] = globRes[idof++];
    }
    ofs += nb[k-1];
  }
}


void ASMmxBase::injectNodeVecMx (Vector& globRes, const Vector& nodeVec,
                                 int basis) const
{
  if (basis > (int)nfx.size())
    basis = 0;

  size_t i, j;
  int idof, ldof = 0;
  size_t k=basis==0?1:basis;
  size_t ofs = std::accumulate(nb.begin(), nb.begin()+k-1, 0);
  for (; k < (basis==0?nfx.size()+1:(size_t)basis+1); ++k) {
    for (i = ofs; i < nb[k-1]+ofs; i++)
    {
      idof = MADOF[i];
      for (j = 0; j < nfx[k-1]; j++, ldof++)
	globRes[idof++] = nodeVec[ldof];
    }
  }
}


bool ASMmxBase::getSolutionMx (Matrix& sField, const Vector& locSol,
			       const std::vector<int>& nodes) const
{
  if (nodes.empty()) return true;

  int low, high, nvar;
  if ((size_t)nodes.front() <= nb[0])
  {
    nvar = nfx[0];
    low  = 1;
    high = nb[0];
  }
  else
  {
    nvar = nfx[1];
    low  = nb[0]+1;
    high = nb[0]+nb[1];
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
      if (low > 1) idof += nfx[0]*nb[0];
      for (int j = 0; j < nvar; j++)
	sField(j+1,i+1) = locSol[idof++];
    }

  return true;
}
