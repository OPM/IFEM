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
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "Integrand.h"

#include "GoTools/geometry/GeomObject.h"


ASMstruct::ASMstruct (unsigned char n_p, unsigned char n_s, unsigned char n_f)
  : ASMbase(n_p,n_s,n_f)
{
  geomB = projB = projB2 = nullptr;
}


ASMstruct::ASMstruct (const ASMstruct& patch, unsigned char n_f)
  : ASMbase(patch,n_f)
{
  geomB = patch.geomB;
  projB = patch.projB;
  projB2 = patch.projB2;
}


ASMstruct::~ASMstruct ()
{
  if (projB && projB != geomB)
    delete projB;

  delete projB2;

  if (geomB && !shareFE)
    delete geomB;
}


bool ASMstruct::separateProjectionBasis () const
{
  return projB && projB != geomB;
}


bool ASMstruct::addXNodes (unsigned short int dim, size_t nXn, IntVec& nodes)
{
  if (dim != ndim-1)
  {
    std::cerr <<" *** ASMstruct::addXNodes: Invalid boundary dimension "<< dim
              <<", only "<< (int)ndim-1 <<" is allowed."<< std::endl;
    return false;
  }
  else if (!geomB || shareFE == 'F')
    return false; // logic error

  else if (MNPC.size() == nel && MLGE.size() == nel)
  {
    // Extend the element number and element topology arrays to double size,
    // to account for extra-ordinary elements associated with contact surfaces
    myMLGE.resize(2*nel,0);
    myMNPC.resize(2*nel);
  }
  else if (MLGE.size() != 2*nel || MNPC.size() != 2*nel)
  {
    // Already added interface elements, currently not allowed
    std::cerr <<" *** ASMstruct::addXNodes: Already have interface elements."
              << std::endl;
    return false;
  }

  // Add nXn extra-ordinary nodes to the list of global node numbers
  for (size_t i = 0; i < nXn; i++)
  {
    if (nodes.size() == i)
      nodes.push_back(++gNod);
    myMLGN.push_back(nodes[i]);
  }

  return true;
}


bool ASMstruct::checkThreadGroups (const std::vector<std::set<int>>& nodes,
                                   int group, bool ignoreGlobalLM)
{
  if (group < 0 || group >= (int)nodes.size())
    return false;

#if SP_DEBUG > 1
  int n = 0;
  std::cout <<"\n\t   nodes:";
  for (int node : nodes[group])
    std::cout << ((++n)%10 ? " " : "\n\t          ") << node;
#endif

  bool ok = true;
  for (int k = 0; k < group; k++)
    for (int node : nodes[group])
      if ((this->getLMType(node+1) != 'G' || !ignoreGlobalLM) &&
          nodes[k].find(node) != nodes[k].end()) {
        std::cout <<"\n  ** Warning: Node "<< node <<" is present on both"
                  <<" thread "<< k+1 <<" and thread "<< group+1;
        ok = false;
      }
  if (!ok) std::cout << std::endl;

  return ok;
}


bool ASMstruct::diracPoint (Integrand& integr, GlobalIntegral& glInt,
                            const double* u, const Vec3& pval)
{
  int iel = this->findElementContaining(u);
  if (iel < 1 || iel > (int)nel)
  {
    std::cerr <<" *** ASMstruct::diracPoint: The point";
    for (unsigned char i = 0; i < ndim; i++) std::cerr <<" "<< u[i];
    std::cerr <<" is outside the patch domain."<< std::endl;
    return false;
  }

#ifdef INDEX_CHECK
  double uElm[6];
  this->getElementBorders(iel,uElm);
  for (unsigned char i = 0; i < ndim; i++)
    if (u[i] < uElm[2*i] || u[i] > uElm[2*i+1])
    {
      unsigned char d;
      std::cerr <<" *** ASMstruct::diracPoint: The point";
      for (d = 0; d < ndim; d++) std::cerr <<" "<< u[d];
      std::cerr <<" is not within the domain of the found element "<< iel
                <<" which is";
      for (d = 0; d < ndim; d++)
        std::cerr <<(d ? "x[":" [") << uElm[2*d] <<","<< uElm[2*d+1] <<"]";
      std::cerr << std::endl;
      return false;
    }
#if SP_DEBUG > 1
    else
    {
      unsigned char d;
      std::cerr <<"   * The point";
      for (d = 0; d < ndim; d++) std::cerr <<" "<< u[d];
      std::cerr <<" is within the domain of element "<< iel <<" which is";
      for (d = 0; d < ndim; d++)
        std::cerr << (d ? "x[":" [") << uElm[2*d] <<","<< uElm[2*d+1] <<"]";
      std::cerr << std::endl;
    }
#endif
#endif

  FiniteElement fe;
  fe.iel = MLGE[iel-1];
  if (ndim > 0) fe.u = u[0];
  if (ndim > 1) fe.v = u[1];
  if (ndim > 2) fe.w = u[2];
  this->evaluateBasis(fe.u,fe.v,fe.w,fe.N);

  LocalIntegral* A = integr.getLocalIntegral(MNPC[iel-1].size(),fe.iel,true);
  bool ok = (integr.evalPoint(*A,fe,pval) &&
             integr.finalizeElement(*A,fe,TimeDomain(),0) &&
             glInt.assemble(A,fe.iel));
  A->destruct();

  return ok;
}


void ASMstruct::swapProjectionBasis ()
{
  if (projB2)
    std::swap(projB, projB2);
}
