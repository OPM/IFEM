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
  geo = nullptr;
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


void ASMstruct::resetNumbering (int n)
{
  gEl = 0;
  gNod = n;
  xNode.clear();
}


bool ASMstruct::addXNodes (unsigned short int dim, size_t nXn, IntVec& nodes)
{
  if (dim != ndim-1)
  {
    std::cerr <<" *** ASMstruct::addXNodes: Invalid boundary dimension "<< dim
              <<", only "<< (int)ndim-1 <<" is allowed."<< std::endl;
    return false;
  }
  else if (!geo || shareFE == 'F')
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
#if SP_DEBUG > 1
  auto nit = nodes[group].begin();
  std::cout <<"\n\t   nodes: "<< *(nit++);
  for (int k = 1; nit != nodes[group].end(); ++nit, k++)
    std::cout << (k%10 > 0 ? " " : "\n\t          ") << *nit;
#endif

  bool ok = true;
  for (int k = 0; k < group; k++)
    for (const int& node : nodes[group])
      if ((this->getLMType(node+1) != 'G' || !ignoreGlobalLM) &&
          nodes[k].find(node) != nodes[k].end()) {
        std::cout <<"\n  ** Warning: Node "<< node <<" is present on both"
                  <<" thread "<< k+1 <<" and thread "<< group+1;
        ok = false;
      }

  return ok;
}
