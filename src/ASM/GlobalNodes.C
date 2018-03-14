// $Id$
//==============================================================================
//!
//! \file GlobalNodes.C
//!
//! \date Mar 13 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Simple global node establishment for unstructured FE models.
//!
//==============================================================================

#include "GlobalNodes.h"
#include "ASMunstruct.h"
#include "Utilities.h"


GlobalNodes::IntVec GlobalNodes::getBoundaryNodes(const LR::LRSpline& lr,
                                                  int dim, int lidx, int orient)
{
  LR::parameterEdge edge;
  if (dim == 0) {
    if (lr.nVariate() == 2) {
      switch (lidx) {
        case 1: edge = LR::WEST | LR::SOUTH; break;
        case 2: edge = LR::EAST | LR::SOUTH; break;
        case 3: edge = LR::WEST | LR::NORTH; break;
        case 4: edge = LR::EAST | LR::NORTH; break;
      }
    } else {
      switch (lidx) {
        case 1: edge = LR::WEST | LR::SOUTH | LR::BOTTOM; break;
        case 2: edge = LR::EAST | LR::SOUTH | LR::BOTTOM; break;
        case 3: edge = LR::WEST | LR::NORTH | LR::BOTTOM; break;
        case 4: edge = LR::EAST | LR::NORTH | LR::BOTTOM; break;
        case 5: edge = LR::WEST | LR::SOUTH | LR::TOP; break;
        case 6: edge = LR::EAST | LR::SOUTH | LR::TOP; break;
        case 7: edge = LR::WEST | LR::NORTH | LR::TOP; break;
        case 8: edge = LR::EAST | LR::NORTH | LR::TOP; break;
      }
    }
  } else if (dim == 1) {
    if (lr.nVariate() == 2) {
      switch (lidx) {
        case 1: edge = LR::WEST; break;
        case 2: edge = LR::EAST; break;
        case 3: edge = LR::SOUTH; break;
        case 4: edge = LR::NORTH; break;
        default: break;
      }
    } else {
      switch (lidx) {
        case 1: edge = LR::BOTTOM | LR::SOUTH; break;
        case 2: edge = LR::BOTTOM | LR::NORTH; break;
        case 3: edge = LR::TOP | LR::SOUTH; break;
        case 4: edge = LR::TOP | LR::NORTH; break;
        case 5: edge = LR::BOTTOM | LR::WEST; break;
        case 6: edge = LR::BOTTOM | LR::EAST; break;
        case 7: edge = LR::TOP | LR::WEST; break;
        case 8: edge = LR::TOP | LR::WEST; break;
        case 9: edge = LR::SOUTH | LR::WEST; break;
        case 10: edge = LR::SOUTH | LR::EAST; break;
        case 11: edge = LR::NORTH | LR::WEST; break;
        case 12: edge = LR::NORTH | LR::EAST; break;
      }
    }
  } else if (dim == 2) {
    switch (lidx) {
      case 1: edge = LR::WEST; break;
      case 2: edge = LR::EAST; break;
      case 3: edge = LR::SOUTH; break;
      case 4: edge = LR::NORTH; break;
      case 5: edge = LR::BOTTOM; break;
      case 6: edge = LR::TOP; break;
    }
  }

  std::vector<LR::Basisfunction*> edgeFunctions;
  lr.getEdgeFunctions(edgeFunctions, edge);

  if (dim == 1) {
    if (lr.nVariate() == 2) {
      int v = (lidx == 1 || lidx == 2) ? 0 : 1;
      int u = 1-v;
      ASMunstruct::Sort(u, v, orient, edgeFunctions);
    } else {
      int dir = (lidx-1)/4;
      int u = dir == 0;
      int v = 1 + (dir != 2);
      ASMunstruct::Sort(u, v, orient, edgeFunctions);
    }
  } else if (dim == 2) {
    int dir = (lidx-1)/2;
    int u = dir == 0;
    int v = 1 + (dir != 2);
    ASMunstruct::Sort(u, v, orient, edgeFunctions);
  }

  GlobalNodes::IntVec lNodes;
  lNodes.reserve(edgeFunctions.size());
  for (const LR::Basisfunction* func : edgeFunctions)
    lNodes.push_back(func->getId());

  return lNodes;
}


/*!
  \brief Class for ordering interfaces in processing order.
*/

class InterfaceOrder {
public:
  //! \brief Comparison operator for interfaces
  //! \param A First interface
  //! \param B Second interface
  bool operator()(const ASM::Interface& A, const ASM::Interface& B) const
  {
    if (A.master != B.master)
      return A.master < B.master;

    if (A.slave != B.slave)
      return A.slave < B.slave;

    if (A.dim != B.dim)
      return A.dim < B.dim;

    return A.midx < B.midx;
  }
};


std::vector<GlobalNodes::IntVec>
GlobalNodes::calcGlobalNodes(const GlobalNodes::LRSplineVec& pchs,
                             const GlobalNodes::InterfaceVec& interfaces)
{
  // count total number of nodes
  size_t nNodes = 0;
  std::vector<GlobalNodes::IntVec> result(pchs.size());
  auto it = result.begin();
  for (const LR::LRSpline* pch : pchs) {
    it->resize(pch->nBasisFunctions());
    std::iota(it->begin(), it->end(), nNodes);
    nNodes += pch->nBasisFunctions();
    ++it;
  }

  // remap common nodes
  InterfaceOrder ifOrder;
  std::set<ASM::Interface, InterfaceOrder> ifset(ifOrder);
  for (const ASM::Interface& it : interfaces)
    ifset.insert(it);
  for (size_t i = 0; i < pchs.size(); ++i) {
    std::map<int,int> old2new;
    for (const ASM::Interface& it : ifset) {
      if (it.master != (int)i+1)
        continue;

      IntVec mNodes = getBoundaryNodes(*pchs[i], it.dim, it.midx, 0);
      IntVec sNodes = getBoundaryNodes(*pchs[it.slave-1], it.dim, it.sidx, it.orient);
      for (size_t n = 0; n < mNodes.size(); ++n)
        old2new[result[it.slave-1][sNodes[n]]] = result[i][mNodes[n]];
    }

    // renumber
    for (size_t j = i; j < pchs.size(); ++j)
      for (int& it : result[j])
        utl::renumber(it, old2new, false);

    // compress
    int maxNode = *std::max_element(result[i].begin(), result[i].end());
    for (size_t j = i+1; j < pchs.size(); ++j)
      for (int& n : result[j])
        if (n > maxNode)
          n = ++maxNode;
  }

  return result;
}
