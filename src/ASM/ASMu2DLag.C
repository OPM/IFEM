// $Id$
//==============================================================================
//!
//! \file ASMu2DLag.C
//!
//! \date Aug 12 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of unstructured 2D Lagrange FE models.
//!
//==============================================================================

#include "ASMu2DLag.h"
#include "ElementBlock.h"
#include "GaussQuadrature.h"
#include "Integrand.h"
#include "Lagrange.h"
#include <numeric>


ASMu2DLag::ASMu2DLag (unsigned char n_s,
                      unsigned char n_f, char fType) : ASMs2DLag(n_s,n_f)
{
  fileType = fType;
}


ASMu2DLag::ASMu2DLag (const ASMu2DLag& p, unsigned char n_f) :
  ASMs2DLag(p,n_f), nodeSets(p.nodeSets)
{
  fileType = 0;
}


ASMu2DLag::ASMu2DLag (const ASMu2DLag& p) :
  ASMs2DLag(p), nodeSets(p.nodeSets)
{
  fileType = 0;
}


bool ASMu2DLag::read (std::istream& is)
{
  switch (fileType) {
  case 'm':
  case 'M':
    return ASM::readMatlab(is,myMNPC,myCoord,nodeSets);
  case 'x':
  case 'X':
    return ASM::readXML(is,myMNPC,myCoord,nodeSets);
  default:
    std::cerr <<" *** ASMu2DLag::read: Undefined file format."<< std::endl;
    return false;
  }
}


bool ASMu2DLag::generateFEMTopology ()
{
  p1 = p2 = 2; // So far only linear elements supported

  nnod = myCoord.size();
  nel  = myMNPC.size();

  myMLGN.resize(nnod);
  myMLGE.resize(nel);

  std::iota(myMLGN.begin(),myMLGN.end(),gNod+1);
  std::iota(myMLGE.begin(),myMLGE.end(),gEl+1);

  gNod += nnod;
  gEl  += nel;

  myCache.emplace_back(std::make_unique<BasisFunctionCache>(*this, cachePolicy));

  return true;
}


int ASMu2DLag::getNodeSetIdx (const std::string& setName) const
{
  int idx = 1;
  for (const ASM::NodeSet& ns : nodeSets)
    if (ns.first == setName)
      return idx;
    else
      ++idx;

  return 0;
}


const IntVec& ASMu2DLag::getNodeSet (int idx) const
{
  int count = 0;
  for (const ASM::NodeSet& ns : nodeSets)
    if (++count == idx)
      return ns.second;

  return this->ASMbase::getNodeSet(idx);
}


IntVec& ASMu2DLag::getNodeSet (const std::string& setName, int& idx)
{
  idx = 1;
  for (ASM::NodeSet& ns : nodeSets)
    if (ns.first == setName)
      return ns.second;
    else if (idx)
      ++idx;

  nodeSets.push_back(std::make_pair(setName,IntVec()));
  return nodeSets.back().second;
}


void ASMu2DLag::getBoundaryNodes (int lIndex, IntVec& nodes,
                                  int, int, int, bool local) const
{
  nodes = this->getNodeSet(lIndex);
  if (!local)
    for (int& node : nodes)
      node = this->getNodeID(node);
}


void ASMu2DLag::generateThreadGroups (const Integrand&, bool, bool)
{
  // TODO: Add some coloring scheme later
  threadGroups.oneGroup(nel);
}


bool ASMu2DLag::tesselate (ElementBlock& grid, const int*) const
{
  grid.unStructResize(nel,nnod);

  size_t i, j, k;
  for (i = 0; i < nnod; i++)
    grid.setCoor(i,this->getCoord(1+i));

  for (i = k = 0; i < nel; i++)
    for (j = 0; j < MNPC[i].size(); j++)
      if (j > 1 && MNPC[i].size() == 4)
        grid.setNode(k++,MNPC[i][5-j]);
      else
        grid.setNode(k++,MNPC[i][j]);

  return true;
}


ASMu2DLag::BasisFunctionCache::BasisFunctionCache (const ASMu2DLag& pch,
                                                  ASM::CachePolicy plcy) :
  ::BasisFunctionCache<2>(plcy),
  patch(pch)
{
}


bool ASMu2DLag::BasisFunctionCache::internalInit ()
{
  if (!mainQ->xg[0])
    this->setupQuadrature();

  nTotal = patch.nel*mainQ->ng[0]*mainQ->ng[1];
  if (reducedQ->xg[0])
    nTotalRed = patch.nel*reducedQ->ng[0]*reducedQ->ng[1];

  return true;
}


void ASMu2DLag::BasisFunctionCache::internalCleanup ()
{
  mainQ->reset();
  reducedQ->reset();
}


bool ASMu2DLag::BasisFunctionCache::setupQuadrature ()
{
  // Get Gaussian quadrature points and weights
  for (int d = 0; d < 2; d++)
  {
    mainQ->ng[d] = patch.getNoGaussPt(d == 0 ? patch.p1 : patch.p2);
    mainQ->xg[d] = GaussQuadrature::getCoord(mainQ->ng[d]);
    mainQ->wg[d] = GaussQuadrature::getWeight(mainQ->ng[d]);
    if (!mainQ->xg[d] || !mainQ->wg[d]) return false;
  }

  // Get the reduced integration quadrature points, if needed
  int nRed = integrand ? integrand->getReducedIntegration(mainQ->ng[0]) : 0;
  if (nRed > 0)
  {
    reducedQ->xg[0] = reducedQ->xg[1] = GaussQuadrature::getCoord(nRed);
    reducedQ->wg[0] = reducedQ->wg[1] = GaussQuadrature::getWeight(nRed);
    if (!reducedQ->xg[0] || !reducedQ->wg[0]) return false;
  } else if (nRed < 0)
    nRed = mainQ->ng[0]; // The integrand needs to know nGauss

  reducedQ->ng[0] = reducedQ->ng[1] = nRed;

  return true;
}


BasisFunctionVals ASMu2DLag::BasisFunctionCache::calculatePt (size_t el,
                                                              size_t gp,
                                                              bool reduced) const
{
  std::array<size_t,2> gpIdx = this->gpIndex(gp,reduced);
  const Quadrature& q = reduced ? *reducedQ : *mainQ;

  const ASMu2DLag& pch = static_cast<const ASMu2DLag&>(patch);

  BasisFunctionVals result;
  if (nderiv == 1)
    Lagrange::computeBasis(result.N,result.dNdu,
                           pch.p1,q.xg[0][gpIdx[0]],
                           pch.p2,q.xg[1][gpIdx[1]]);

  return result;
}


void ASMu2DLag::BasisFunctionCache::calculateAll ()
{
  // Evaluate basis function values and derivatives at all integration points.
  // We do this before the integration point loop to exploit multi-threading
  // in the integrand evaluations, which may be the computational bottleneck.
  size_t iel, jp, rp;
  for (iel = jp = rp = 0; iel < patch.nel; iel++)
  {
    for (int j = 0; j < mainQ->ng[1]; j++)
      for (int i = 0; i < mainQ->ng[0]; i++, jp++)
        values[jp] = this->calculatePt(iel,j*mainQ->ng[0]+i,false);

    if (reducedQ->xg[0])
      for (int j = 0; j < reducedQ->ng[1]; j++)
        for (int i = 0; i < reducedQ->ng[0]; i++, rp++)
          valuesRed[rp] = this->calculatePt(iel,j*reducedQ->ng[0]+i,true);
  }
}
