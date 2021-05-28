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


IntVec& ASMu2DLag::getNodeSet (const std::string& setName)
{
  for (ASM::NodeSet& ns : nodeSets)
    if (ns.first == setName)
      return ns.second;

  nodeSets.push_back(std::make_pair(setName,IntVec()));
  return nodeSets.back().second;
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
