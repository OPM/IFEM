// $Id$
//==============================================================================
//!
//! \file ASMu2DLag.C
//!
//! \date Aug 12 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of unstructured 2D %Lagrange FE models.
//!
//==============================================================================

#include "ASMu2DLag.h"
#include "ElementBlock.h"
#include "Vec3Oper.h"
#include "IFEM.h"
#include <numeric>
#include <sstream>


ASMu2DLag::ASMu2DLag (unsigned char n_s,
                      unsigned char n_f, char fType) : ASMs2DLag(n_s,n_f)
{
  fileType = fType;
  swapNode34 = false;
}


ASMu2DLag::ASMu2DLag (const ASMu2DLag& p, unsigned char n_f) :
  ASMs2DLag(p,n_f), nodeSets(p.nodeSets)
{
  fileType = 0;
  swapNode34 = p.swapNode34;
}


ASMu2DLag::ASMu2DLag (const ASMu2DLag& p) :
  ASMs2DLag(p), nodeSets(p.nodeSets)
{
  fileType = 0;
  swapNode34 = p.swapNode34;
}


bool ASMu2DLag::read (std::istream& is)
{
  switch (fileType) {
  case 'm':
  case 'M':
    swapNode34 = true;
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

  bool ok = true;
  if (myMLGN.empty())
  {
    myMLGN.resize(nnod);
    std::iota(myMLGN.begin(),myMLGN.end(),gNod+1);
  }
  else
    ok = myMLGN.size() == nnod;

  if (myMLGE.empty())
  {
    myMLGE.resize(nel);
    std::iota(myMLGE.begin(),myMLGE.end(),gEl+1);
  }
  else if (ok)
    ok = myMLGE.size() == nel;

  if (!ok)
    std::cerr <<" *** ASMu2DLag::generateFEMTopology: Array mismatch, "
              <<" size(coord)="<< myCoord.size() <<" size(MLGN)="<< MLGN.size()
              <<" size(MNPC)="<< MNPC.size() <<" size(MLGE)="<< MLGE.size()
              << std::endl;

  gNod += nnod;
  gEl  += nel;

  myCache.emplace_back(std::make_unique<BasisFunctionCache>(*this));

  return ok;
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

  nodeSets.emplace_back(setName,IntVec());
  return nodeSets.back().second;
}


int ASMu2DLag::parseNodeBox (const std::string& setName, const char* data)
{
  if (myCoord.empty()) return 0; // No nodes yet

  Vec3 X0, X1;
  std::istringstream(data) >> X0 >> X1;

  // Lambda function for checking if a point is within the bounding box
  auto&& isInside=[&X0,&X1](const Vec3& X)
  {
    for (int i = 0; i < 3; i++)
      if (X[i] < X0[i] || X[i] > X1[i])
        return false;
    return true;
  };

  IntVec nodes;
  for (size_t inod = 0; inod < myCoord.size(); inod++)
    if (isInside(myCoord[inod]))
      nodes.push_back(1+inod);

  IFEM::cout <<"\tBounding Box: "<< X0 <<" - "<< X1
             <<": "<< nodes.size() <<" nodes"<< std::endl;

  if (nodes.empty()) return 0; // No nodes are within the given box

  size_t idx = 0;
  while (idx < nodeSets.size())
    if (nodeSets[idx].first == setName)
    {
      nodeSets[idx].second.insert(nodeSets[idx].second.end(),
                                  nodes.begin(),nodes.end());
      return idx+1;
    }

  nodeSets.push_back(std::make_pair(setName,nodes));
  return nodeSets.size();
}


int ASMu2DLag::getElementSetIdx (const std::string& setName) const
{
  int idx = 1;
  for (const ASM::NodeSet& es : elemSets)
    if (es.first == setName)
      return idx;
    else
      ++idx;

  return 0;
}


const IntVec& ASMu2DLag::getElementSet (int idx) const
{
  int count = 0;
  for (const ASM::NodeSet& es : elemSets)
    if (++count == idx)
      return es.second;

  return this->ASMbase::getElementSet(idx);
}


IntVec& ASMu2DLag::getElementSet (const std::string& setName, int& idx)
{
  idx = 1;
  for (ASM::NodeSet& es : nodeSets)
    if (es.first == setName)
      return es.second;
    else if (idx)
      ++idx;

  elemSets.push_back(std::make_pair(setName,IntVec()));
  return elemSets.back().second;
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
  size_t nmnpc = 0;
  for (const IntVec& mnpc : MNPC)
    nmnpc += mnpc.size();
  grid.unStructResize(nel,nnod,nmnpc);

  size_t i, j, k;
  for (i = 0; i < nnod; i++)
    grid.setCoor(i,this->getCoord(1+i));

  for (i = k = 0; i < nel; i++)
  {
    for (j = 0; j < MNPC[i].size(); j++)
      if (j > 1 && swapNode34 && MNPC[i].size() == 4)
        grid.setNode(k++,MNPC[i][5-j]);
      else
        grid.setNode(k++,MNPC[i][j]);
    grid.endOfElm(k);
  }
  return true;
}
