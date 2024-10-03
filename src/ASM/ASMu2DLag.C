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
#include "Utilities.h"
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
  if (idx > 0 && idx <= static_cast<int>(nodeSets.size()))
    return nodeSets[idx-1].second;

  return this->ASMbase::getNodeSet(idx);
}


int ASMu2DLag::parseNodeSet (const std::string& setName, const char* cset)
{
  int idx = this->getNodeSetIdx(setName)-1;
  if (idx < 0)
  {
    idx = nodeSets.size();
    nodeSets.emplace_back(setName,IntVec());
  }

  IntVec& mySet = nodeSets[idx].second;
  size_t ifirst = mySet.size();
  utl::parseIntegers(mySet,cset);

  int inod; // Transform to internal node indices
  for (size_t i = ifirst; i < mySet.size(); i++)
    if ((inod = this->getNodeIndex(mySet[i])) > 0)
      mySet[i] = inod;
    else
      IFEM::cout <<"  ** Warning: Non-existing node "<< mySet[i]
                 <<" in node set \""<< setName <<"\""<< std::endl;

  return 1+idx;
}


int ASMu2DLag::parseNodeBox (const std::string& setName, const char* data)
{
  if (myCoord.empty())
    return 0; // No nodes yet

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

  if (nodes.empty())
    return 0; // No nodes are within the given box

  for (size_t idx = 0; idx < nodeSets.size(); idx++)
    if (nodeSets[idx].first == setName)
    {
      if (nodeSets[idx].second.empty())
        nodeSets[idx].second.swap(nodes);
      else for (int inod : nodes)
        if (utl::findIndex(nodeSets[idx].second,inod) < 0)
          nodeSets[idx].second.push_back(inod);
      return idx+1;
    }

  nodeSets.emplace_back(setName,nodes);
  return nodeSets.size();
}


void ASMu2DLag::addToNodeSet (const std::string& setName, int inod)
{
  int idx = this->getNodeSetIdx(setName);
  if (idx < 1)
    nodeSets.emplace_back(setName,IntVec{inod});
  else
    nodeSets[idx-1].second.push_back(inod);
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
  if (idx > 0 && idx <= static_cast<int>(elemSets.size()))
    return elemSets[idx-1].second;

  return this->ASMbase::getElementSet(idx);
}


bool ASMu2DLag::isInElementSet (int idx, int iel) const
{
  if (idx < 1 || idx > static_cast<int>(elemSets.size()))
    return false;

  return utl::findIndex(elemSets[idx-1].second,iel) >= 0;
}


int ASMu2DLag::parseElemSet (const std::string& setName, const char* cset)
{
  int idx = this->getElementSetIdx(setName)-1;
  if (idx < 0)
  {
    idx = elemSets.size();
    elemSets.emplace_back(setName,IntVec());
  }

  IntVec& mySet = elemSets[idx].second;
  size_t ifirst = mySet.size();
  utl::parseIntegers(mySet,cset);

  int iel; // Transform to internal element indices
  for (size_t i = ifirst; i < mySet.size(); i++)
    if ((iel = this->getElmIndex(mySet[i])) > 0)
      mySet[i] = iel;
    else
      IFEM::cout <<"  ** Warning: Non-existing element "<< mySet[i]
                 <<" in element set \""<< setName <<"\""<< std::endl;

  return 1+idx;
}


int ASMu2DLag::parseElemBox (const std::string& setName, const char* data)
{
  Vec3 X0, X1;
  std::istringstream(data) >> X0 >> X1;

  // Lambda function for checking if an element is within the bounding box
  auto&& isInside=[this,&X0,&X1](size_t iel)
  {
    double nelnod = MNPC[iel].size();
    for (size_t j = 0; j < nsd; j++)
    {
      double X = 0.0;
      for (int inod : MNPC[iel])
        X += coord[inod][j];
      X /= nelnod;
      if (X < X0[j] || X > X1[j])
        return false;
    }
    return true;
  };

  Matrix Xnod;
  IntVec elems;
  for (size_t iel = 0; iel < nel; iel++)
    if (isInside(iel))
      elems.push_back(1+iel);

  IFEM::cout <<"\tBounding Box: "<< X0 <<" - "<< X1
             <<": "<< elems.size() <<" elements"<< std::endl;

  if (elems.empty())
    return 0; // No elements are within the given box

  for (size_t idx = 0; idx < elemSets.size(); idx++)
    if (elemSets[idx].first == setName)
    {
      if (elemSets[idx].second.empty())
        elemSets[idx].second.swap(elems);
      else for (int iel : elems)
        if (utl::findIndex(elemSets[idx].second,iel) < 0)
          elemSets[idx].second.push_back(iel);
      return idx+1;
    }

  elemSets.emplace_back(setName,elems);
  return elemSets.size();
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
