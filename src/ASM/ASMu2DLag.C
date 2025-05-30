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
  int iset = 1;
  for (const ASM::NodeSet& ns : nodeSets)
    if (ns.first == setName)
      return iset;
    else
      ++iset;

  return 0;
}


const IntVec& ASMu2DLag::getNodeSet (int iset) const
{
  if (iset > 0 && iset <= static_cast<int>(nodeSets.size()))
    return nodeSets[iset-1].second;

  return this->ASMbase::getNodeSet(iset);
}


bool ASMu2DLag::getNodeSet (int iset, std::string& name) const
{
  if (iset < 0 || iset > static_cast<int>(nodeSets.size()))
    return false;
  else if (nodeSets[iset-1].second.empty())
    return false;

  name = nodeSets[iset-1].first;
  return true;
}


/*!
  If \a inod is negative, the absolute value is taken as the external node ID.
  Otherwise, it is taken as the 1-based internal node index within the patch.
*/

bool ASMu2DLag::isInNodeSet (int iset, int inod) const
{
  if (iset < 1 || iset > static_cast<int>(nodeSets.size()))
    return false;

  if (inod < 0)
    inod = this->getNodeIndex(-inod);

  return utl::findIndex(nodeSets[iset-1].second,inod) >= 0;
}


int ASMu2DLag::parseNodeSet (const std::string& setName, const char* cset)
{
  int iset = this->getNodeSetIdx(setName)-1;
  if (iset < 0)
  {
    iset = nodeSets.size();
    nodeSets.emplace_back(setName,IntVec());
  }

  IntVec& mySet = nodeSets[iset].second;
  size_t ifirst = mySet.size();
  utl::parseIntegers(mySet,cset);

  int inod; // Transform to internal node indices
  for (size_t i = ifirst; i < mySet.size(); i++)
    if ((inod = this->getNodeIndex(mySet[i])) > 0)
      mySet[i] = inod;
    else
      IFEM::cout <<"  ** Warning: Non-existing node "<< mySet[i]
                 <<" in node set \""<< setName <<"\""<< std::endl;

  return 1+iset;
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

  for (size_t iset = 0; iset < nodeSets.size(); iset++)
    if (nodeSets[iset].first == setName)
    {
      if (nodeSets[iset].second.empty())
        nodeSets[iset].second.swap(nodes);
      else for (int inod : nodes)
        if (utl::findIndex(nodeSets[iset].second,inod) < 0)
          nodeSets[iset].second.push_back(inod);
      return iset+1;
    }

  nodeSets.emplace_back(setName,nodes);
  return nodeSets.size();
}


int ASMu2DLag::addToNodeSet (const std::string& setName, int inod, bool extId)
{
  int iset = this->getNodeSetIdx(setName)-1;
  if (iset < 0)
  {
    iset = nodeSets.size();
    nodeSets.emplace_back(setName,IntVec());
  }

  int nId = inod;
  if (extId) // Transform to internal node index
    inod = this->getNodeIndex(inod);

  if (inod > 0 && inod <= static_cast<int>(nnod))
    nodeSets[iset].second.push_back(inod);
  else
  {
    iset = -1;
    IFEM::cout <<"  ** Warning: Non-existing node "<< nId
               <<" in node set \""<< setName <<"\""<< std::endl;
  }

  return iset;
}


int ASMu2DLag::getElementSetIdx (const std::string& setName) const
{
  int iset = 1;
  for (const ASM::NodeSet& es : elemSets)
    if (es.first == setName)
      return iset;
    else
      ++iset;

  return 0;
}


const IntVec& ASMu2DLag::getElementSet (int iset) const
{
  if (iset > 0 && iset <= static_cast<int>(elemSets.size()))
    return elemSets[iset-1].second;

  return this->ASMbase::getElementSet(iset);
}


bool ASMu2DLag::getElementSet (int iset, std::string& name) const
{
  if (iset < 1 || iset > static_cast<int>(elemSets.size()))
    return false;
  else if (elemSets[iset-1].second.empty())
    return false;

  name = elemSets[iset-1].first;
  return true;
}


/*!
  If \a iel is negative, the absolute value is taken as the external element ID.
  Otherwise, it is taken as the 1-based internal element index within the patch.
*/

bool ASMu2DLag::isInElementSet (int iset, int iel) const
{
  if (iset < 1 || iset > static_cast<int>(elemSets.size()))
    return false;

  if (iel < 0)
    iel = this->getElmIndex(-iel);

  return utl::findIndex(elemSets[iset-1].second,iel) >= 0;
}


int ASMu2DLag::parseElemSet (const std::string& setName, const char* cset)
{
  int iset = this->getElementSetIdx(setName)-1;
  if (iset < 0)
  {
    iset = elemSets.size();
    elemSets.emplace_back(setName,IntVec());
  }

  IntVec& mySet = elemSets[iset].second;
  size_t ifirst = mySet.size();
  utl::parseIntegers(mySet,cset);

  int iel; // Transform to internal element indices
  for (size_t i = ifirst; i < mySet.size(); i++)
    if ((iel = this->getElmIndex(mySet[i])) > 0)
      mySet[i] = iel;
    else
      IFEM::cout <<"  ** Warning: Non-existing element "<< mySet[i]
                 <<" in element set \""<< setName <<"\""<< std::endl;

  return 1+iset;
}


int ASMu2DLag::parseElemBox (const std::string& setName,
                             const std::string& unionSet, const char* data)
{
  int iuSet = unionSet.empty() ? 0 : this->getElementSetIdx(unionSet);

  Vec3 X0, X1;
  std::istringstream(data) >> X0 >> X1;

  // Lambda function for checking if an element is within the bounding box
  auto&& isInside=[this,iuSet,&X0,&X1](size_t iel)
  {
    if (iuSet > 0 && !this->isInElementSet(iuSet,1+iel))
      return false; // Filter out all elements not in the set unionSet

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

  for (size_t iset = 0; iset < elemSets.size(); iset++)
    if (elemSets[iset].first == setName)
    {
      if (elemSets[iset].second.empty())
        elemSets[iset].second.swap(elems);
      else for (int iel : elems)
        if (utl::findIndex(elemSets[iset].second,iel) < 0)
          elemSets[iset].second.push_back(iel);
      return iset+1;
    }

  elemSets.emplace_back(setName,elems);
  return elemSets.size();
}


int ASMu2DLag::addToElemSet (const std::string& setName, int iel, bool extId)
{
  int iset = this->getElementSetIdx(setName)-1;
  if (iset < 0)
  {
    iset = elemSets.size();
    elemSets.emplace_back(setName,IntVec());
  }

  int eId = iel;
  if (extId) // Transform to internal element index
    iel = this->getElmIndex(iel);

  if (iel > 0 && iel <= static_cast<int>(nel))
    elemSets[iset].second.push_back(iel);
  else
  {
    iset = -1;
    IFEM::cout <<"  ** Warning: Non-existing element "<< eId
               <<" in element set \""<< setName <<"\""<< std::endl;
  }

  return iset;
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
  size_t nelms = nel;
  for (const IntVec& mnpc : MNPC)
    if (mnpc.size() > 1) // ignore 1-noded elements (point masses, etc.)
      nmnpc += mnpc.size();
    else
      --nelms;
  grid.unStructResize(nelms,nnod,nmnpc);

  size_t i, j, k, e;
  for (i = 0; i < nnod; i++)
    grid.setCoor(i,this->getCoord(1+i));

  for (i = k = e = 0; i < nel; i++)
    if (MNPC[i].size() > 1) // ignore 1-noded elements
    {
      for (j = 0; j < MNPC[i].size(); j++)
        if (j > 1 && swapNode34 && MNPC[i].size() == 4)
          grid.setNode(k++,MNPC[i][5-j]);
        else
          grid.setNode(k++,MNPC[i][j]);
      grid.endOfElm(k);
      grid.setElmId(++e,MLGE[i]);
    }

  return true;
}
