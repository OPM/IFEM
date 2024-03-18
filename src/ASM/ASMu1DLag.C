// $Id$
//==============================================================================
//!
//! \file ASMu1DLag.C
//!
//! \date Aug 26 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of unstructured 1D %Lagrange FE models.
//!
//==============================================================================

#include "ASMu1DLag.h"
#include "ElementBlock.h"
#include <numeric>


ASMu1DLag::ASMu1DLag (unsigned char n_s,
                      unsigned char n_f, char fType) : ASMs1DLag(n_s,n_f)
{
  fileType = fType;
}


ASMu1DLag::ASMu1DLag (const ASMu1DLag& p, unsigned char n_f) : ASMs1DLag(p,n_f)
{
  fileType = 0;
}


ASMu1DLag::ASMu1DLag (const ASMu1DLag& p) : ASMs1DLag(p)
{
  fileType = 0;
}


bool ASMu1DLag::read (std::istream& is)
{
  switch (fileType) {
  case 'm':
  case 'M':
    return ASM::readMatlab(is,myMNPC,myCoord,nodeSets);
  case 'x':
  case 'X':
    return ASM::readXML(is,myMNPC,myCoord,nodeSets,&elemSets);
  default:
    std::cerr <<" *** ASMu1DLag::read: Undefined file format."<< std::endl;
    return false;
  }
}


bool ASMu1DLag::generateOrientedFEModel (const Vec3& Zaxis)
{
  p1 = 2; // So far only linear elements supported

  nnod = myCoord.size();
  nel  = myMNPC.size();

  myMLGN.resize(nnod);
  myMLGE.resize(nel);
  if (nsd == 3 && nf == 6)
  {
    // This is a 3D beam problem, allocate the nodal/element rotation tensors.
    // The nodal rotations are updated during the simulation according to the
    // deformation state, whereas the element tensors are kept constant.
    myCS.resize(nel,Tensor(3));
    myT.resize(nnod,Tensor(3,true)); // Initialize nodal rotations to unity
  }

  std::iota(myMLGN.begin(),myMLGN.end(),gNod+1);
  std::iota(myMLGE.begin(),myMLGE.end(),gEl+1);

  gNod += nnod;
  gEl  += nel;

  return myCS.empty() ? true : this->initLocalElementAxes(Zaxis);
}


int ASMu1DLag::getElementSetIdx (const std::string& setName) const
{
  int idx = 1;
  for (const ASM::NodeSet& es : elemSets)
    if (es.first == setName)
      return idx;
    else
      ++idx;

  return 0;
}


const IntVec& ASMu1DLag::getElementSet (int idx) const
{
  int count = 0;
  for (const ASM::NodeSet& es : elemSets)
    if (++count == idx)
      return es.second;

  return this->ASMbase::getElementSet(idx);
}


IntVec& ASMu1DLag::getElementSet (const std::string& setName, int& idx)
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


void ASMu1DLag::getBoundaryNodes (int lIndex, IntVec& nodes,
                                  int, int, int, bool local) const
{
  nodes = this->getNodeSet(lIndex);
  if (!local)
    for (int& node : nodes)
      node = this->getNodeID(node);
}


void ASMu1DLag::shiftGlobalElmNums (int eshift)
{
  this->ASMs1DLag::shiftGlobalElmNums(eshift);

  for (ASM::NodeSet& es : elemSets)
    for (int& iel : es.second)
      iel += eshift;
}


bool ASMu1DLag::tesselate (ElementBlock& grid, const int*) const
{
  grid.unStructResize(nel,nnod);

  size_t i, k;
  for (i = 0; i < nnod; i++)
    grid.setCoor(i,this->getCoord(1+i));

  for (i = k = 0; i < nel; i++)
    for (int j : MNPC[i])
      grid.setNode(k++,j);

  return true;
}
