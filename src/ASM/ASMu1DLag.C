// $Id$
//==============================================================================
//!
//! \file ASMu1DLag.C
//!
//! \date Aug 26 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of unstructured 1D Lagrange FE models.
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


ASMu1DLag::ASMu1DLag (const ASMu1DLag& p, unsigned char n_f) :
  ASMs1DLag(p,n_f), nodeSets(p.nodeSets)
{
  fileType = 0;
}


ASMu1DLag::ASMu1DLag (const ASMu1DLag& p) :
  ASMs1DLag(p), nodeSets(p.nodeSets)
{
  fileType = 0;
}


bool ASMu1DLag::read (std::istream& is)
{
  switch (fileType) {
  case 'm':
  case 'M':
    return ASM::readMatlab(is,myMNPC,myCoord,nodeSets);
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


int ASMu1DLag::getNodeSetIdx (const std::string& setName) const
{
  int idx = 1;
  for (const ASM::NodeSet& ns : nodeSets)
    if (ns.first == setName)
      return idx;
    else
      ++idx;

  return 0;
}


const IntVec& ASMu1DLag::getNodeSet (int idx) const
{
  int count = 0;
  for (const ASM::NodeSet& ns : nodeSets)
    if (++count == idx)
      return ns.second;

  return this->ASMbase::getNodeSet(idx);
}


IntVec& ASMu1DLag::getNodeSet (const std::string& setName)
{
  for (ASM::NodeSet& ns : nodeSets)
    if (ns.first == setName)
      return ns.second;

  nodeSets.push_back(std::make_pair(setName,IntVec()));
  return nodeSets.back().second;
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
