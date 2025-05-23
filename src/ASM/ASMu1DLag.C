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
#include "Utilities.h"
#include "IFEM.h"
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

  bool ok = true;
  if (myMLGN.empty())
  {
    myMLGN.resize(nnod);
    std::iota(myMLGN.begin(),myMLGN.end(),gNod+1);
  }
  else
    ok = myMLGN.size() == nnod && gNod == 0;

  if (myMLGE.empty())
  {
    myMLGE.resize(nel);
    std::iota(myMLGE.begin(),myMLGE.end(),gEl+1);
  }
  else if (ok)
    ok = myMLGE.size() == nel && gEl == 0;

  if (!ok)
  {
    std::cerr <<" *** ASMu1DLag::generateOrientedFEModel: Array mismatch, "
              <<" size(coord)="<< myCoord.size() <<" size(MLGN)="<< MLGN.size()
              <<" size(MNPC)="<< MNPC.size() <<" size(MLGE)="<< MLGE.size()
              << std::endl;
    return false;
  }

  if (nsd == 3 && nf == 6)
  {
    // This is a 3D beam problem, allocate the nodal/element rotation tensors.
    // The nodal rotations are updated during the simulation according to the
    // deformation state, whereas the element tensors are kept constant.
    myCS.resize(nel,Tensor(3));
    myT.resize(nnod,Tensor(3,true)); // Initialize nodal rotations to unity
  }

  gNod += nnod;
  gEl  += nel;

  return myCS.empty() ? true : this->initLocalElementAxes(Zaxis);
}


int ASMu1DLag::getElementSetIdx (const std::string& setName) const
{
  int iset = 1;
  for (const ASM::NodeSet& es : elemSets)
    if (es.first == setName)
      return iset;
    else
      ++iset;

  return 0;
}


const IntVec& ASMu1DLag::getElementSet (int iset) const
{
  if (iset > 0 && iset <= static_cast<int>(elemSets.size()))
    return elemSets[iset-1].second;

  return this->ASMbase::getElementSet(iset);
}


/*!
  If \a iel is negative, the absolute value is taken as the external element ID.
  Otherwise, it is taken as the 1-based internal element index within the patch.
*/

bool ASMu1DLag::isInElementSet (int iset, int iel) const
{
  if (iset < 1 || iset > static_cast<int>(elemSets.size()))
    return false;

  if (iel < 0)
    iel = this->getElmIndex(-iel);

  return utl::findIndex(elemSets[iset-1].second,iel) >= 0;
}


int ASMu1DLag::parseElemSet (const std::string& setName, const char* cset)
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
  {
    for (int j : MNPC[i])
      grid.setNode(k++,j);
    grid.setElmId(1+i,MLGE[i]);
  }

  return true;
}
