// $Id$
//==============================================================================
//!
//! \file ASMsupel.C
//!
//! \date Mar 30 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of general superelements.
//!
//==============================================================================

#include "ASMsupel.h"
#include "GlobalIntegral.h"
#include "ElementBlock.h"
#include "Vec3Oper.h"


bool ASMsupel::read (std::istream& is)
{
  int numNod = 0;
  is >> numNod;
  myNodes.resize(numNod);
  myElmMat.resize(1,1);
  for (Vec3& Xn : myNodes) is >> Xn;
  is >> myElmMat.A.front() >> myElmMat.b.front();
  return is.good();
}


bool ASMsupel::generateFEMTopology ()
{
  nnod = myNodes.size();
  nel  = 1;

  myMLGE = { ++gEl };
  myMLGN.resize(nnod,0);
  myMNPC.resize(1);
  myMNPC.front().resize(nnod,0);
  std::iota(myMLGN.begin(),myMLGN.end(),gNod+1);
  std::iota(myMNPC.front().begin(),myMNPC.front().end(),0);
  gNod += nnod;

  return true;
}


int ASMsupel::getNodeSetIdx (const std::string& setName) const
{
  int idx = 1;
  for (const NodeSet& ns : nodeSets)
    if (ns.first == setName)
      return idx;
    else
      ++idx;

  return 0;
}


const IntVec& ASMsupel::getNodeSet (int idx) const
{
  int count = 0;
  for (const NodeSet& ns : nodeSets)
    if (++count == idx)
      return ns.second;

  return this->ASMbase::getNodeSet(idx);
}


IntVec& ASMsupel::getNodeSet (const std::string& setName)
{
  for (NodeSet& ns : nodeSets)
    if (ns.first == setName)
      return ns.second;

  nodeSets.push_back(std::make_pair(setName,IntVec()));
  return nodeSets.back().second;
}


void ASMsupel::getBoundaryNodes (int lIndex, IntVec& nodes,
                                 int, int, int, bool local) const
{
  nodes = this->getNodeSet(lIndex);
  if (!local)
    for (int& node : nodes)
      node = this->getNodeID(node);
}


Vec3 ASMsupel::getCoord (size_t inod) const
{
  if (inod > 0 && inod <= nnod)
    return myNodes[inod-1];

  return Vec3();
}


void ASMsupel::getNodalCoordinates (Matrix& X) const
{
  X.resize(3,nnod);
  for (size_t inod = 1; inod <= nnod; inod++)
  {
    X(1,inod) = myNodes[inod-1].x;
    X(2,inod) = myNodes[inod-1].y;
    X(3,inod) = myNodes[inod-1].z;
  }
}


bool ASMsupel::getElementCoordinates (Matrix& X, int iel) const
{
  if (iel != 1)
    return false;

  this->getNodalCoordinates(X);
  return true;
}


bool ASMsupel::integrate (Integrand& integrand, GlobalIntegral& glbInt,
                          const TimeDomain&)
{
  return glbInt.assemble(&myElmMat,MLGE.front());
}


/*!
  This method creates a simplified visualization of the superelement
  as lines extending from the centroid to each of the supernodes.
*/

bool ASMsupel::tesselate (ElementBlock& grid, const int*) const
{
  Vec3 XC;
  for (const Vec3& X : myNodes)
    XC += X;
  XC /= nnod;

  grid.unStructResize(nnod,1+nnod);
  grid.setCoor(0,XC);
  for (size_t i = 1; i <= nnod; i++)
  {
    grid.setCoor(i,myNodes[i-1]);
    grid.setNode(2*i-2,0);
    grid.setNode(2*i-1,i);
    grid.setElmId(i,i);
  }

  return true;
}


bool ASMsupel::evalSolution (Matrix& sField, const Vector& locSol,
                             const int*, int) const
{
  size_t i, j, k, nComp = locSol.size() / nnod;
  sField.resize(nComp,1+nnod,true);
  for (j = k = 1; j <= nnod; j++)
    for (i = 1; i <= nComp; i++)
    {
      sField(i,1+j) = locSol(k++);
      sField(i,1) += sField(i,1+j) / (double)nnod;
    }

  return true;
}
