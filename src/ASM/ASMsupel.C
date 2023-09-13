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
#include <numeric>


bool ASMsupel::read (std::istream& is)
{
  // Lambda function for reading the supernode coordinates.
  auto&& readCoord = [&is](Vec3Vec& Xsup)
  {
    int numNod = 0;
    is >> numNod;
    Xsup.resize(numNod);
    for (Vec3& Xn : Xsup) is >> Xn;
  };

  myElmMat.resize(1,1);

  char c;
  int readMat = 0;
  while (readMat < 7 && is.get(c))
    switch (c) {
    case 'K':
    case 'k':
      is >> myElmMat.A.front();
      if (c == 'K') // assume stored column-wise
        myElmMat.A.front().transpose();
      readMat |= 1;
      break;
    case 'R':
      is >> myElmMat.b.front();
      readMat |= 2;
      break;
    case 'L':
      {
        // The load vector is stored as a ndof x 1 matrix and not a vector
        Matrix tmpMat;
        is >> tmpMat;
        myElmMat.b.front() = tmpMat.getColumn(1);
      }
      readMat |= 2;
      break;
    case 'G':
      readCoord(myNodes);
      readMat |= 4;
      break;
    case '\n':
      break;
    default:
      is.putback(c);
      if (readMat)
      {
        std::cerr <<" *** ASMsupel::read: Unknown label "<< c << std::endl;
        return false;
      }
      else
      {
        // Assuming the order G, K, L but without the labels
        readCoord(myNodes);
        is >> myElmMat.A.front() >> myElmMat.b.front();
        readMat = 7;
      }
    }

  return readMat == 7 && is.good();
}


bool ASMsupel::write (std::ostream& os, int) const
{
  // Write out the spider as a lagrangian mesh
  os <<"# LAGRANGIAN nodes="<< 1+nnod <<" elements="<< nnod
     <<" type=superelement\n";
  os << this->getCoord(0) <<"\n";
  for (const Vec3& X : myNodes)
    os << X <<"\n";
  for (size_t i = 1; i <= nnod; i++)
    os <<"0 "<< i <<"\n";

  return os.good();
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

  if (nodeSets.size() == 1 && nodeSets.front().second.empty())
  {
    // Define the set of supernodes for this patch
    nodeSets.front().second.resize(nnod);
    std::iota(nodeSets.front().second.begin(),nodeSets.front().second.end(),1);
  }

  return true;
}


int ASMsupel::getNodeSetIdx (const std::string& setName) const
{
  int idx = 1;
  for (const ASM::NodeSet& ns : nodeSets)
    if (ns.first == setName)
      return idx;
    else
      ++idx;

  return 0;
}


const IntVec& ASMsupel::getNodeSet (int idx) const
{
  int count = 0;
  for (const ASM::NodeSet& ns : nodeSets)
    if (++count == idx)
      return ns.second;

  return this->ASMbase::getNodeSet(idx);
}


IntVec& ASMsupel::getNodeSet (const std::string& setName, int& idx)
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
  Vec3 Xn;
  if (inod == 0)
  {
    // Calculate patch center
    for (const Vec3& X : myNodes)
      Xn += X;
    Xn /= nnod;
  }
  else if (inod <= nnod)
    Xn = myNodes[inod-1];

  return Xn;
}


void ASMsupel::getNodalCoordinates (Matrix& X, bool) const
{
  X.resize(3,nnod);
  for (size_t inod = 1; inod <= nnod; inod++)
  {
    X(1,inod) = myNodes[inod-1].x;
    X(2,inod) = myNodes[inod-1].y;
    X(3,inod) = myNodes[inod-1].z;
  }
}


bool ASMsupel::getElementCoordinates (Matrix& X, int iel, bool) const
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
  grid.unStructResize(nnod,1+nnod);
  grid.setCoor(0,this->getCoord(0));
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



bool ASMsupel::transform (const Matrix& Tlg)
{
  if (Tlg.rows() < 3 || Tlg.cols() < 3)
    return true; // No transformation defined

  // Transform nodal coordinates to global system
  for (Vec3& X : myNodes)
    X = Tlg*X;

#ifdef SP_DEBUG
  std::cout <<"\nGlobal coordinates for superelement "<< idx+1;
  for (const Vec3& X : myNodes)
    std::cout <<"\n" << X;
  std::cout << std::endl;
#endif

  // Transform the element matrices to global system
  for (Matrix& A : myElmMat.A)
    for (size_t k = 1; k < A.cols(); k += 3)
      if (!utl::transform(A,Tlg,k))
        return false;

  // Transform the element force vectors to global system
  for (Vector& b : myElmMat.b)
    for (size_t k = 1; k < b.size(); k += 3)
      if (!utl::transform(b,Tlg,k))
        return false;

  return true;
}
