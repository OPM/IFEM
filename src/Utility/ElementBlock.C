// $Id$
//==============================================================================
//!
//! \file ElementBlock.C
//!
//! \date Jun 03 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Storage of a standard FE grid block of uniform element type.
//!
//==============================================================================

#include "ElementBlock.h"
#include "Vec3Oper.h"
#include "Point.h"
#include <algorithm>
#include <numeric>


ElementBlock::ElementBlock (size_t nenod)
{
  if (nenod != 2 && nenod != 3 && nenod != 4 && nenod != 6 && nenod != 9 &&
      nenod != 8 && nenod != 27)
  {
    std::cout <<"ElementBlock: Invalid number of element nodes "<< nenod
              <<" reset to 8"<< std::endl;
    nenod = 8;
  }
  nen = nenod;
}


void ElementBlock::resize (size_t nI, size_t nJ, size_t nK)
{
  coord.resize(nI*nJ*nK);
  param.resize(nI*nJ*nK);
  if (nen == 2 && nJ < 2 && nK < 2)
    MMNPC.resize(2*(nI-1));
  else if (nen == 3 && nK < 2)
    MMNPC.resize(6*(nI-1)*(nJ-1));
  else if (nen == 4 && nK < 2)
    MMNPC.resize(4*(nI-1)*(nJ-1));
  else if (nen == 6 && nK < 2)
    MMNPC.resize(12*(nI-1)*(nJ-1));
  else if (nen == 9 && nK < 2)
    MMNPC.resize(9*(nI-1)*(nJ-1)/4);
  else if (nen == 8)
    MMNPC.resize(8*(nI-1)*(nJ-1)*(nK-1));
  else if (nen == 27)
    MMNPC.resize(27*(nI-1)*(nJ-1)*(nK-1)/8);

  MINEX.resize(MMNPC.size()/nen,0);
  std::iota(MINEX.begin(),MINEX.end(),1);
}


void ElementBlock::unStructResize (size_t nEl, size_t nPts)
{
  coord.resize(nPts);
  param.resize(nPts);
  MMNPC.resize(nen*nEl);
  MINEX.resize(MMNPC.size()/nen,0);
  std::iota(MINEX.begin(),MINEX.end(),1);
}


bool ElementBlock::setCoor (size_t i, size_t j, Real x)
{
  if (i >= coord.size() || j >= 3) return false;

  coord[i][j] = x;
  return true;
}


bool ElementBlock::setCoor (size_t i, const Vec3& X)
{
  if (i >= coord.size()) return false;

  coord[i] = X;
  return true;
}


bool ElementBlock::setParams (size_t i, Real u, Real v, Real w)
{
  if (i >= param.size()) return false;

  param[i][0] = u;
  param[i][1] = v;
  param[i][2] = w;
  return true;
}


bool ElementBlock::setNode (size_t i, int nodeNumb)
{
  if (i >= MMNPC.size()) return false;

  MMNPC[i] = nodeNumb;
  return true;
}


bool ElementBlock::addLine (Real x1, Real y1, Real z1,
                            Real x2, Real y2, Real z2)
{
  if (nen != 2) return false;

  coord.push_back(Vec3(x1,y1,z1));
  coord.push_back(Vec3(x2,y2,z2));
  MMNPC.push_back(coord.size()-2);
  MMNPC.push_back(coord.size()-1);
  MINEX.push_back(MINEX.size()+1);
  return true;
}


size_t ElementBlock::addLine (size_t i1, const Vec3& X2)
{
  if (nen != 2 || i1 >= coord.size()) return 0;

  coord.push_back(X2);
  MMNPC.push_back(i1);
  MMNPC.push_back(coord.size()-1);
  MINEX.push_back(MINEX.size()+1);
  return coord.size()-1;
}


void ElementBlock::merge (const ElementBlock* other,
                          std::vector<int>& nodeNums, bool uniqNodes)
{
  nodeNums.clear();
  nodeNums.reserve(other->coord.size());

  std::vector<Vec3>::const_iterator cit;
  for (const Vec3& X : other->coord)
    if (!uniqNodes || (cit = find(coord.begin(),coord.end(),X)) == coord.end())
    {
      nodeNums.push_back(coord.size());
      coord.push_back(X);
    }
    else
      nodeNums.push_back(cit - coord.begin());

  for (int inod : other->MMNPC)
    MMNPC.push_back(nodeNums[inod]);

  MINEX.insert(MINEX.end(),other->MINEX.begin(),other->MINEX.end());
}


void ElementBlock::merge (const ElementBlock& other, bool uniqNodes)
{
  std::vector<int> nodes;
  this->merge(&other,nodes,uniqNodes);
}


void ElementBlock::transform (const Matrix& Tlg)
{
  if (!Tlg.empty())
    for (Vec3& X : coord)
      X = Tlg * X;
}


utl::Point ElementBlock::getCenter (size_t i) const
{
  if (i < 1 || i > MINEX.size())
    return utl::Point();

  utl::Point XC;
  for (size_t j = 0; j < nen; j++)
  {
    const Prm3& uu = param[MMNPC[nen*(i-1)+j]];
    XC += utl::Point(coord[MMNPC[nen*(i-1)+j]],{uu[0],uu[1],uu[2]});
  }
  XC /= nen;

  return XC;
}


CubeBlock::CubeBlock (const Vec3& X0, double dX) : ElementBlock(8)
{
  dX *= 0.5;
  this->resize(2,2,2);
  this->setCoor(0, X0.x-dX, X0.y-dX, X0.z-dX);
  this->setCoor(1, X0.x+dX, X0.y-dX, X0.z-dX);
  this->setCoor(2, X0.x+dX, X0.y+dX, X0.z-dX);
  this->setCoor(3, X0.x-dX, X0.y+dX, X0.z-dX);
  this->setCoor(4, X0.x-dX, X0.y-dX, X0.z+dX);
  this->setCoor(5, X0.x+dX, X0.y-dX, X0.z+dX);
  this->setCoor(6, X0.x+dX, X0.y+dX, X0.z+dX);
  this->setCoor(7, X0.x-dX, X0.y+dX, X0.z+dX);
  for (int i = 0; i < 8; i++)
    this->setNode(i,i);
  this->setElmId(1,1);
}
