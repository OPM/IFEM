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
  params.resize(nI*nJ*nK);
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
  params.resize(nPts);
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
  if (i >= params.size()) return false;

  params[i] = {{u,v,w}};
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


void ElementBlock::merge (const ElementBlock* other, std::vector<int>& nodeNums)
{
  nodeNums.resize(other->coord.size());

  size_t i;
  std::vector<Vec3>::const_iterator cit;
  for (i = 0; i < nodeNums.size(); i++)
    if ((cit = find(coord.begin(),coord.end(),other->coord[i])) == coord.end())
    {
      nodeNums[i] = coord.size();
      coord.push_back(other->coord[i]);
    }
    else
      nodeNums[i] = cit - coord.begin();

  for (i = 0; i < other->MMNPC.size(); i++)
    MMNPC.push_back(nodeNums[other->MMNPC[i]]);

  MINEX.insert(MINEX.end(),other->MINEX.begin(),other->MINEX.end());
}
