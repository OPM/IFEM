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


ElementBlock::ElementBlock (size_t nenod)
{
  if (nenod != 2 && nenod != 4 && nenod != 9 && nenod != 8 && nenod != 27)
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
  if (nen == 2 && nJ < 2 && nK < 2)
    MMNPC.resize(2*(nI-1));
  else if (nen == 4 && nK < 2)
    MMNPC.resize(4*(nI-1)*(nJ-1));
  else if (nen == 9 && nK < 2)
    MMNPC.resize(9*(nI-1)*(nJ-1)/4);
  else if (nen == 8)
    MMNPC.resize(8*(nI-1)*(nJ-1)*(nK-1));
  else if (nen == 27)
    MMNPC.resize(27*(nI-1)*(nJ-1)*(nK-1)/8);

  MINEX.resize(MMNPC.size()/nen,0);
}


bool ElementBlock::setCoor (size_t i, real x, real y, real z)
{
  if (i >= coord.size()) return false;

  coord[i] = Vec3(x,y,z);
  return true;
}


bool ElementBlock::setCoor (size_t i, size_t j, real x)
{
  if (i >= coord.size() || j >= 3) return false;

  coord[i][j] = x;
  return true;
}


bool ElementBlock::setNode (size_t i, int nodeNumb)
{
  if (i >= MMNPC.size()) return false;

  MMNPC[i] = nodeNumb;
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
