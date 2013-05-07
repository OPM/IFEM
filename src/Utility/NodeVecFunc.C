// $Id$
//==============================================================================
//!
//! \file NodeVecFunc.C
//!
//! \date Apr 2 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Vector function wrapper for a nodal field.
//!
//==============================================================================

#include "NodeVecFunc.h"
#include "SIMbase.h"
#include "ASMbase.h"
#include "Vec3.h"


Vec3 NodeVecFunc::evaluate (const Vec3& xp) const
{
  int idx = this->getPointIndex(xp);
  if (idx < 0) return Vec3();

  return Vec3(&value.front()+model.getNoFields()*idx);
}


int NodeVecFunc::getPointIndex (const Vec3& xp) const
{
  // Search among the earlier nodes found
  std::map<Vec3,int>::const_iterator it = ptMap.find(xp);
  if (it != ptMap.end()) return it->second;

  // Not found, search among all nodes in the model
  const ASMbase* pch = NULL;
  for (size_t i = 0; (pch = model.getPatch(i)); i++)
    for (size_t inod = 1; inod <= pch->getNoNodes(); inod++)
      if (xp.equal(pch->getCoord(inod)))
      {
	int node = pch->getNodeID(inod);
	ptMap[xp] = node;
	return node;
      }

  std::cerr <<" *** NodeVecFunc::getPointIndex: No nodes matches the point "
            << xp << std::endl;
  return -1;
}
