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
  int idx = this->getPointIndex(xp) - 1;
  if (idx < 0) return Vec3();

  size_t nf = model.getNoFields();
  if (nf*idx+nf > value.size())
  {
    std::cerr <<" *** NodeVecFunc::evaluate: Index "<< nf*idx+nf
              <<" is out of range [1,"<< value.size() <<"]."<< std::endl;
    return Vec3();
  }

  return Vec3(&value.front()+nf*idx,nf);
}


int NodeVecFunc::getPointIndex (const Vec3& xp) const
{
  // Check if the nodal index is stored in the Vec3 object itself
  const Vec4* x4 = dynamic_cast<const Vec4*>(&xp);
  if (x4 && x4->idx > 0)
    if (idMap.empty())
      return x4->idx; // No index map provided, assume 1:1 mapping
    else
    {
      std::map<int,int>::const_iterator it = idMap.find(x4->idx);
      if (it != idMap.end()) return it->second;

      std::cerr <<" *** NodeVecFunc::getPointIndex: Point "<< xp
                <<" is not present in the index map."<< std::endl;
    }

  // Search among the earlier nodes found
  std::map<Vec3,int>::const_iterator it = ptMap.find(xp);
  if (it != ptMap.end()) return it->second;

  // Not found, search among all nodes in the model
  const ASMbase* pch = NULL;
  for (size_t i = 0; (pch = model.getPatch(i)); i++)
    for (size_t inod = 1; inod <= pch->getNoNodes(); inod++)
      if (xp.equal(pch->getCoord(inod)))
        return ptMap[xp] = pch->getNodeID(inod);

  std::cerr <<" *** NodeVecFunc::getPointIndex: No nodes matches the point "
            << xp << std::endl;
  return -1;
}
