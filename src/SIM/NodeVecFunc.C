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
#include "Vec3Oper.h"


bool NodeVecFunc::isZero () const
{
  if (value && value->empty()) return true;
  if (&idMap != &dummy) return idMap.empty();
  return false;
}


Vec3 NodeVecFunc::evaluate (const Vec3& X) const
{
  int idx = this->getPointIndex(X).second - 1;
  if (idx < 0 || !value) return Vec3();

  size_t nf = model.getNoFields();
  if (nf*idx+nf > value->size())
  {
    std::cerr <<" *** NodeVecFunc::evaluate: Index "<< nf*idx+nf
              <<" is out of range [1,"<< value->size() <<"]."<< std::endl;
    return Vec3();
  }

  return Vec3(&value->front()+nf*idx,nf);
}


std::pair<int,int> NodeVecFunc::getPointIndex (const Vec3& Xp) const
{
  // Check if the nodal index is stored in the Vec3 object itself
  const Vec4* x4 = dynamic_cast<const Vec4*>(&Xp);
  if (x4 && x4->idx > 0)
  {
    if (idMap.empty())
      return std::make_pair(x4->idx,x4->idx); // Assume 1:1 mapping
    else
    {
      std::map<int,int>::const_iterator it = idMap.find(x4->idx);
      if (it != idMap.end()) return *it;

      std::cerr <<" *** NodeVecFunc::getPointIndex: Point "<< Xp
                <<" is not present in the index map."<< std::endl;
    }
  }

  // Search among the earlier nodes found
  std::map<Vec3,int>::const_iterator it = ptMap.find(Xp);
  if (it != ptMap.end()) return std::make_pair(x4 ? x4->idx : 0, it->second);

  // Not found, search among all nodes in the model
  const ASMbase* pch = nullptr;
  for (size_t i = 1; (pch = model.getPatch(i)); i++)
    for (size_t inod = 1; inod <= pch->getNoNodes(); inod++)
      if (Xp.equal(pch->getCoord(inod)))
        return std::make_pair(x4 ? x4->idx : 0,
                              ptMap[Xp] = pch->getNodeID(inod));

  std::cerr <<" *** NodeVecFunc::getPointIndex: No nodes matches the point "
            << Xp << std::endl;
  return std::make_pair(0,0);
}
