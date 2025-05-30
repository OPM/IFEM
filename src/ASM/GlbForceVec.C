// $Id$
//==============================================================================
//!
//! \file GlbForceVec.C
//!
//! \date Dec 13 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Storage of a global nodal force vector for a FEM problem.
//!
//==============================================================================

#include "GlbForceVec.h"
#include "ElmMats.h"
#include "Vec3.h"
#include "SAM.h"
#if SP_DEBUG > 1
#include "Vec3Oper.h"
#endif


bool GlbForceVec::initNodeMap (const std::vector<int>& globalNodes, size_t nfc)
{
  nodeMap.clear();
  nodeNum = globalNodes;

  int illegal = 0, nnod = sam.getNoNodes();
  for (size_t i = 0; i < globalNodes.size(); i++)
    if (globalNodes[i] > 0 && globalNodes[i] <= nnod)
      nodeMap[globalNodes[i]] = i+1;
    else
      illegal++;

  F.resize(nfc,nodeMap.size());
  if (illegal == 0) return true;

  std::cerr <<" *** GlbForceVec::initNodeMap: "<< illegal
            <<" node numbers (out of "<< globalNodes.size()
            <<") are out of range [1,"<< nnod <<"]."<< std::endl;
  return false;
}


void GlbForceVec::initialize (char)
{
  F.fill(0.0);
}


bool GlbForceVec::assemble (const LocalIntegral* elmObj, int elmId)
{
  const ElmMats* elm = dynamic_cast<const ElmMats*>(elmObj);
  if (!elm || elm->b.size() < 1) return false;

  const Vector& ES = elm->getRHSVector();
  const size_t nfc = F.rows();
  std::vector<int> mnpc;
  if (!sam.getElmNodes(mnpc,elmId))
    return false;
/*  else if (ES.size() < nfc*mnpc.size())
  {
    std::cerr <<" *** GlbForceVec::assemble: Invalid element force vector,"
              <<" size="<< ES.size() <<" should be (at least) "
              << nfc*mnpc.size() << std::endl;
    return false;
  }
*/
  // Assemble the nodal forces into the Matrix F
  size_t i, j, k, ninod = 0;
  std::map<int,size_t>::const_iterator nit;
  for (i = k = 0; i < mnpc.size(); i++) {
    if ((nit = nodeMap.find(mnpc[i])) == nodeMap.end())
      ninod++;
    else for (j = 0; j < nfc; j++)
      F(j+1,nit->second) -= ES[k+j];
    auto dofs = sam.getNodeDOFs(mnpc[i]);
    k += dofs.second-dofs.first+1;
  }

  if (ninod < mnpc.size())
    return true;

  std::cerr <<" *** GlbForceVec::assemble: Element "<< elmId
            <<" has no nodal force contributions on this boundary"<< std::endl;
  return false;
}


bool GlbForceVec::finalize (bool)
{
  // TODO: Add MPI reduction of the global F here
  return true;
}


Vec3 GlbForceVec::getForce (int node) const
{
  Vec3 force;
  std::map<int,size_t>::const_iterator nit = nodeMap.find(node);
  if (nit != nodeMap.end())
    force = Vec3(F.getColumn(nit->second));

#if SP_DEBUG > 1
  std::cout <<"Force in node "<< node <<": "<< force << std::endl;
#endif
  return force;
}


Vec3 GlbForceVec::getTotalForce () const
{
  Vec3 force;
  for (size_t j = 1; j <= F.cols(); j++)
    force += Vec3(F.getColumn(j));

  return force;
}


int GlbForceVec::getForce (size_t indx, Vec3& force) const
{
  if (indx < F.cols() && indx < nodeNum.size())
    force = Vec3(F.getColumn(indx+1));
  else
    return 0;

#if SP_DEBUG > 1
  std::cout <<"Force in node "<< nodeNum[indx] <<": "<< force << std::endl;
#endif
  return nodeNum[indx];
}
