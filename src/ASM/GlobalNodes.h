// $Id$
//==============================================================================
//!
//! \file GlobalNodes.h
//!
//! \date Mar 13 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Simple global node establishment for unstructured FE models.
//!
//==============================================================================

#ifndef _GLOBAL_NODES_H_
#define _GLOBAL_NODES_H_

#include "Interface.h"
#include <LRSpline/LRSplineSurface.h>

#include <vector>

class DomainDecomposition;
class ProcessAdm;
class SIMbase;


/*!
  \brief Class establishing global node numbers for unstructed FE models.
*/

class GlobalNodes
{
public:
  typedef std::vector<int> IntVec; //!< Convenience typedef
  typedef std::vector<const LR::LRSpline*> LRSplineVec; //!< Convenience typedef
  typedef std::vector<ASM::Interface> InterfaceVec; //!< Convenience typedef

  //! \brief Extract local boundary nodes for a LR spline.
  //! \param lr The LR spline to extract boundary nodes for
  //! \param dim The dimension of the boundary to extract
  //! \param lidx The local index of the boundary to extract
  //! \param orient Orientation of nodes on boundary
  static IntVec getBoundaryNodes(const LR::LRSpline& lr,
                                 int dim, int lidx, int orient);

  //! \brief Calculate global node numbers for a FE model.
   //! \param pchs The spline patches in the model
  //! \param interfaces The topological connections for the spline patches
  static std::vector<IntVec> calcGlobalNodes(const LRSplineVec& pchs,
                                             const InterfaceVec& interfaces);

  //! \brief Calculate parallel global node numbers for a FE model.
  //! \param pchs The spline patches in the model
  //! \param MLGN Process-global node numbers
  //! \param sim The simulator holding patch owner information
  //! \param adm The parallel process administrator
  //! \param dd The domain decomposition holding ghost connections
  //! \param[out] nNodes Total number of nodes in the model
  static IntVec calcDDMapping(const LRSplineVec& pchs,
                              const std::vector<GlobalNodes::IntVec>& MLGN,
                              const SIMbase& sim, int& nNodes);
};

#endif
