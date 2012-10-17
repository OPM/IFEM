// $Id$
//==============================================================================
//!
//! \file SAMpatch.h
//!
//! \date Dec 10 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of FE matrices into system matrices for multi-patch models.
//!
//==============================================================================

#ifndef _SAM_PATCH_H
#define _SAM_PATCH_H

#include "SAM.h"
#include "ASMbase.h"

/*!
  \brief This is a sub-class of SAM with added functionality for spline models.
  \details It contains some additional functions for initializing the SAM-data
  for an FE model based on ASMbase patches.
*/

class SAMpatch : public SAM
{
public:
  //! \brief Empty default constructor.
  SAMpatch() {}
  //! \brief Empty destructor.
  virtual ~SAMpatch() {}

  //! \brief Allocates the dynamic arrays and populates them with data.
  //! \param[in] model All spline patches in the model
  //! \param[in] numNod Total number of unique nodes in the model
  bool init(const std::vector<ASMbase*>& model, int numNod = 0);

  //! \brief Updates the multi-point constraint array \a TTCC.
  //! \param[in] model All spline patches in the model
  //! \param[in] prevSol Previous primary solution vector in DOF-order
  virtual bool updateConstraintEqs(const std::vector<ASMbase*>& model,
				   const Vector* prevSol = 0);

  //! \brief Split the local dofs into a number of unique subdomains
  //! \param[out] locSubds Global node number for each of the subdomains
  //! \param[in]  nx, ny, nz Number of paritionings in each direction
  virtual bool getLocalSubdomains(std::vector<IntVec>& locSubds,
				  int nx = 1, int ny = 1, int nz = 1) const;

  //! \brief Split the local dofs into a number of unique subdomains
  //! \param[out] subds Global node number for each of the subdomains
  //! \param[in]  overlap Overlap between subdomains
  //! \param[in]  nx, ny, nz Number of paritionings in each direction
  virtual bool getSubdomains(std::vector<IntVec>& locSubds, int overlap = 1,
			     int nx = 1, int ny = 1, int nz = 1) const;
  
protected:
  ASMVec patch; //!< Patches

  //! \brief Initializes the nodal arrays \a MINEX, \a MADOF and \a MSC.
  bool initNodeDofs(const std::vector<ASMbase*>& model);
  //! \brief Initializes the element topology arrays \a MPMNPC and \a MMNPC.
  bool initElementConn(const std::vector<ASMbase*>& model);
  //! \brief Initializes the multi-point constraint arrays
  //! \a MPMCEQ, \a MMCEQ and \a TTCC.
  virtual bool initConstraintEqs(const std::vector<ASMbase*>& model);

  //! \brief Split the local dofs into a number of unique subdomains (1D problem)
  //! \param[in] nx Split each patch in npart smaller subdomains
  //! \param[in] minNodeId Minimum node number for each patch
  //! \param[in] maxNodeId Maximum node number for each patch
  //! \param[out] locSubds Global node number for each of the subdomains
  bool getLocalSubdomains1D(IntVec& nx, IntVec& minNodeId, IntVec& maxNodeId, 
			    std::vector<IntVec>& locSubds) const;

  //! \brief Split the local dofs into a number of unique subdomains (1D problem)
  //! \param[in] nx, ny Split each patch in npart smaller subdomains
  //! \param[in] minNodeId Minimum node number for each patch
  //! \param[in] maxNodeId Maximum node number for each patch
  //! \param[out] locSubds Global node number for each of the subdomains
  bool getLocalSubdomains2D(IntVec& nx, IntVec& ny, IntVec& minNodeId,
			    IntVec& maxNodeId, std::vector<IntVec>& locSubds) const;

  //! \brief Split the local dofs into a number of unique subdomains (1D problem)
  //! \param[in] nx, ny, nz Split each patch in npart smaller subdomains
  //! \param[in] minNodeId Minimum node number for each patch
  //! \param[in] maxNodeId Maximum node number for each patch
  //! \param[out] locSubds Global node number for each of the subdomains
  bool getLocalSubdomains3D(IntVec& nx, IntVec& ny, IntVec& nz, IntVec& minNodeId,
			    IntVec& maxNodeId, std::vector<IntVec>& locSubds) const;

   //! \brief Split the dofs into a number of overlapping subdomains (1D problem)
  //! \param[in] nx Split each patch in npart smaller subdomains
  //! \param[in] overlap Overlap for subdomains
  //! \param[out] locSubds Global node number for each of the subdomains
  bool getSubdomains1D(IntVec& nx, int overlap, std::vector<IntVec>& locSubds) const;

  //! \brief Split the local dofs into a number of unique subdomains (1D problem)
  //! \param[in] nx, ny Split each patch in npart smaller subdomains
  //! \param[in] overlap Overlap for subdomains
  //! \param[out] locSubds Global node number for each of the subdomains
  bool getSubdomains2D(IntVec& nx, IntVec& ny, int overlap, 
		       std::vector<IntVec>& locSubds) const;

  //! \brief Split the local dofs into a number of unique subdomains (1D problem)
  //! \param[in] nx, ny, nz Split each patch in npart smaller subdomains
  //! \param[in] overlap Overlap for subdomains
  //! \param[out] locSubds Global node number for each of the subdomains
  bool getSubdomains3D(IntVec& nx, IntVec& ny, IntVec& nz, int overlap,
		       std::vector<IntVec>& locSubds) const;
};

#endif
