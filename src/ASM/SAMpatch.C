// $Id$
//==============================================================================
//!
//! \file SAMpatch.C
//!
//! \date Dec 10 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of FE matrices into system matrices for multi-patch models.
//!
//==============================================================================

#include "SAMpatch.h"
#include "ASMbase.h"
#include "ASMs1D.h"
#include "ASMs2D.h"
#include "ASMs3D.h"
#include "MPC.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"

bool SAMpatch::init (const ASMVec& model, int numNod)
{
  // Get local 2 global node mapping for each patch
  patch = model;

  // Initialize some model size parameters
  nnod = numNod;
  for (size_t i = 0; i < model.size(); i++)
  {
    nel  += model[i]->getNoElms();
    nceq += model[i]->getNoMPCs();
    if (numNod == 0) nnod += model[i]->getNoNodes();
  }

  // Initialize the node/dof arrays (madof,msc) and compute ndof
  if (!this->initNodeDofs(model))
    return false;

  std::cout <<"\n\n >>> SAM model summary <<<"
	    <<"\nNumber of elements    "<< nel
	    <<"\nNumber of nodes       "<< nnod
	    <<"\nNumber of dofs        "<< ndof << std::endl;

  if (!nodeType.empty())
  {
    // Count the number of DOFs in each basis
    std::map<char,size_t>::const_iterator it;
    std::map<char,size_t> ndofs;
    ndofs['D'] = ndofs['L'] = ndofs['P'] = ndofs['X'] = 0;
    for (size_t n = 0; n < nodeType.size(); n++)
      ndofs[nodeType[n]] += madof[n+1] - madof[n];
    for (it = ndofs.begin(); it != ndofs.end(); it++)
      if (it->second > 0)
	std::cout <<"Number of "<< it->first <<"-dofs      "
		  << it->second << std::endl;
  }

  // Initialize the element connectivity arrays (mpmnpc,mmnpc)
  if (!this->initElementConn(model))
    return false;

  // Initialize the constraint equation arrays (mpmceq,mmceq,ttcc)
  if (!this->initConstraintEqs(model))
    return false;
  else if (nceq > 0)
    std::cout <<"Number of constraints "<< nceq << std::endl;

  // Initialize the dof-to-equation connectivity array (meqn)
  bool status = this->initSystemEquations();
  std::cout <<"Number of unknowns    "<< neq << std::endl;
  return status;
}


bool SAMpatch::initNodeDofs (const ASMVec& model)
{
  if (nnod < 1) return true;

  int n;
  char t;
  size_t i, j;

  // Initialize the array of accumulated DOFs for the nodes
  madof = new int[nnod+1];
  memset(madof,0,(nnod+1)*sizeof(int));

  bool ierr = 0;
  for (i = 0; i < model.size(); i++)
    for (j = 0; j < model[i]->getNoNodes(); j++)
      if ((n = model[i]->getNodeID(j+1)) > 0 && n <= nnod)
      {
	if (madof[n] == 0)
	  madof[n] = model[i]->getNodalDOFs(j+1);
	else if (madof[n] != model[i]->getNodalDOFs(j+1))
	  ierr++;

	// Define the node type for mixed methods (used by norm evaluations)
	t = model[i]->getNodeType(j+1);
	if (t != 'D')
	{
	  if (nodeType.empty()) nodeType.resize(nnod,'D');
	  nodeType[n-1] = t;
	}
      }

  madof[0] = 1;
  for (n = 0; n < nnod; n++)
    madof[n+1] += madof[n];

  for (i = 0; i < model.size(); i++)
    model[i]->initMADOF(madof);

  // Initialize the array of DOF status codes
  ndof = madof[nnod]-1;
  msc = new int[ndof];
  for (n = 0; n < ndof; n++)
    msc[n] = 1;

  ASMbase::BCVec::const_iterator bit;
  for (j = 0; j < model.size(); j++)
    for (bit = model[j]->begin_BC(); bit != model[j]->end_BC(); bit++)
    {
      int idof1 = madof[bit->node-1];
      int nndof = madof[bit->node] - idof1;
      if (nndof > 0) msc[idof1-1] *= bit->CX;
      if (nndof > 1) msc[idof1  ] *= bit->CY;
      if (nndof > 2) msc[idof1+1] *= bit->CZ;
    }

  if (ierr == 0) return true;

  std::cerr <<" *** SAMpatch::initNodeDOFs: Detected "<< ierr <<" nodes"
	    <<" nodes with conflicting number of DOFs in adjacent patches."
	    << std::endl;
  return false;
}


bool SAMpatch::initElementConn (const ASMVec& model)
{
  if (nel < 1) return true;

  // Find the size of the element connectivity array
  size_t i, j;
  IntMat::const_iterator eit;
  IntVec::const_iterator nit;
  for (j = 0; j < model.size(); j++)
    for (i = 1, eit = model[j]->begin_elm(); eit != model[j]->end_elm(); eit++)
      if (model[j]->getElmID(i++) > 0)
	nmmnpc += eit->size();

  IntVec elmId;
  elmId.reserve(nel);
  int id, outOfOrder = 0;

  // Initialize the element connectivity arrays
  mpmnpc = new int[nel+1];
  mmnpc  = new int[nmmnpc];
  int ip = mpmnpc[0] = 1;
  for (j = 0; j < model.size(); j++)
    for (i = 1, eit = model[j]->begin_elm(); eit != model[j]->end_elm(); eit++)
      if ((id = model[j]->getElmID(i++)) > 0)
      {
	mpmnpc[ip] = mpmnpc[ip-1];
	for (nit = eit->begin(); nit != eit->end(); nit++)
	  if (*nit == -2147483648) // Hack for node 0: Using -maxint as flag
	    mmnpc[(mpmnpc[ip]++)-1] = -model[j]->getNodeID(1);
	  else if (*nit < 0)
	    mmnpc[(mpmnpc[ip]++)-1] = -model[j]->getNodeID(1-(*nit));
	  else
	    mmnpc[(mpmnpc[ip]++)-1] =  model[j]->getNodeID(1+(*nit));

	// Check that the elements are in consequtive order
	if ((ip++) > 1 && id <= elmId.back())
	  outOfOrder++;

	elmId.push_back(id);
      }

  if (outOfOrder == 0) return true;

  // We need to sort the elements in increasing external element numbers
  std::cout <<"Detected "<< outOfOrder <<" elements out of order, reordering..."
	    << std::endl;

  std::map<int, std::pair<int,int> > sortedElms;
  for (i = 0; i < elmId.size(); i++)
    if (sortedElms.find(elmId[i]) == sortedElms.end())
      sortedElms[elmId[i]] = std::make_pair(mpmnpc[i]-1,mpmnpc[i+1]-mpmnpc[i]);
    else
    {
      std::cerr <<" *** SAMpatch::initElementConn: Multiple elements with "
		<<" external ID "<< elmId[i] <<" detected."<< std::endl;
      return false;
    }

  // Create new element connectivity arrays
  int* new_mpmnpc = new int[nel+1];
  int* new_mmnpc  = new int[nmmnpc];
  ip = new_mpmnpc[0] = 1;
  std::map<int, std::pair<int,int> >::const_iterator it;
  for (it = sortedElms.begin(); it != sortedElms.end(); it++, ip++)
  {
    int nen = it->second.second;
    new_mpmnpc[ip] = new_mpmnpc[ip-1] + nen;
    memcpy(new_mmnpc+new_mpmnpc[ip-1]-1,mmnpc+it->second.first,nen*sizeof(int));
  }

  // Replace the old ones...
  delete[] mpmnpc;
  delete[] mmnpc;
  mpmnpc = new_mpmnpc;
  mmnpc  = new_mmnpc;
  return true;
}


bool SAMpatch::initConstraintEqs (const ASMVec& model)
{
  // Estimate the size of the constraint equation array
  size_t j;
  MPCIter cit;
  for (j = 0; j < model.size(); j++)
    for (cit = model[j]->begin_MPC(); cit != model[j]->end_MPC(); cit++)
      nmmceq += 1 + (*cit)->getNoMaster();

  // Initialize the constraint equation arrays
  mpmceq = new int[nceq+1];
  int ip = mpmceq[0] = 1;
  if (nceq < 1) return true;

  mmceq  = new int[nmmceq];
  ttcc   = new Real[nmmceq];
  for (j = 0; j < model.size(); j++)
    for (cit = model[j]->begin_MPC(); cit != model[j]->end_MPC(); cit++, ip++)
    {
      mpmceq[ip] = mpmceq[ip-1];

      // Slave dof ...
      int idof = madof[(*cit)->getSlave().node-1] + (*cit)->getSlave().dof - 1;
      if (msc[idof-1] == 0)
      {
	std::cerr <<"SAM: Ignoring constraint equation for dof "
		  << idof <<" ("<< (*cit)->getSlave()
		  <<").\n     This dof is already marked as FIXED."<< std::endl;
	ip--;
	nceq--;
	continue;
      }
      else if (msc[idof-1] < 0)
      {
	std::cerr <<"SAM: Ignoring constraint equation for dof "
		  << idof <<" ("<< (*cit)->getSlave()
		  <<").\n     This dof is already marked as SLAVE."<< std::endl;
	ip--;
	nceq--;
	continue;
      }

      int ipslv = (mpmceq[ip]++) - 1;
      mmceq[ipslv] = idof;
      ttcc[ipslv] = (*cit)->getSlave().coeff;
      msc[idof-1] = -ip;

      // Master dofs ...
      for (size_t i = 0; i < (*cit)->getNoMaster(); i++)
      {
	idof = madof[(*cit)->getMaster(i).node-1] + (*cit)->getMaster(i).dof-1;
	if (msc[idof-1] > 0)
	{
	  int ipmst = (mpmceq[ip]++) - 1;
	  mmceq[ipmst] = idof;
	  ttcc[ipmst] = (*cit)->getMaster(i).coeff;
	}
	else if (msc[idof-1] < 0)
	{
	  // Master dof is constrained (unresolved chaining)
	  std::cerr <<"SAM: Chained MPCs detected"
		    <<", slave "<< (*cit)->getSlave()
		    <<", master "<< (*cit)->getMaster(i)
		    <<" (ignored)."<< std::endl;
	  mpmceq[ip] = mpmceq[ip-1];
	  ip--;
	  nceq--;
	  break;
	}
      }

      (*cit)->iceq = ip-1; // index into mpmceq for this MPC equation
    }

  // Reset the negative values in msc before calling SYSEQ
  for (ip = 0; ip < ndof; ip++)
    if (msc[ip] < 0) msc[ip] = 0;

  return true;
}


bool SAMpatch::updateConstraintEqs (const ASMVec& model, const Vector* prevSol)
{
  if (nceq < 1) return true; // No constraints in this model

  MPCIter cit;
  for (size_t j = 0; j < model.size(); j++)
    for (cit = model[j]->begin_MPC(); cit != model[j]->end_MPC(); cit++)
    {
      if ((*cit)->iceq < 0) continue; // Skip the ignored constraint equations

      // Slave dof ...
      int idof = madof[(*cit)->getSlave().node-1] + (*cit)->getSlave().dof - 1;
      int ipeq = mpmceq[(*cit)->iceq] - 1;
      if (msc[idof-1] > 0 || mmceq[ipeq] != idof)
      {
	std::cerr <<" *** Corrupted SAM arrays detected in update."<< std::endl;
	return false;
      }
      else if (!prevSol)
	ttcc[ipeq] = 0.0;
      else if (idof <= (int)prevSol->size())
	ttcc[ipeq] = (*cit)->getSlave().coeff - (*prevSol)(idof);
      else
	ttcc[ipeq] = (*cit)->getSlave().coeff;

      // Master dofs ...
      for (size_t i = 0; prevSol && i < (*cit)->getNoMaster(); i++)
      {
	idof = madof[(*cit)->getMaster(i).node-1] + (*cit)->getMaster(i).dof-1;
	if (msc[idof-1] > 0 && mmceq[++ipeq] == idof)
	  ttcc[ipeq] = (*cit)->getMaster(i).coeff;
      }
    }

  return true;
}


bool SAMpatch::getLocalSubdomains(std::vector<IntVec>& locSubds,
				  int nx, int ny, int nz) const
{
  // Define some parameters
  const int npatch = patch.size();
  const int nsd    = patch[0]->getNoSpaceDim();

  // Find min and max node for each patch on this processor
  IntVec maxNodeId, minNodeId;
  maxNodeId.resize(npatch,true);
  minNodeId.resize(npatch,true);

  for (int n = 0;n < npatch;n++) {
    const IntVec& MLGN = patch[n]->getGlobalNodeNums();
  
    int min = 0;
    int max = 0;
    for (size_t i = 0;i < MLGN.size();i++) {
      if (MLGN[i] > max)
	max = MLGN[i];
      if (MLGN[i] < min)
	min = MLGN[i];
    }
    minNodeId[n] = min;
    maxNodeId[n] = max;
  }
  
  minNodeId[0] = 1;
  for (int n = 1;n < npatch;n++)
    minNodeId[n] = maxNodeId[n-1] + 1;
  
  switch (nsd) {
  case 1:
  {
    IntVec nxVec; 
    nxVec.assign(npatch,nx);
    return this->getLocalSubdomains1D(nxVec,minNodeId,maxNodeId,locSubds);
  }
  case 2:
  {
    IntVec nxVec(npatch); nxVec.assign(npatch,nx);
    IntVec nyVec(npatch); nyVec.assign(npatch,ny);
    return this->getLocalSubdomains2D(nxVec,nyVec,minNodeId,maxNodeId,locSubds);
  }
  case 3:
  {
    IntVec nxVec(npatch); nxVec.assign(npatch,nx);
    IntVec nyVec(npatch); nyVec.assign(npatch,ny);
    IntVec nzVec(npatch); nzVec.assign(npatch,nz);
    return this->getLocalSubdomains3D(nxVec,nyVec,nzVec,minNodeId,maxNodeId,locSubds);
  }
  default:
    return false;
  }
}


bool SAMpatch::getSubdomains(std::vector<IntVec>& subds, int overlap,
				  int nx, int ny, int nz) const
{
  // Define some parameters
  const int npatch = patch.size();
  const int nsd    = patch[0]->getNoSpaceDim();

  switch (nsd) {
  case 1:
  {
    IntVec nxVec(npatch); nxVec.assign(npatch,nx);
    return this->getSubdomains1D(nxVec,overlap,subds);
  }
  case 2:
  {
    IntVec nxVec(npatch); nxVec.assign(npatch,nx);
    IntVec nyVec(npatch); nyVec.assign(npatch,ny);
    return this->getSubdomains2D(nxVec,nyVec,overlap,subds);
  }
  case 3:
  {
    IntVec nxVec(npatch); nxVec.assign(npatch,nx);
    IntVec nyVec(npatch); nyVec.assign(npatch,ny);
    IntVec nzVec(npatch); nzVec.assign(npatch,nz);
    return this->getSubdomains3D(nxVec,nyVec,nzVec,overlap,subds);
  }
  default:
    return false;
  }
}


bool SAMpatch::getLocalSubdomains1D(IntVec& nxvec, IntVec& minNodeId, IntVec& maxNodeId, 
				    std::vector<IntVec>& locSubds)  const
{
  // Define some parameters
  const size_t npatch = patch.size();

  if (nxvec.size() != npatch)
    return false;

  // Split the patches into smaller subdomains
  for (size_t n = 0;n < npatch;n++) {
    int nnod;
    ASMs1D* pch1 = dynamic_cast<ASMs1D*>(patch[n]);
    if (pch1)
	nnod = pch1->getCurve()->numCoefs();
      else
	return false;

    int nx = nxvec[n];
    int n1 = nnod/nx;
    
    const IntVec& MLGN = patch[n]->getGlobalNodeNums();
    for (int p1 = 0;p1 < nx;p1++) {
      IntVec subdDofs;
      
      int i1 = p1*n1;
      int i2 = (p1+1)*n1;
      if ((p1 == nx-1) && (i2 < nnod)) i2 = nnod;
      
      for (int i = i1;i < i2;i++) {
	int globNode = MLGN[i];
	
	if (globNode >= minNodeId[n] && globNode <= maxNodeId[n]) {
	  int nodedof = madof[globNode-1];
	  int nnodedof = patch[n]->getNodalDOFs(i+1);
	  for (int m = 0;m < nnodedof;m++) {
	    int ieq = meqn[nodedof+m-1];
	    if (ieq > 0)
	      subdDofs.push_back(ieq-1);
	  }
	}
      }
      
      locSubds.push_back(subdDofs);
    }
  }
  
  return true;
}


bool SAMpatch::getLocalSubdomains2D(IntVec& nxvec, IntVec& nyvec, IntVec& minNodeId, 
				    IntVec& maxNodeId, std::vector<IntVec>& locSubds) const 
{
  // Define some parameters
  const size_t npatch = patch.size();

  if (nxvec.size() != npatch || nyvec.size() != npatch)
    return false;

  // Split the patches into smaller subdomains
  for (size_t n = 0;n < npatch;n++) {
    int nnod1, nnod2;
    ASMs2D* pch2 = dynamic_cast<ASMs2D*>(patch[n]);
    if (pch2) {
      nnod1 = pch2->getSurface()->numCoefs_u();
      nnod2 = pch2->getSurface()->numCoefs_v();
    }
    else
      return false;

    int nx = nxvec[n];
    int ny = nyvec[n];
    int n1 = nnod1/nx;
    int n2 = nnod2/ny;
    
    const IntVec& MLGN = patch[n]->getGlobalNodeNums();
    for (int p2 = 0;p2 < ny;p2++) {
      int j1 = p2*n2;
      int j2 = (p2+1)*n2;
      if ((p2 == ny-1) && (j2 < nnod2)) j2 = nnod2;
      
      for (int p1 = 0;p1 < nx;p1++) {
	IntVec subdDofs;
	
	int i1 = p1*n1;
	int i2 = (p1+1)*n1;
	if ((p1 == nx-1) && (i2 < nnod1)) i2 = nnod1;

	for (int j = j1;j < j2;j++)
	  for (int i = i1;i < i2;i++) {
	    int locNode = j*nnod1 + i;
	    int globNode = MLGN[locNode];
	    
	    if (globNode >= minNodeId[n] && globNode <= maxNodeId[n]) {
	      int nodedof = madof[globNode-1];
	      int nnodedof = patch[n]->getNodalDOFs(locNode+1);
	      for (int m = 0;m < nnodedof;m++) {
		int ieq = meqn[nodedof+m-1];
		if (ieq > 0)
		  subdDofs.push_back(ieq-1);
	      }
	    }
	  }
	
	locSubds.push_back(subdDofs);
      }
    }
  }

  return true;
}


bool SAMpatch::getLocalSubdomains3D(IntVec& nxvec, IntVec& nyvec, IntVec& nzvec,
				    IntVec& minNodeId, IntVec& maxNodeId, 
				    std::vector<IntVec>& locSubds)  const  
{
  // Define some parameters
  const size_t npatch = patch.size();

  if (nxvec.size() != npatch || nyvec.size() != npatch || nzvec.size() != npatch)
    return false;

  // Split the patches into smaller subdomains
  for (size_t n = 0;n < npatch;n++) {
    int nnod1, nnod2, nnod3;
    ASMs3D* pch3 = dynamic_cast<ASMs3D*>(patch[n]);
    if (pch3) {
      nnod1 = pch3->getVolume()->numCoefs(0);
      nnod2 = pch3->getVolume()->numCoefs(1);
      nnod3 = pch3->getVolume()->numCoefs(2);
    }
    else
      return false;
	
    int nx = nxvec[n];
    int ny = nyvec[n];
    int nz = nzvec[n];
    int n1 = nnod1/nx;
    int n2 = nnod2/ny;
    int n3 = nnod3/nz;
    
    const IntVec& MLGN = patch[n]->getGlobalNodeNums();
    for (int p3 = 0;p3 < nz;p3++) {
      int k1 = p3*n3;
      int k2 = (p3+1)*n3;
      if ((p3 == nz-1) && (k2 < nnod3)) k2 = nnod3;

      for (int p2 = 0;p2 < ny;p2++) {
	int j1 = p2*n2;
	int j2 = (p2+1)*n2;
	if ((p2 == ny-1) && (j2 < nnod2)) j2 = nnod2;

	for (int p1 = 0;p1 < nx;p1++) {
	  IntVec subdDofs;

	  int i1 = p1*n1;
	  int i2 = (p1+1)*n1;
	  if ((p1 == nx-1) && (i2 < nnod1)) i2 = nnod1;

	  for (int k = k1;k < k2;k++)
	    for (int j = j1;j < j2;j++)
	      for (int i = i1;i < i2;i++) {
		int locNode = k*nnod3*nnod2 + j*nnod1 + i;
		int globNode = MLGN[locNode];
		
		if (globNode >= minNodeId[n] && globNode <= maxNodeId[n]) {
		  int nodedof = madof[globNode-1];
		  int nnodedof = patch[n]->getNodalDOFs(locNode+1);
		  for (int m = 0;m < nnodedof;m++) {
		    int ieq = meqn[nodedof+m-1];
		    if (ieq > 0)
		      subdDofs.push_back(ieq-1);
		  }
		}
	      }
	    
	    locSubds.push_back(subdDofs);
	}
      }
    }
  }
  
  return true;
}


bool SAMpatch::getSubdomains1D(IntVec& nxvec, int overlap, std::vector<IntVec>& subds)  const
{
  // Define some parameters
  const size_t npatch = patch.size();

  if (nxvec.size() != npatch)
    return false;

  // Overlap
  int olow  = overlap/2;
  int ohigh = overlap/2 + overlap%2; 
    
  // Split the patches into smaller subdomains
  for (size_t n = 0;n < npatch;n++) {
    int nnod;
    ASMs1D* pch1 = dynamic_cast<ASMs1D*>(patch[n]);
    if (pch1)
	nnod = pch1->getCurve()->numCoefs();
      else
	return false;

    int nx = nxvec[n];
    int n1 = nnod/nx;
    
    const IntVec& MLGN = patch[n]->getGlobalNodeNums();
    for (int p1 = 0;p1 < nx;p1++) {
      IntVec subdDofs;
      
      int min = 0;
      int i1 = std::max(p1*n1 - olow,min);
      int i2 = std::min((p1+1)*n1+ohigh,nnod);
      if ((p1 == nx-1) && (i2 < nnod)) i2 = nnod;
      
      for (int i = i1;i < i2;i++) {
	int globNode = MLGN[i];
	
	int nodedof = madof[globNode-1];
	int nnodedof = patch[n]->getNodalDOFs(i+1);
	for (int m = 0;m < nnodedof;m++) {
	  int ieq = meqn[nodedof+m-1];
	  if (ieq > 0)
	    subdDofs.push_back(ieq-1);
	}
      }
    
      subds.push_back(subdDofs);
    }
  }

  return true;
}


bool SAMpatch::getSubdomains2D(IntVec& nxvec, IntVec& nyvec, int overlap,
			       std::vector<IntVec>& subds) const 
{
  // Define some parameters
  const size_t npatch = patch.size();

  if (nxvec.size() != npatch || nyvec.size() != npatch)
    return false;

  // Overlap
  int olow  = overlap/2;
  int ohigh = overlap/2 + overlap%2; 

  // Split the patches into smaller subdomains
  for (size_t n = 0;n < npatch;n++) {
    int nnod1, nnod2;
    ASMs2D* pch2 = dynamic_cast<ASMs2D*>(patch[n]);
    if (pch2) {
      nnod1 = pch2->getSurface()->numCoefs_u();
      nnod2 = pch2->getSurface()->numCoefs_v();
    }
    else
      return false;
 	
    int nx = nxvec[n];
    int ny = nyvec[n];
    int n1 = nnod1/nx;
    int n2 = nnod2/ny;
    
    const IntVec& MLGN = patch[n]->getGlobalNodeNums();
    for (int p2 = 0;p2 < ny;p2++) {
      int jmin = 0;
      int j1 = std::max(p2*n2-olow,jmin);
      int j2 = std::min((p2+1)*n2+ohigh,nnod2);
      if ((p2 == ny-1) && (j2 < nnod2)) j2 = nnod2;

      for (int p1 = 0;p1 < nx;p1++) {
	IntVec subdDofs;
	
	int imin = 0;
	int i1 = std::max(p1*n1-olow,imin);
	int i2 = std::min((p1+1)*n1+ohigh,nnod1);
	if ((p1 == nx-1) && (i2 < nnod1)) i2 = nnod1;
	
	for (int j = j1;j < j2;j++)
	  for (int i = i1;i < i2;i++) {
	    int locNode = j*nnod1 + i;
	    int globNode = MLGN[locNode];
	    
	    int nodedof = madof[globNode-1];
	    int nnodedof = patch[n]->getNodalDOFs(locNode+1);
	    for (int m = 0;m < nnodedof;m++) {
	      int ieq = meqn[nodedof+m-1];
	      if (ieq > 0)
		subdDofs.push_back(ieq-1);
	    }
	  }
      
	subds.push_back(subdDofs);
      }
    }
    
  }

  return true;
}


bool SAMpatch::getSubdomains3D(IntVec& nxvec, IntVec& nyvec, IntVec& nzvec,
			       int overlap, std::vector<IntVec>& subds)  const  
{
  // Define some parameters
  const size_t npatch = patch.size();

  if (nxvec.size() != npatch || nyvec.size() != npatch || nzvec.size() != npatch)
    return false;

  // Overlap
  int olow  = overlap/2;
  int ohigh = overlap/2 + overlap%2; 

  // Split the patches into smaller subdomains
  for (size_t n = 0;n < npatch;n++) {
    int nnod1, nnod2, nnod3;
    ASMs3D* pch3 = dynamic_cast<ASMs3D*>(patch[n]);
    if (pch3) {
      nnod1 = pch3->getVolume()->numCoefs(0);
      nnod2 = pch3->getVolume()->numCoefs(1);
      nnod3 = pch3->getVolume()->numCoefs(2);
    }
    else
      return false;
	
    int nx = nxvec[n];
    int ny = nyvec[n];
    int nz = nzvec[n];
    int n1 = nnod1/nx;
    int n2 = nnod2/ny;
    int n3 = nnod3/nz;
    
    const IntVec& MLGN = patch[n]->getGlobalNodeNums();
    for (int p3 = 0;p3 < nz;p3++) {
      int kmin = 0;
      int k1 = std::max(p3*n3-olow,kmin);
      int k2 = std::min((p3+1)*n3+ohigh,nnod3);
      if ((p3 == nz-1) && (k2 < nnod3)) k2 = nnod3;

      for (int p2 = 0;p2 < ny;p2++) {
	int jmin = 0;
	int j1 = std::max(p2*n2-olow,jmin);
	int j2 = std::min((p2+1)*n2+ohigh,nnod2);
	if ((p2 == ny-1) && (j2 < nnod2)) j2 = nnod2;

	for (int p1 = 0;p1 < nx;p1++) {
	  IntVec subdDofs;

	  int imin = 0;
	  int i1 = std::max(p1*n1-olow,imin);
	  int i2 = std::min((p1+1)*n1+ohigh,nnod1);
	  if ((p1 == nx-1) && (i2 < nnod1)) i2 = nnod1;

	  for (int k = k1;k < k2;k++)
	    for (int j = j1;j < j2;j++)
	      for (int i = i1;i < i2;i++) {
		int locNode = k*nnod3*nnod2 + j*nnod1 + i;
		int globNode = MLGN[locNode];
		
		int nodedof = madof[globNode-1];
		int nnodedof = patch[n]->getNodalDOFs(locNode+1);
		for (int m = 0;m < nnodedof;m++) {
		  int ieq = meqn[nodedof+m-1];
		  if (ieq > 0)
		    subdDofs.push_back(ieq-1);
		}
	      }
	  
	  subds.push_back(subdDofs);
	}
      }
    }
  }
  
  return true;
}
