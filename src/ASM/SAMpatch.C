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
#include "IFEM.h"


bool SAMpatch::init (const ASMVec& model, int numNod)
{
  patches = model;

  // Initialize some model size parameters
  nnod = numNod;
  for (size_t i = 0; i < model.size(); i++)
  {
    nel  += model[i]->getNoElms(false,true);
    nceq += model[i]->getNoMPCs();
    if (numNod == 0) nnod += model[i]->getNoNodes();
  }

  // Initialize the node/dof arrays (madof,msc) and compute ndof
  if (!this->initNodeDofs(model))
    return false;

  IFEM::cout <<"\n\n >>> SAM model summary <<<"
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
        IFEM::cout <<"Number of "<< it->first <<"-dofs      "
                   << it->second << std::endl;
  }

  // Initialize the element connectivity arrays (mpmnpc,mmnpc)
  if (!this->initElementConn(model))
    return false;

  // Initialize the constraint equation arrays (mpmceq,mmceq,ttcc)
  if (!this->initConstraintEqs(model))
    return false;
  else if (nceq > 0)
    IFEM::cout <<"Number of constraints "<< nceq << std::endl;

  // Initialize the dof-to-equation connectivity array (meqn)
  bool status = this->initSystemEquations();
  IFEM::cout <<"Number of unknowns    "<< neq << std::endl;
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

  int ierr = 0;
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
      if (nndof > 3) msc[idof1+2] *= bit->RX;
      if (nndof > 4) msc[idof1+3] *= bit->RY;
      if (nndof > 5) msc[idof1+4] *= bit->RZ;
    }

  if (ierr == 0) return true;

  std::cerr <<" *** SAMpatch::initNodeDOFs: Detected "<< ierr <<" nodes with"
	    <<" conflicting number of DOFs in adjacent patches."<< std::endl;
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
  IFEM::cout <<"Detected "<< outOfOrder <<" elements out of order, reordering..."
             << std::endl;

  std::map<int, std::pair<int,int> > sortedElms;
  for (i = 0; i < elmId.size(); i++)
    if (sortedElms.find(elmId[i]) == sortedElms.end())
      sortedElms[elmId[i]] = std::make_pair(mpmnpc[i]-1,mpmnpc[i+1]-mpmnpc[i]);
    else
    {
      std::cerr <<" *** SAMpatch::initElementConn: Multiple elements with"
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
      nmmceq += 1 + (*cit)->getNoMaster(true);

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
	std::cerr <<"SAM: Ignoring constraint equation for dof "<< idof
		  <<" ("<< (*cit)->getSlave() <<").\n"
		  <<"     This dof is already marked as SLAVE."<< std::endl;
	ip--;
	nceq--;
	continue;
      }

      (*cit)->iceq = ip-1; // index into mpmceq for this MPC equation

      int ipslv = (mpmceq[ip]++) - 1;
      mmceq[ipslv] = idof;
      ttcc[ipslv] = (*cit)->getSlave().coeff;
      msc[idof-1] = -ip;

      // Master dofs ...
      for (size_t i = 0; i < (*cit)->getNoMaster(); i++)
	if (!this->initConstraintEqMaster((*cit)->getMaster(i),
					  (*cit)->getSlave(),
					  ttcc[ipslv],ip))
	{
	  // Something is wrong, ignore this constraint equation
	  (*cit)->iceq = -1;
	  mpmceq[ip] = mpmceq[ip-1];
	  ip--;
	  nceq--;
	  break;
	}
    }

  // Reset the negative values in msc before calling SYSEQ
  for (ip = 0; ip < ndof; ip++)
    if (msc[ip] < 0) msc[ip] = 0;

  return true;
}


bool SAMpatch::initConstraintEqMaster (const MPC::DOF& master,
				       const MPC::DOF& slave,
				       Real& offset, int ip, Real scale)
{
  int idof = madof[master.node-1] + master.dof - 1;
  if (msc[idof-1] > 0)
  {
    int ipmst = (mpmceq[ip]++) - 1;
    mmceq[ipmst] = idof;
    ttcc[ipmst] = master.coeff*scale;
  }
  else if (msc[idof-1] < 0)
    // This master dof is constrained (unresolved chaining)
    if (!master.nextc)
    {
      std::cerr <<" SAM: Chained MPCs detected, slave "<< slave
		<<", master "<< master <<" (ignored)."<< std::endl;
      return false;
    }
    else
    {
      scale *= master.coeff;
      offset += master.nextc->getSlave().coeff*scale;
      for (size_t i = 0; i < master.nextc->getNoMaster(); i++)
	if (!this->initConstraintEqMaster(master.nextc->getMaster(i),
					  master.nextc->getSlave(),
					  offset,ip,scale))
	  return false;
    }

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
	std::cerr <<" *** SAM: Failed to update constraint equations from "
		  << **cit << std::endl;
	return false;
      }

      int jpeq = ipeq;
      Real c0 = (*cit)->getSlave().coeff;

      // Master dofs ...
      for (size_t i = 0; prevSol && i < (*cit)->getNoMaster(); i++)
	this->updateConstraintEqMaster((*cit)->getMaster(i),c0,ipeq);

      // Update the constant term of the constraint equation
      if (!prevSol)
	ttcc[jpeq] = 0.0;
      else if (idof <= (int)prevSol->size())
	ttcc[jpeq] = c0 - (*prevSol)(idof);
      else
	ttcc[jpeq] = c0;

#if SP_DEBUG > 1
      std::cout <<"Updated constraint value for dof="<< idof <<": "<< ttcc[jpeq]
                <<" from "<< **cit;
#endif
    }

  return true;
}


void SAMpatch::updateConstraintEqMaster (const MPC::DOF& master,
					 Real& offset, int& ipeq, Real scale)
{
  int idof = madof[master.node-1] + master.dof - 1;
  if (msc[idof-1] > 0 && mmceq[++ipeq] == idof)
    ttcc[ipeq] = master.coeff*scale;
  else if (master.nextc)
  {
    scale *= master.coeff;
    offset += master.nextc->getSlave().coeff*scale;
    for (size_t i = 0; i < master.nextc->getNoMaster(); i++)
      this->updateConstraintEqMaster(master.nextc->getMaster(i),
				     offset,ipeq,scale);
  }
}
