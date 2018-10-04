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


bool SAMpatch::init (const ASMVec& model, int numNod,
                     const std::vector<char>& dTypes)
{
  patches = model;

  // Initialize some model size parameters
  nnod = numNod;
  for (const ASMbase* pch : model)
  {
    nel += pch->getNoElms(false,true);
    nceq += pch->getNoMPCs();
    if (numNod == 0)
      nnod += pch->getNoNodes();
  }

  // Initialize the node/dof arrays (madof,msc) and compute ndof
  if (!this->initNodeDofs(model,dTypes))
    return false;

  IFEM::cout <<"\n\n >>> SAM model summary <<<"
             <<"\nNumber of elements    "<< nel
             <<"\nNumber of nodes       "<< nnod
             <<"\nNumber of dofs        "<< ndof << std::endl;

  // Count the number of DOFs of each type
  std::map<char,size_t> ndofs;
  ndofs['D'] = ndofs['L'] = ndofs['P'] = ndofs['Q'] = ndofs['X'] = 0;
  for (size_t n = 0; n < nodeType.size(); n++)
    ndofs[nodeType[n]] += madof[n+1] - madof[n];
  for (size_t d = 0; d < dof_type.size(); d++)
    ndofs[dof_type[d]] ++;
  for (const std::pair<char,size_t>& dof : ndofs)
    if (dof.second > 0)
      IFEM::cout <<"Number of "<< dof.first <<"-dofs      "
                 << dof.second << std::endl;

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


bool SAMpatch::initNodeDofs (const ASMVec& model,
                             const std::vector<char>& dTypes)
{
  if (nnod < 1) return true;

  // Initialize the array of accumulated DOFs for the nodes
  madof = new int[nnod+1];
  memset(madof,0,(nnod+1)*sizeof(int));

  int n, ierr = 0;
  for (const ASMbase* pch : model)
    for (size_t j = 1; j <= pch->getNoNodes(); j++)
      if ((n = pch->getNodeID(j)) > 0 && n <= nnod)
      {
        if (madof[n] == 0)
          madof[n] = pch->getNodalDOFs(j);
        else if (madof[n] != pch->getNodalDOFs(j))
          ierr++;

        // Define the node type for mixed methods (used by norm evaluations)
        char nt = pch->getNodeType(j);
        if (nt != 'D')
        {
          if (nodeType.empty()) nodeType.resize(nnod,'D');
          nodeType[n-1] = nt;
        }
      }

  if (!dTypes.empty())
  {
    dof_type.reserve(nnod*dTypes.size());
    for (n = 1; n <= nnod; n++)
      if (nodeType.empty() || nodeType[n-1] == 'D')
        dof_type.insert(dof_type.end(),dTypes.begin(),dTypes.begin()+madof[n]);
      else
        dof_type.insert(dof_type.end(),madof[n],' ');
  }

  madof[0] = 1;
  for (n = 0; n < nnod; n++)
    madof[n+1] += madof[n];

  for (ASMbase* pch : model)
    pch->initMADOF(madof);

  // Initialize the array of DOF status codes
  ndof = madof[nnod]-1;
  msc = new int[ndof];
  for (n = 0; n < ndof; n++)
    msc[n] = 1;

  ASMbase::BCVec::const_iterator bit;
  for (const ASMbase* pch : model)
    for (bit = pch->begin_BC(); bit != pch->end_BC(); ++bit)
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
  size_t i;
  IntMat::const_iterator eit;
  for (const ASMbase* pch : model)
    for (i = 1, eit = pch->begin_elm(); eit != pch->end_elm(); ++eit)
      if (pch->getElmID(i++) > 0)
	nmmnpc += eit->size();

  IntVec elmId;
  elmId.reserve(nel);
  int id, outOfOrder = 0;

  // Initialize the element connectivity arrays
  mpmnpc = new int[nel+1];
  mmnpc  = new int[nmmnpc];
  int ip = mpmnpc[0] = 1;
  for (const ASMbase* pch : model)
    for (i = 1, eit = pch->begin_elm(); eit != pch->end_elm(); ++eit)
      if ((id = pch->getElmID(i++)) > 0)
      {
	mpmnpc[ip] = mpmnpc[ip-1];
	for (int inod : *eit)
	  if (inod == -2147483648) // Hack for node 0: Using -maxint as flag
	    mmnpc[(mpmnpc[ip]++)-1] = -pch->getNodeID(1);
	  else if (inod < 0)
	    mmnpc[(mpmnpc[ip]++)-1] = -pch->getNodeID(1-inod);
	  else
	    mmnpc[(mpmnpc[ip]++)-1] =  pch->getNodeID(1+inod);

	// Check that the elements are in consequtive order
	if ((ip++) > 1 && id <= elmId.back())
	  outOfOrder++;

	elmId.push_back(id);
      }

  if (outOfOrder == 0) return true;

  // We need to sort the elements in increasing external element numbers
  IFEM::cout <<"Detected "<< outOfOrder
             <<" elements out of order, reordering..."<< std::endl;

  typedef std::pair<int,int> Ipair;
  std::map<int,Ipair> sortedElms;
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
  for (const std::pair<int,Ipair>& elm : sortedElms)
  {
    int nen = elm.second.second;
    new_mpmnpc[ip] = new_mpmnpc[ip-1] + nen;
    memcpy(new_mmnpc+new_mpmnpc[ip-1]-1,mmnpc+elm.second.first,nen*sizeof(int));
    ip++;
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
  MPCIter cit;
  for (const ASMbase* pch : model)
    for (cit = pch->begin_MPC(); cit != pch->end_MPC(); ++cit)
      nmmceq += 1 + (*cit)->getNoMaster(true);

  // Initialize the constraint equation arrays
  mpmceq = new int[nceq+1];
  int ip = mpmceq[0] = 1;
  if (nceq < 1) return true;

  mmceq  = new int[nmmceq];
  ttcc   = new Real[nmmceq];
  for (const ASMbase* pch : model)
    for (cit = pch->begin_MPC(); cit != pch->end_MPC(); ++cit, ip++)
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
  for (const ASMbase* pch : model)
    for (cit = pch->begin_MPC(); cit != pch->end_MPC(); ++cit)
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
      if (!prevSol || c0 == Real(0)) // Note: c0=0 means no constant offset
        ttcc[jpeq] = Real(0);
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
