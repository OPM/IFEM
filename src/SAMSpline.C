// $Id: SAMSpline.C,v 1.3 2009-08-18 08:49:56 kmo Exp $
//==============================================================================
//!
//! \file SAMSpline.C
//!
//! \date Dec 10 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of FE matrices into system matrices for SplineVolume models.
//!
//==============================================================================

#include "SAMSpline.h"
#include "VolumePatch.h"
#include "LinEqSystem.h"
#include "MPC.h"


bool SAMSpline::init (const VolumePatch& patch)
{
  // Initialize some model size parameters
  nnod = patch.getNoNodes();
  nel  = patch.getNoElms();
  nceq = patch.getNoMPCs();
  ndof = 3*nnod;

  std::vector<VolumePatch*> smod(1,(VolumePatch*)&patch);

  // Initialize the node/dof arrays (madof,msc)
  initNodeDofs(smod);

  // Initialize the element connectivity arrays (mpmnpc,mmnpc)
  initElementConn(smod);

  // Initialize the constraint equation arrays (mpmceq,mmceq,ttcc)
  initConstraintEqs(smod);

  // Initialize the dof-to-equation connectivity array (meqn)
  return initSystemEquations();
}


bool SAMSpline::init (const std::vector<VolumePatch*>& smod, int numNod)
{
  // Initialize some model size parameters
  nnod = numNod;
  for (size_t i = 0; i < smod.size(); i++)
  {
    if (numNod == 0) nnod += smod[i]->getNoNodes();
    nel  += smod[i]->getNoElms();
    nceq += smod[i]->getNoMPCs();
  }
  ndof = 3*nnod;
  std::cout <<"\n\n >>> SAM model summary <<<"
	    <<"\nNumber of elements "<< nel
	    <<"\nNumber of nodes    "<< nnod
	    <<"\nNumber of dofs     "<< ndof << std::endl;

  // Initialize the node/dof arrays (madof,msc)
  initNodeDofs(smod);

  // Initialize the element connectivity arrays (mpmnpc,mmnpc)
  initElementConn(smod);

  // Initialize the constraint equation arrays (mpmceq,mmceq,ttcc)
  initConstraintEqs(smod);

  // Initialize the dof-to-equation connectivity array (meqn)
  bool status = initSystemEquations();
  std::cout <<"Number of unknowns "<< neq << std::endl;
  return status;
}


void SAMSpline::initNodeDofs (const std::vector<VolumePatch*>& smod)
{
  if (nnod < 1) return;

  // Initialize the node and dof arrays
  madof = new int[nnod+1];
  msc   = new int[ndof];

  int n; madof[0] = 1;
  for (n = 0; n < nnod; n++)
    madof[n+1] = madof[n] + 3;

  for (n = 0; n < ndof; n++)
    msc[n] = 1;

  std::vector<VolumePatch::BC>::const_iterator bit;
  for (size_t j = 0; j < smod.size(); j++)
    for (bit = smod[j]->begin_BC(); bit != smod[j]->end_BC(); bit++)
    {
      n = bit->node;
      msc[3*n-3] *= bit->CX;
      msc[3*n-2] *= bit->CY;
      msc[3*n-1] *= bit->CZ;
    }
}


void SAMSpline::initElementConn (const std::vector<VolumePatch*>& smod)
{
  if (nel < 1) return;

  // Find the size of the element connectivity array
  size_t j;
  IntMat::const_iterator eit;
  for (j = 0; j < smod.size(); j++)
    for (eit = smod[j]->begin_elm(); eit != smod[j]->end_elm(); eit++)
      nmmnpc += eit->size();

  // Initialize the element connectivity arrays
  mpmnpc = new int[nel+1];
  mmnpc  = new int[nmmnpc];
  int ip = mpmnpc[0] = 1;
  for (j = 0; j < smod.size(); j++)
    for (eit = smod[j]->begin_elm(); eit != smod[j]->end_elm(); eit++, ip++)
    {
      mpmnpc[ip] = mpmnpc[ip-1];
      for (size_t i = 0; i < eit->size(); i++)
	mmnpc[(mpmnpc[ip]++)-1] = smod[j]->getNodeID(1+(*eit)[i]);
    }
}


void SAMSpline::initConstraintEqs (const std::vector<VolumePatch*>& smod)
{
  // Estimate the size of the constraint equation array
  size_t j;
  MPCIter cit;
  for (j = 0; j < smod.size(); j++)
    for (cit = smod[j]->begin_MPC(); cit != smod[j]->end_MPC(); cit++)
      nmmceq += 1 + (*cit)->getNoMaster();

  // Initialize the constraint equation arrays
  mpmceq = new int[nceq+1];
  int ip = mpmceq[0] = 1;
  if (nceq < 1) return;
  mmceq  = new int[nmmceq];
  ttcc   = new real[nmmceq];
  for (j = 0; j < smod.size(); j++)
    for (cit = smod[j]->begin_MPC(); cit != smod[j]->end_MPC(); cit++, ip++)
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
    }

  // Reset the negative values in msc before calling SYSEQ
  for (ip = 0; ip < ndof; ip++)
    if (msc[ip] < 0) msc[ip] = 0;
}


bool SAMSpline::initForAssembly (LinEqSystem& sys, bool withRHS) const
{
  if (sys.K) sys.K->initAssembly(*this);
  if (sys.M) sys.M->initAssembly(*this);
  if (withRHS) sys.RHS.resize(neq,true);
  return neq > 0 ? true : false;
}
