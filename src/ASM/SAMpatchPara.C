// $Id: SAMpatchPara.C,v 1.5 2011-02-08 11:59:51 rho Exp $
//==============================================================================
//!
//! \file SAMpatchPara.C
//!
//! \date Sep 6 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Assembly of FE matrices into system matrices for multi-patch models.
//!
//==============================================================================

#include "SAMpatchPara.h"
#include "PETScMatrix.h"
#include "ASMbase.h"
#include "MPC.h"


SAMpatchPara::SAMpatchPara(const IntVec& l2gn_mp)
{
  l2gn = l2gn_mp;
#ifdef PARALLEL_PETSC
  MPI_Comm_size(PETSC_COMM_WORLD,&nProc);
#else
  nProc = 1;
#endif
}


SAMpatchPara::~SAMpatchPara ()
{
#ifdef PARALLEL_PETSC
  ISDestroy(iglob);
  ISDestroy(iloc);
#endif
}


bool SAMpatchPara::getNoDofCouplings (int ifirst, int ilast,
				      IntVec& d_nnz, IntVec& o_nnz) const
{
#ifdef PARALLEL_PETSC
  int i, e;
  int d_ldof, d_gdof, o_ldof, o_gdof;
  size_t j, k;

  // Find number of dof couplings for each node
  std::vector<IntSet> d_dofc(ndof), o_dofc(ndof);
  for (e = 1; e <= nel; e++)
  {
    IntVec meen;
    if (!this->getElmEqns(meen,e))
      return false;

    for (j = 0; j < meen.size(); j++)
      if (meen[j] > 0) {
	d_ldof = meen[j]-1;
	d_gdof = meqn[d_ldof]-1;
	if (d_gdof >= ifirst && d_gdof < ilast) {
	  for (k = 0; k < meen.size(); k++)
	    if (meen[k] > 0) {
	      o_ldof = meen[k]-1;
	      o_gdof = meqn[o_ldof]-1;
	      if (o_gdof >= ifirst && o_gdof < ilast)
		d_dofc[d_ldof].insert(o_ldof);
	      else
		o_dofc[d_ldof].insert(o_ldof);
	    }
	}
	else {
	  for (k = 0; k < meen.size(); k++)
	    if (meen[k] > 0) {
	      o_ldof = meen[k]-1;
	      o_gdof = meqn[o_ldof]-1; 
	      if (o_gdof >= ifirst && o_gdof < ilast)
		o_dofc[d_ldof].insert(o_ldof);
	    }
	}
      }
  }

  // Generate nnz for diagonal block
  int locsize = ilast-ifirst;
  d_nnz.resize(locsize,0);
  for (i = 0;i < ndof;i++) {
    d_gdof = meqn[i]-1;
    if (d_gdof >= ifirst && d_gdof < ilast)
      d_nnz[d_gdof-ifirst] = d_dofc[i].size();
  }

  // Generate nnz for off-diagonal block
  IntVec l2g;
  Vector nnz;
  l2g.resize(ndof);
  nnz.resize(ndof);
  for (i = 0;i < ndof;i++) {
    l2g[i] = meqn[i]-1;
    nnz[i] = o_dofc[i].size();
  }

  Vec x;
  VecCreate(PETSC_COMM_WORLD,&x);
  VecSetSizes(x,locsize,PETSC_DECIDE);
  VecSetFromOptions(x);
  VecSet(x,0.0);
  VecSetValues(x,ndof,&(l2g[0]),&(nnz[0]),ADD_VALUES);
  VecAssemblyBegin(x);
  VecAssemblyEnd(x);

  PetscScalar *vec;
  VecGetArray(x,&vec);

  o_nnz.resize(locsize);
  for (i = 0;i < locsize;i++)
    o_nnz[i] = ceil(vec[i]);

  VecRestoreArray(x,&vec);
  VecDestroy(x);

#else
  this->SAM::getNoDofCouplings(d_nnz);
  o_nnz = IntVec(ndof,0);
#endif
  return true;
}


bool SAMpatchPara::assembleSystem (SystemVector& sysRHS,
				   const Matrix& eK, int iel,
				   Vector* reactionForces) const
{
  return this->SAM::assembleSystem(sysRHS,eK,iel,reactionForces);
}


bool SAMpatchPara::assembleSystem (SystemVector& sysRHS,
				   const RealArray& eS, int iel,
				   Vector* reactionForces) const
{
#ifdef PARALLEL_PETSC
  IntVec l2g;
  if (!this->getElmEqns(l2g,iel,eS.size()))
    return false;
  
  RealArray eSv(eS);
  for (size_t i = 0; i < l2g.size(); i++) {
    if (mpmceq[--l2g[i]] != 0)
      eSv[i] = real(0);
    l2g[i] = meqn[l2g[i]]-1;
  }

  PETScVector* pvec = dynamic_cast<PETScVector*>(&sysRHS);
  if (!pvec) return false;

  // Add contributions to SV (righthand side)
  VecSetValues(pvec->getVector(),nedof,&l2g[0],&eSv[0],ADD_VALUES);

  // Add contributions to reaction forces
  if (reactionForces)
    this->assembleReactions(*reactionForces,eS,iel);

  return true;
#else
  return this->SAM::assembleSystem(sysRHS,eS,iel,reactionForces);
#endif
}


bool SAMpatchPara::getElmEqns (IntVec& meen, int iel, int nedof) const
{
  if (iel < 1 || iel > nel) return false;

  int ip = mpmnpc[iel-1];
  int nenod = mpmnpc[iel] - ip;
  if (nedof < 1) nedof = nenod*ndof/nnod;

#ifdef USE_F77SAM
  int neldof, neslv, neprd;
  meen.resize(nedof,0);
  elmeq_(madof,mmnpc+ip-1,mpmceq,meqn,nenod,&meen.front(),neldof,neslv,neprd);
#else  
  meen.clear();
  meen.reserve(nedof);
  for (int i = 0; i < nenod; i++, ip++)
  {
    int node = mmnpc[ip-1];
    for (int j = madof[node-1]; j < madof[node]; j++)
      meen.push_back(j);
  }
  int neldof = meen.size();
#endif
  if (neldof == nedof) return true;

  std::cerr <<"SAMpatchPara::getElmEqns: Invalid element matrix dimension "
	    << nedof <<" (should have been "<< neldof <<")"<< std::endl;
  return false;
}


bool SAMpatchPara::expandSolution (const SystemVector& solVec,
				   Vector& dofVec) const
{
  if (solVec.dim() < (size_t)nleq) return false;

#ifdef PARALLEL_PETSC
  Vec solution;
  VecScatter ctx;

  dofVec.resize(ndof,true);

  SystemVector* sv  = const_cast<SystemVector*>(&solVec);
  PETScVector* svec = dynamic_cast<PETScVector*>(sv);
  if (!svec) return false;

  VecCreateSeqWithArray(PETSC_COMM_SELF,dofVec.size(),&(dofVec[0]),&solution);
  VecScatterCreate(svec->getVector(),iglob,solution,iloc,&ctx);
  VecScatterBegin(ctx,svec->getVector(),solution,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,svec->getVector(),solution,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterDestroy(ctx);
  VecDestroy(solution);

  return true;
#else
  return this->expandVector(solVec.getRef(),dofVec);
#endif
}


real SAMpatchPara::dot (const Vector& x, const Vector& y) const
{
  real globVal = this->SAM::dot(x,y);

#ifdef PARALLEL_PETSC
  if (nProc > 1)
  {
    real locVal = globVal;
    for (size_t i = 0; i < ghostNodes.size(); i++)
    {
      int inod = ghostNodes[i];
      if (madof[inod] - madof[inod-1] == mpar[17])
	for (int j = madof[inod-1]; j < madof[inod]; j++)
	  locVal -= x(j)*y(j);
    }

    MPI_Allreduce(&locVal,&globVal,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  }
#endif

  return globVal;
}


real SAMpatchPara::normL2 (const Vector& x) const
{
#ifdef PARALLEL_PETSC
  if (nProc > 1 && nnodGlob > 1) 
    return this->norm2(x)/sqrt(mpar[17]*nnodGlob);
#endif
  return this->SAM::normL2(x);
}


real SAMpatchPara::normInf (const Vector& x, size_t& comp) const
{
  real locmax = this->SAM::normInf(x,comp);
#ifdef PARALLEL_PETSC
  int nProc, myRank;
  MPI_Comm_size(PETSC_COMM_WORLD,&nProc);
  MPI_Comm_rank(PETSC_COMM_WORLD,&myRank);

  if (nProc > 1) {
    Vector locval, globval;
    locval.resize(2*nProc,0.0);
    globval.resize(2*nProc,0.0);
    
    comp = meqn[(comp-1)*mpar[17]]/mpar[17]+1;
    
    locval[2*myRank]   = locmax;
    locval[2*myRank+1] = 1.0*comp;
    
    MPI_Allreduce(&locval[0],&globval[0],2*nProc,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);
    
    for (size_t n = 0;n < nProc;n++)
      if (globval[2*n] > locmax) {
	locmax = globval[2*n];
	comp   = (size_t) globval[2*n+1];
      }
  }
#endif
  return locmax;
}


void SAMpatchPara::initConstraintEqs (const std::vector<ASMbase*>& model)
{
  // TODO: Rewrite this calling the parent-class method, and then replacing
  // the mpmceq array into the ndof-sized array used here (why is it needed?)

  // Estimate the size of the constraint equation array
  size_t j;
  MPCIter cit;
  for (j = 0; j < model.size(); j++)
    for (cit = model[j]->begin_MPC(); cit != model[j]->end_MPC(); cit++)
      nmmceq += 1 + (*cit)->getNoMaster();

  // Initialize the constraint equation arrays
  mpmceq = new int[ndof];
  memset(mpmceq,0,ndof*sizeof(int));
  if (nceq < 1) return;
  mmceq  = new int[nmmceq];
  ttcc   = new real[nmmceq];
  int ip = 1;
  for (j = 0; j < model.size(); j++)
    for (cit = model[j]->begin_MPC(); cit != model[j]->end_MPC(); cit++, ip++)
    {
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

      mpmceq[idof-1] = ip;
      int ipslv = ip - 1;

      mmceq[ipslv] = idof;
      ttcc[ipslv] = (*cit)->getSlave().coeff;
      msc[idof-1] = -ip;

      // Master dofs ...
      for (size_t i = 0; i < (*cit)->getNoMaster(); i++)
      {
	idof = madof[(*cit)->getMaster(i).node-1] + (*cit)->getMaster(i).dof-1;
	if (msc[idof-1] > 0)
	{
	  int ipmst = (mpmceq[idof]++) - 1;
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
	  mpmceq[idof] = mpmceq[idof-1];
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


bool SAMpatchPara::updateConstraintEqs (const std::vector<ASMbase*>& model,
					const Vector* prevSol)
{
  if (nceq < 1) return true; // No constraints in this model

  MPCIter cit;
  for (size_t j = 0; j < model.size(); j++)
    for (cit = model[j]->begin_MPC(); cit != model[j]->end_MPC(); cit++)
    {
      // Slave dof ...
      int idof = madof[(*cit)->getSlave().node-1] + (*cit)->getSlave().dof - 1;
      int ipeq = mpmceq[idof-1] - 1;
      if (msc[idof-1] > 0 || mmceq[ipeq] != idof)
      {
        std::cerr <<" *** Corrupted SAM arrays detected in update."<< std::endl;
        return false;
      }
      else if (!prevSol)
        ttcc[ipeq] = 0.0;
      else if (idof <= prevSol->size())
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


bool SAMpatchPara::initSystemEquations ()
{
  int i, j;

  // Initialize matrix-of-equation-numbers
  meqn = new int[ndof];
  if (!l2gn.empty())
  {
    int min = l2gn.front();
    int max = min;

    for (i = 1; i < nnod; i++)
      if (l2gn[i] < min)
	min = l2gn[i];
      else if (max < l2gn[i])
	max = l2gn[i];

#ifdef PARALLEL_PETSC
    int myRank, nProc;
    MPI_Status status;
    MPI_Comm_rank(PETSC_COMM_WORLD,&myRank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nProc);
    if (myRank < nProc-1)
      MPI_Send(&max,1,MPI_INT,myRank+1,101,PETSC_COMM_WORLD);
    if (myRank > 0) {
      MPI_Recv(&min,1,MPI_INT,myRank-1,101,PETSC_COMM_WORLD,&status);
      min++;
    }

    // Find number of global nodes
    // RUNAR
    //MPI_Scatter(&max,1,MPI_INT,&nnodGlob,1,MPI_INT,nProc-1,PETSC_COMM_WORLD);
    MPI_Allreduce(&max,&nnodGlob,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);

    // Generate list of ghost nodes
    for (size_t k = 0; k < l2gn.size(); k++)
      if (l2gn[k] < min) ghostNodes.push_back(k+1);
#endif

    // TODO: Fix this for mixed field interpolations (varying DOFs per node)
    int nndof = ndof/nnod;
    nleq = (max-min+1)*nndof;
    for (i = 0; i < nnod; i++)
      for (j = 0; j < nndof; j++)
	meqn[i*nndof+j] = (l2gn[i]-1)*nndof + j + 1;
  }
  else {
    nleq = ndof;
    for (i = 0; i < ndof; i++)
      meqn[i] = i+1;
  }

#ifdef PARALLEL_PETSC
  // Generate 0-based local-to-global dof mapping
  IntVec l2g(ndof);
  for (i = 0; i < ndof; i++)
    l2g[i] = meqn[i]-1;

  // Generate global and local index sets
  ISCreateGeneral(PETSC_COMM_WORLD,ndof,&(l2g[0]),&iglob);
  ISCreateStride(PETSC_COMM_WORLD,ndof,0,1,&iloc);
#endif

  // Number of equations equals number of dofs
  neq = ndof;

  // Initialize the number of nodal DOFs
  int inod, nndof;
  mpar[16] = mpar[17] = madof[1]-madof[0];
  for (int inod = 1; inod < nnod && madof; inod++)
    if ((nndof = madof[1]-madof[0]) < mpar[16])
      mpar[16] = nndof;
    else if (nndof > mpar[17])
      mpar[17] = nndof;

  return true;
}
