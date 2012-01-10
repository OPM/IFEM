// $Id$
//==============================================================================
//!
//! \file SAMpatchPara.C
//!
//! \date Sep 6 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Assembly of FE matrices into system matrices for distributed models.
//!
//==============================================================================

#include "SAMpatchPara.h"
#include "PETScMatrix.h"
#include "LinAlgInit.h"
#include "ASMbase.h"
#include "MPC.h"

#ifdef HAS_PETSC
#include "petscversion.h"

#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2
#define PETSCMANGLE(x) &x
#else
#define PETSCMANGLE(x) x
#endif

#endif


SAMpatchPara::SAMpatchPara (const std::map<int,int>& g2ln)
{
  l2gn.resize(g2ln.size(),0);
  std::map<int,int>::const_iterator it;
  for (it = g2ln.begin(); it != g2ln.end(); it++)
    l2gn[it->second-1] = it->first;

#ifdef PARALLEL_PETSC
  MPI_Comm_size(PETSC_COMM_WORLD,&nProc);
#else
  nProc = 1;
#endif
  LinAlgInit::increfs();
}


SAMpatchPara::~SAMpatchPara ()
{
#ifdef PARALLEL_PETSC
  ISDestroy(PETSCMANGLE(iglob));
  ISDestroy(PETSCMANGLE(iloc));
#endif
  LinAlgInit::decrefs();
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
    if (!this->getElmEqns(meen,e,nelmdof))
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
  for (i = 0; i < ndof; i++) {
    d_gdof = meqn[i]-1;
    if (d_gdof >= ifirst && d_gdof < ilast)
      d_nnz[d_gdof-ifirst] = d_dofc[i].size();
  }

  // Generate nnz for off-diagonal block
  std::vector<PetscInt> l2g(ndof);
  Vector nnz(ndof);
  for (i = 0; i < ndof; i++) {
    l2g[i] = meqn[i]-1;
    nnz[i] = o_dofc[i].size();
  }

  Vec x;
  VecCreate(PETSC_COMM_WORLD,&x);
  VecSetSizes(x,locsize,PETSC_DECIDE);
  VecSetFromOptions(x);
  VecSet(x,0.0);
  VecSetValues(x,ndof,&l2g[0],&nnz[0],ADD_VALUES);
  VecAssemblyBegin(x);
  VecAssemblyEnd(x);

  PetscScalar* vec;
  VecGetArray(x,&vec);

  o_nnz.resize(locsize);
  for (i = 0; i < locsize; i++)
    o_nnz[i] = ceil(vec[i]);

  VecRestoreArray(x,&vec);
  VecDestroy(PETSCMANGLE(x));

#else
  this->SAM::getNoDofCouplings(d_nnz);
  o_nnz = IntVec(ndof,0);
#endif
  return true;
}


bool SAMpatchPara::initForAssembly (SystemVector& sysRHS,
				    Vector* reactionForces) const
{
  sysRHS.redim(nleq);
  sysRHS.init();
  if (reactionForces)
    reactionForces->resize(nspdof,true);

  return nleq > 0 ? true : false;
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

  std::vector<PetscInt> L2g;
  L2g.resize(l2g.size());
  for (size_t i=0;i<l2g.size();++i)
    L2g[i] = l2g[i];

  // Add contributions to SV (righthand side)
  VecSetValues(pvec->getVector(),eSv.size(),&L2g[0],&eSv[0],ADD_VALUES);

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

  meen.clear();
  meen.reserve(nedof);
  for (int i = 0; i < nenod; i++, ip++)
  {
    int node = mmnpc[ip-1];
    for (int j = madof[node-1]; j < madof[node]; j++)
      meen.push_back(j);
  }  
  if ((int)meen.size() == nedof) return true;   

  std::cerr <<"SAMpatchPara::getElmEqns: Invalid element matrix dimension "
	    << nedof <<" (should have been "<< meen.size() <<")"<< std::endl;
  return false;
}


bool SAMpatchPara::expandSolution (const SystemVector& solVec,
				   Vector& dofVec, real scaleSD) const
{
  if (solVec.dim() < (size_t)nleq) return false;

#ifdef PARALLEL_PETSC
  Vec solution;
  VecScatter ctx;

  dofVec.resize(ndof,true);

  SystemVector* sv  = const_cast<SystemVector*>(&solVec);
  PETScVector* svec = dynamic_cast<PETScVector*>(sv);
  if (!svec) return false;

  VecCreateSeqWithArray(PETSC_COMM_SELF,dofVec.size(),&dofVec[0],&solution);
  VecScatterCreate(svec->getVector(),iglob,solution,iloc,&ctx);
  VecScatterBegin(ctx,svec->getVector(),solution,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,svec->getVector(),solution,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterDestroy(PETSCMANGLE(ctx));
  VecDestroy(PETSCMANGLE(solution));

  return true;
#else
  return this->expandVector(solVec.getRef(),dofVec,scaleSD);
#endif
}


real SAMpatchPara::dot (const Vector& x, const Vector& y, char dofType) const
{
  real globVal = this->SAM::dot(x,y,dofType);

#ifdef PARALLEL_PETSC
  if (nProc > 1)
  {
    real locVal = globVal;

    for (size_t i = 0; i < ghostNodes.size(); i++) {
      int inod = ghostNodes[i];
      if (nodeType.empty() || nodeType[inod-1] == dofType || dofType == 'A')
	for (int j = madof[inod-1]; j < madof[inod]; j++)
	  locVal -= x(j)*y(j);
    }

    MPI_Allreduce(&locVal,&globVal,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  }
#endif

  return globVal;
}


real SAMpatchPara::normL2 (const Vector& x, char dofType) const
{
#ifdef PARALLEL_PETSC
  if (nProc > 1 && nnodGlob > 1)
    return this->norm2(x,dofType)/sqrt((madof[1]-madof[0])*nnodGlob);
  // TODO,kmo: The above is not correct for mixed methods. We need to find the
  // global number of DOFs of type dofType and use that in the denominator.
#endif
  return this->SAM::normL2(x,dofType);
}


real SAMpatchPara::normInf (const Vector& x, size_t& comp, char dofType) const
{
  real locmax = this->SAM::normInf(x,comp,dofType);
#ifdef PARALLEL_PETSC
  if (nProc > 1)
  {
    int myRank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&myRank);

    int nndof = madof[1]-madof[0];
    for (size_t i = 0; i < nodeType.size(); i++)
      if (nodeType[i] == dofType)
      {
	nndof = madof[i+1]-madof[i];
	break;
      }

    // TODO,kmo: Don't think this is correct in case of mixed methods
    comp = meqn[(comp-1)*nndof]/nndof+1;

    RealArray locval(2*nProc,0.0), globval(2*nProc,0.0);
    locval[2*myRank]   = locmax;
    locval[2*myRank+1] = 1.0*comp;
    MPI_Allreduce(&locval[0],&globval[0],2*nProc,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);

    // TODO,kmo: Is this calculation of comp correct? I doubth it...
    for (int n = 0; n < nProc; n++)
      if (globval[2*n] > locmax) {
	locmax = globval[2*n];
	comp   = (size_t)globval[2*n+1];
      }
  }
#endif
  return locmax;
}


bool SAMpatchPara::initConstraintEqs (const std::vector<ASMbase*>& model)
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
  if (nceq < 1) return true;

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

  return true;
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
      if (madof[i] < madof[i+1])
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
    MPI_Allreduce(&max,&nnodGlob,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);

    // Generate list of ghost nodes
    int m = 1;
    for (size_t k = 0; k < l2gn.size(); k++)
      if (madof[k] < madof[k+1]) 
	if (l2gn[k] < min) ghostNodes.push_back(m++);

#endif

    // TODO: Fix this for mixed methods (varying DOFs per node)
    int nndof = madof[1]-madof[0];
    nleq = (max-min+1)*nndof;
    int l = 0;
    for (i = 0; i < nnod; i++) 
      if (madof[i] < madof[i+1]) {
	for (j = 0; j < nndof; j++)
	  meqn[l*nndof+j] = (l2gn[i]-1)*nndof + j + 1;
	l++;
      }
  }
  else {
    nleq = ndof;
    for (i = 0; i < ndof; i++)
      meqn[i] = i+1;
  }

#ifdef PARALLEL_PETSC
  // Generate 0-based local-to-global dof mapping
  std::vector<PetscInt> l2g(ndof);
  for (i = 0; i < ndof; i++)
    l2g[i] = meqn[i]-1;

  // Generate global and local index sets
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2
  ISCreateGeneral(PETSC_COMM_WORLD,ndof,&l2g[0],PETSC_COPY_VALUES,&iglob);
#else
  ISCreateGeneral(PETSC_COMM_WORLD,ndof,&l2g[0],&iglob);
#endif
  ISCreateStride(PETSC_COMM_WORLD,ndof,0,1,&iloc);
#endif

  // Number of equations equals number of dofs
  neq = ndof;

  return true;
}
