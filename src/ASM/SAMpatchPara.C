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
#include "SystemMatrix.h"
#include "LinAlgInit.h"
#include "ASMbase.h"
#include "ASMs1D.h"
#include "ASMs2D.h"
#include "ASMs3D.h"
#include "MPC.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"

#ifdef HAS_PETSC
// RUNAR
//#include "PETScMatrix.h"
#include "PETScBlockMatrix.h"
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


bool SAMpatchPara::init (const ASMVec& model, int numNod)
{
  // Get local 2 global node mapping for each patch
  patch = model;

  return SAMpatch::init(model,numNod);
}


bool SAMpatchPara::getEqns(IntVec& eqns, int f1, int f2) const
{
  int nf    = ndof/nnod;    
  int nbf   = f2-f1+1;
  int nbdof = nnod*nbf;

  if (f1 > f2 || f1 < 0 || f2 > nf-1)
    return false;

  eqns.resize(nbdof);
  for (int n = 0;n < nnod;n++) {
    int dof  = nnod*nf + f1;
    int bdof = nnod*nbf;
    
    for (int k = 0;k < nbf;k++, dof++, bdof++)
      eqns[bdof] = meqn[dof];
  }

  return true;
}


bool SAMpatchPara::getNoDofCouplings (int ifirst, int ilast,
				      IntVec& d_nnz, IntVec& o_nnz) const
{
#ifdef PARALLEL_PETSC
  int d_ldof, d_gdof, o_ldof, o_gdof;
  size_t j, k;

  // Find number of dof couplings for each node
  std::vector<IntSet> d_dofc(ndof), o_dofc(ndof);
  for (int iel = 1; iel <= nel; iel++)
  {
    IntVec meen;
    this->getElmEqns(meen,iel);

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
  int i, locsize = ilast-ifirst;
  d_nnz.resize(locsize,0);
  for (i = 0; i < ndof; i++) {
    d_gdof = meqn[i]-1;
    if (d_gdof >= ifirst && d_gdof < ilast)
      d_nnz[d_gdof-ifirst] = d_dofc[i].size();
  }

  // Generate nnz for off-diagonal block
  std::vector<PetscInt> l2g(ndof);
  RealArray nnz(ndof);
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



bool SAMpatchPara::getNoDofCouplings(int ifirst, int ilast, IntVec ncomps,
				     std::vector<std::vector<IntVec> >& d_nnz, 
				     std::vector<std::vector<IntVec> >& o_nnz) const
{
  // Find number of dofs per row for each submatrix
  size_t nblock = ncomps.size();
  size_t nf = 0;
  for (size_t n = 0;n < nblock;n++)
    nf += ncomps[n];
  size_t nlocnode = (ilast-ifirst)/nf;

  // Must have same number of dofs per node, i.e. not mixed
  if (ndof%nf > 0)
    return false;
   
  std::vector<std::vector<std::vector<IntSet> > > d_dofc, o_dofc;

  // Resize vectors for each submatrix
  d_dofc.resize(nblock); 
  for (size_t i = 0;i < nblock;i++) {
    d_dofc[i].resize(nblock);
    for (size_t j = 0;j < nblock;j++) {
      d_dofc[i][j].resize(nnod*ncomps[i]);
    }
  }

#ifdef PARALLEL_CODE
  o_dofc.resize(nblock); 
  for (size_t i = 0;i < nblock;i++) {
    o_dofc[i].resize(nblock);
    for (size_t j = 0;j < nblock;j++) 
      o_dofc[i][j].resize(nnod*ncomps[i]);
  }
#endif

  IntVec meenI, meenJ;
  for (int iel = 1;iel <= nel;iel++) {
    size_t if1 = 0;
    for (size_t i = 0;i < nblock;i++)  {
      size_t if2 = if1 + ncomps[i];
    
      this->getElmEqns(meenI,iel,if1,if2-1,true);

      size_t jf1 = 0;
      for (size_t j = 0;j < nblock;j++) {
	size_t jf2 = jf1 + ncomps[j];
	
	this->getElmEqns(meenJ,iel,jf1,jf2-1,true);
	
	for (size_t k = 0;k < meenI.size();k++) {
	  if (meenI[k] > 0) {
	    int d_ldof = meenI[k]-1;
	    int d_gdof = meqn[d_ldof]-1;
	    
	    int d_lbdof = (d_ldof/nf)*ncomps[i] + (d_ldof%nf)-if1;
     
	    if (d_gdof >= ifirst && d_gdof < ilast) {
	      for (size_t l = 0; l < meenJ.size(); l++)
		if (meenJ[l] > 0) {
		  int o_ldof = meenJ[l]-1;
		  int o_gdof = meqn[o_ldof]-1;

		  int o_lbdof = (o_ldof/nf)*ncomps[j] + (o_ldof%nf)-jf1;		  
		  if (o_gdof >= ifirst && o_gdof < ilast)
		    d_dofc[i][j][d_lbdof].insert(o_lbdof);
		  else
		    o_dofc[i][j][d_lbdof].insert(o_lbdof);
		}
	    }
	    else 
	      for (size_t l = 0; l < meenJ.size(); l++)
		if (meenJ[l] > 0) {
		  int o_ldof = meenJ[l]-1;
		  int o_gdof = meqn[o_ldof]-1;

		  int o_lbdof = (o_ldof/nf)*ncomps[j] + (o_ldof%nf)-jf1;		  
		  if (o_gdof >= ifirst && o_gdof < ilast)
		    o_dofc[i][j][d_lbdof].insert(o_lbdof);
	    }
	  }
	}
      }
    }
  }

  // Generate nnz for diagonal blocks
  d_nnz.resize(nblock);
  size_t if1 = 0;
  for (size_t i = 0;i < nblock;i++) {
    d_nnz[i].resize(nblock);

    for (size_t j = 0;j < nblock;j++) {
      d_nnz[i][j].resize(nlocnode*ncomps[i]);
      for (int k = 0; k < nnod; k++) {
	int bdof = k*ncomps[i];
	int dof  = k*nf + if1;
	for (int l = 0;l < ncomps[i];l++, bdof++, dof++) {
	  int d_gdof = meqn[dof]-1;
	  int beq = (d_gdof-ifirst)/nf*ncomps[i] + l;
	  if (d_gdof >= ifirst && d_gdof < ilast) 
	    d_nnz[i][j][beq] = d_dofc[i][j][bdof].size();
	}
      }
    }      

    if1 += ncomps[i];
  }

#ifdef PARALLEL_CODE
  // Generate nnz for off-diagonal block);
  o_nnz.resize(nblock);
  if1 = 0;
  for (size_t i = 0;i < nblock;i++) {
    o_nnz[i].resize(nblock);
    for (size_t j = 0;j < nblock;j++) {
      o_nnz[i][j].resize(nlocnode*ncomps[i]);
    
      RealArray nnz(nnod*comps[i]);
      std::vector<PetscInt> l2g(nnod*ncomps[i]);

      for (size_t k = 0; k < nnod; k++) {
	size_t bdof = k*ncomps[i];
        size_t dof  = k*nf + f1;
        for (size_t l = 0;l < ncomps[i];l++, bdof++, dof++) {
          size_t d_gdof = meqn[dof]-1;
          if (d_gdof >= ifirst && d_gdof < ilast) { 
            nnz[bdof] = d_dofc[i][j][bdof].size();
	    l2g[bdof] = (d_gdof/nf)*ncomps[i] + l; 
	  }
        }
      }
      
      size_t locsize = nlocnode*ncomps[i];
      size_t nldof   = nnod*ncomps[i];

      Vec x;
      VecCreate(PETSC_COMM_WORLD,&x);
      VecSetSizes(x,locsize,PETSC_DECIDE);
      VecSetFromOptions(x);
      VecSet(x,0.0);
      VecSetValues(x,nldof,&l2g[0],&(nnz[i][j][0]),ADD_VALUES);
      VecAssemblyBegin(x);
      VecAssemblyEnd(x);
      
      PetscScalar* vec;
      VecGetArray(x,&vec);
      
      o_nnz[i][j].resize(locsize);
      for (k = 0; k < locsize; k++)
	o_nnz[i][j][k] = ceil(vec[k]);
      
      VecRestoreArray(x,&vec);
      VecDestroy(PETSCMANGLE(x));
    }
    
    if1 += ncomps[i];
  }
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
#ifdef HAS_PETSC
  IntVec l2g;
  if (!this->getElmEqns(l2g,iel,eS.size()))
    return false;

  size_t i;
  RealArray eSv(eS);
  for (i = 0; i < l2g.size(); i++) {
    if (mpmceq[--l2g[i]] != 0)
      eSv[i] = Real(0);
    l2g[i] = meqn[l2g[i]]-1;
  }

  // RUNAR
  PETScVector* pvec = dynamic_cast<PETScVector*>(&sysRHS);
  //PETScBlockVector* pvec = dynamic_cast<PETScBlockVector*>(&sysRHS);
  if (!pvec) return false;

  std::vector<PetscInt> L2g(l2g.size());
  for (i = 0; i < l2g.size(); i++)
    L2g[i] = l2g[i];

  // RUNAR
  // Add contributions to SV (righthand side)
  VecSetValues(pvec->getVector(),eSv.size(),&L2g[0],&eSv[0],ADD_VALUES);
  //VecSetValues(pvec->getVector(0),eSv.size(),&L2g[0],&eSv[0],ADD_VALUES);

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
  if (iel < 1 || iel > nel)
  {
    std::cerr <<"SAMpatchPara::getElmEqns: Element "<< iel
	      <<" is out of range [1,"<< nel <<"]"<< std::endl;
    return false;
  }

  int ip = mpmnpc[iel-1];
  int nenod = mpmnpc[iel] - ip;
  int oldof = nedof;
  if (nedof < 1) nedof = nenod*ndof/nnod;

  meen.clear();
  meen.reserve(nedof);
  for (int i = 0; i < nenod; i++, ip++)
  {
    int node = mmnpc[ip-1];
    if (node > 0)
      for (int j = madof[node-1]; j < madof[node]; j++)
        meen.push_back(j);
    else if (node < 0)
      meen.insert(meen.end(),madof[-node]-madof[-node-1],0);
  }
  if ((int)meen.size() == nedof || oldof < 1) return true;

  std::cerr <<"SAMpatchPara::getElmEqns: Invalid element matrix dimension "
	    << nedof <<" (should have been "<< meen.size() <<")"<< std::endl;
  return false;
}


bool SAMpatchPara::getElmEqns (IntVec& meen, int iel, int f1, int f2, bool globalEq, int nedof) const
{
  if (iel < 1 || iel > nel)
  {
    std::cerr <<"SAMpatchPara::getElmEqns: Element "<< iel
	      <<" is out of range [1,"<< nel <<"]"<< std::endl;
    return false;
  }

  int ip = mpmnpc[iel-1];
  int nenod = mpmnpc[iel] - ip;
  int oldof = nedof;
  if (nedof < 1) nedof = nenod*ndof/nnod;

  int nf      = ndof/nnod;
  int nbf     = f2-f1+1; 
  int nebdof  = nenod*nbf; 

  // Fields must be correct
  if (f2 < f1 || f1 < 0 || f2 >= nf)
    return false;

  meen.clear();
  meen.reserve(nebdof);
  for (int i = 0; i < nenod; i++, ip++)
  {
    int node = mmnpc[ip-1];
    if (node > 0) {
      int dof = (globalEq) ? madof[node-1] : (node-1)*nbf + 1 + f1;
      for (int j = f1; j <= f2; j++, dof++)
        meen.push_back(dof);
    }
    else if (node < 0) 
      meen.insert(meen.end(),f2-f1+1,0);
  }
  if ((int)meen.size() == nebdof || oldof < 1) return true;

  std::cerr <<"SAMpatchPara::getElmEqns: Invalid element matrix dimension "
	    << nedof <<" (should have been "<< meen.size() <<")"<< std::endl;
  return false;
}


bool SAMpatchPara::expandSolution (const SystemVector& solVec,
				   Vector& dofVec, Real scaleSD) const
{
  if (solVec.dim() < (size_t)nleq) return false;

#ifdef PARALLEL_PETSC
  Vec solution;
  VecScatter ctx;

  dofVec.resize(ndof,true);

  SystemVector* sv  = const_cast<SystemVector*>(&solVec);
  PETScVector* svec = dynamic_cast<PETScVector*>(sv);
  if (!svec) return false;

  VecCreateSeqWithArray(PETSC_COMM_SELF,1,dofVec.size(),&dofVec[0],&solution);
  VecScatterCreate(svec->getVector(),iglob,solution,iloc,&ctx);
  VecScatterBegin(ctx,svec->getVector(),solution,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,svec->getVector(),solution,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterDestroy(PETSCMANGLE(ctx));
  VecDestroy(PETSCMANGLE(solution));

  return true;
#else
  // RUNAR
  return this->expandVector(solVec.getRef(),dofVec,scaleSD);
  //SystemVector* sv  = const_cast<SystemVector*>(&solVec);
  //PETScBlockVector* sb = dynamic_cast<PETScBlockVector*>(sv);  
  //return this->expandVector(sb->getRef(0),dofVec,scaleSD);
#endif
}


Real SAMpatchPara::dot (const Vector& x, const Vector& y, char dofType) const
{
  Real globVal = this->SAM::dot(x,y,dofType);

#ifdef PARALLEL_PETSC
  if (nProc > 1)
  {
    Real locVal = globVal;

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


Real SAMpatchPara::normL2 (const Vector& x, char dofType) const
{
#ifdef PARALLEL_PETSC
  if (nProc > 1 && nnodGlob > 1)
    return this->norm2(x,dofType)/sqrt((madof[1]-madof[0])*nnodGlob);
  // TODO,kmo: The above is not correct for mixed methods. We need to find the
  // global number of DOFs of type dofType and use that in the denominator.
#endif
  return this->SAM::normL2(x,dofType);
}


Real SAMpatchPara::normInf (const Vector& x, size_t& comp, char dofType) const
{
  Real locmax = this->SAM::normInf(x,comp,dofType);
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
  ttcc   = new Real[nmmceq];
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


bool SAMpatchPara::getLocalSubdomains(std::vector<PetscIntVec>& locSubds,
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
  
#ifdef PARALLEL_PETSC
  // Adjust for parallel
  int myRank, nProc;
  MPI_Status status;
  MPI_Comm_rank(PETSC_COMM_WORLD,&myRank);
  MPI_Comm_size(PETSC_COMM_WORLD,&nProc);
  if (myRank < nProc-1)
    MPI_Send(&maxNodeId[npatch-1],1,MPI_INT,myRank+1,101,PETSC_COMM_WORLD);
  if (myRank > 0) {
    MPI_Recv(&minNodeId[0],1,MPI_INT,myRank-1,101,PETSC_COMM_WORLD,&status);
    minNodeId[0]++;
  }
#endif  

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


bool SAMpatchPara::getLocalSubdomainsBlock(std::vector<PetscIntVec>& locSubds, int f1, int f2,
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
  
#ifdef PARALLEL_PETSC
  // Adjust for parallel
  int myRank, nProc;
  MPI_Status status;
  MPI_Comm_rank(PETSC_COMM_WORLD,&myRank);
  MPI_Comm_size(PETSC_COMM_WORLD,&nProc);
  if (myRank < nProc-1)
    MPI_Send(&maxNodeId[npatch-1],1,MPI_INT,myRank+1,101,PETSC_COMM_WORLD);
  if (myRank > 0) {
    MPI_Recv(&minNodeId[0],1,MPI_INT,myRank-1,101,PETSC_COMM_WORLD,&status);
    minNodeId[0]++;
  }
#endif  

  switch (nsd) {
  case 1:
  {
    IntVec nxVec; 
    nxVec.assign(npatch,nx);
    return this->getLocalSubdomains1D(locSubds,nxVec,minNodeId,maxNodeId,f1,f2);
  }
  case 2:
  {
    IntVec nxVec(npatch); nxVec.assign(npatch,nx);
    IntVec nyVec(npatch); nyVec.assign(npatch,ny);
    return this->getLocalSubdomains2D(locSubds,nxVec,nyVec,minNodeId,maxNodeId,f1,f2);
  }
  case 3:
  {
    IntVec nxVec(npatch); nxVec.assign(npatch,nx);
    IntVec nyVec(npatch); nyVec.assign(npatch,ny);
    IntVec nzVec(npatch); nzVec.assign(npatch,nz);
    return this->getLocalSubdomains3D(locSubds,nxVec,nyVec,nzVec,minNodeId,maxNodeId,f1,f2);
  }
  default:
    return false;
  }
}


bool SAMpatchPara::getSubdomains(std::vector<PetscIntVec>& subds, int overlap,
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


bool SAMpatchPara::getSubdomainsBlock(std::vector<PetscIntVec>& subds, int f1, int f2, 
				      int overlap, int nx, int ny, int nz) const
{
  // Define some parameters
  const int npatch = patch.size();
  const int nsd    = patch[0]->getNoSpaceDim();

  switch (nsd) {
  case 1:
  {
    IntVec nxVec(npatch); nxVec.assign(npatch,nx);
    return this->getSubdomains1D(subds,nxVec,overlap,f1,f2);
  }
  case 2:
  {
    IntVec nxVec(npatch); nxVec.assign(npatch,nx);
    IntVec nyVec(npatch); nyVec.assign(npatch,ny);
    return this->getSubdomains2D(subds,nxVec,nyVec,overlap,f1,f2);
  }
  case 3:
  {
    IntVec nxVec(npatch); nxVec.assign(npatch,nx);
    IntVec nyVec(npatch); nyVec.assign(npatch,ny);
    IntVec nzVec(npatch); nzVec.assign(npatch,nz);
    return this->getSubdomains3D(subds,nxVec,nyVec,nzVec,overlap,f1,f2);
  }
  default:
    return false;
  }
}


bool SAMpatchPara::getLocalSubdomains1D(IntVec& nxvec, IntVec& minNodeId, IntVec& maxNodeId, 
					std::vector<PetscIntVec>& locSubds)  const
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
      PetscIntVec subdDofs;
      
      int i1 = p1*n1;
      int i2 = (p1+1)*n1;
      if ((p1 == nx-1) && (i2 < nnod)) i2 = nnod;
      
      for (int i = i1;i < i2;i++) {
	int globNode = MLGN[i];
	
	if (globNode >= minNodeId[n] && globNode <= maxNodeId[n]) {
	  int nodedof = madof[globNode-1];
	  int nnodedof = patch[n]->getNodalDOFs(i+1);
	  for (int m = 0;m < nnodedof;m++) {
	    PetscInt ieq = meqn[nodedof+m-1];
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


bool SAMpatchPara::getLocalSubdomains2D(IntVec& nxvec, IntVec& nyvec, IntVec& minNodeId, 
					IntVec& maxNodeId, std::vector<PetscIntVec>& locSubds) const 
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
	PetscIntVec subdDofs;
	
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
		PetscInt ieq = meqn[nodedof+m-1];
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


bool SAMpatchPara::getLocalSubdomains3D(IntVec& nxvec, IntVec& nyvec, IntVec& nzvec,
					IntVec& minNodeId, IntVec& maxNodeId, 
					std::vector<PetscIntVec>& locSubds)  const  
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
	  PetscIntVec subdDofs;

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
		    PetscInt ieq = meqn[nodedof+m-1];
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


bool SAMpatchPara::getSubdomains1D(IntVec& nxvec, int overlap, std::vector<PetscIntVec>& subds)  const
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
      PetscIntVec subdDofs;
      
      int min = 0;
      int i1 = std::max(p1*n1 - olow,min);
      int i2 = std::min((p1+1)*n1+ohigh,nnod);
      if ((p1 == nx-1) && (i2 < nnod)) i2 = nnod;
      
      for (int i = i1;i < i2;i++) {
	int globNode = MLGN[i];
	
	int nodedof = madof[globNode-1];
	int nnodedof = patch[n]->getNodalDOFs(i+1);
	for (int m = 0;m < nnodedof;m++) {
	  PetscInt ieq = meqn[nodedof+m-1];
	  if (ieq > 0)
	    subdDofs.push_back(ieq-1);
	}
      }
    
      subds.push_back(subdDofs);
    }
  }

  return true;
}


bool SAMpatchPara::getSubdomains2D(IntVec& nxvec, IntVec& nyvec, int overlap,
				   std::vector<PetscIntVec>& subds) const 
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
	PetscIntVec subdDofs;
	
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
	      PetscInt ieq = meqn[nodedof+m-1];
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


bool SAMpatchPara::getSubdomains3D(IntVec& nxvec, IntVec& nyvec, IntVec& nzvec,
				   int overlap, std::vector<PetscIntVec>& subds)  const  
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
	  PetscIntVec subdDofs;

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
		  PetscInt ieq = meqn[nodedof+m-1];
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


bool SAMpatchPara::getLocalSubdomains1D(std::vector<PetscIntVec>& locSubds, IntVec& nxvec, 
					IntVec& minNodeId, IntVec& maxNodeId, int f1, int f2)  const
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
    
    int nfield   = f2-f1+1;

    const IntVec& MLGN = patch[n]->getGlobalNodeNums();
    for (int p1 = 0;p1 < nx;p1++) {
      PetscIntVec subdDofs;
      
      int i1 = p1*n1;
      int i2 = (p1+1)*n1;
      if ((p1 == nx-1) && (i2 < nnod)) i2 = nnod;
      
      for (int i = i1;i < i2;i++) {
	int globNode = MLGN[i];
	
	if (globNode >= minNodeId[n] && globNode <= maxNodeId[n]) {
	  int nodedof   = madof[globNode-1];
	  int gdof      = meqn[nodedof-1];
	  int nnodedof  = patch[n]->getNodalDOFs(i+1);
	  int gbnodedof = (gdof-1)/nnodedof;
	  for (int m = 0;m < nfield;m++) {
	    PetscInt ieq = gbnodedof*nfield+m;
	    if (ieq > -1)
	      subdDofs.push_back(ieq);
	  }
	}
      }
      
      locSubds.push_back(subdDofs);
    }
  }
  
  return true;
}


bool SAMpatchPara::getLocalSubdomains2D(std::vector<PetscIntVec>& locSubds, IntVec& nxvec, IntVec& nyvec, 
					IntVec& minNodeId, IntVec& maxNodeId, int f1, int f2) const 
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
    
    int nfield = f2-f1+1;

    const IntVec& MLGN = patch[n]->getGlobalNodeNums();
    for (int p2 = 0;p2 < ny;p2++) {
      int j1 = p2*n2;
      int j2 = (p2+1)*n2;
      if ((p2 == ny-1) && (j2 < nnod2)) j2 = nnod2;
      
      for (int p1 = 0;p1 < nx;p1++) {
	PetscIntVec subdDofs;
	
	int i1 = p1*n1;
	int i2 = (p1+1)*n1;
	if ((p1 == nx-1) && (i2 < nnod1)) i2 = nnod1;

	for (int j = j1;j < j2;j++)
	  for (int i = i1;i < i2;i++) {
	    int locNode = j*nnod1 + i;
	    int globNode = MLGN[locNode];
	    
	    if (globNode >= minNodeId[n] && globNode <= maxNodeId[n]) {
	      int nodedof  = madof[globNode-1];
	      int gdof      = meqn[nodedof-1];
	      int nnodedof  = patch[n]->getNodalDOFs(locNode+1);
	      int gbnodedof = (gdof-1)/nnodedof;
	      for (int m = 0;m < nfield;m++) {
		PetscInt ieq = gbnodedof*nfield+m;
		if (ieq > -1)
		  subdDofs.push_back(ieq);
	      }
	    }
	  }
	
	locSubds.push_back(subdDofs);
      }
    }
  }

  return true;
}


bool SAMpatchPara::getLocalSubdomains3D(std::vector<PetscIntVec>& locSubds, IntVec& nxvec, IntVec& nyvec, IntVec& nzvec,
					IntVec& minNodeId, IntVec& maxNodeId, int f1, int f2)  const  
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
    
    int nfield = f2-f1+1;

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
	  PetscIntVec subdDofs;

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
		  int gdof      = meqn[nodedof-1];
		  int nnodedof = patch[n]->getNodalDOFs(locNode+1);
		  int gbnodedof = (gdof-1)/nnodedof;
		  for (int m = 0;m < nfield;m++) {
		    PetscInt ieq = gbnodedof*nfield+m;
		    if (ieq > -1)
		      subdDofs.push_back(ieq);
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


bool SAMpatchPara::getSubdomains1D(std::vector<PetscIntVec>& subds, IntVec& nxvec, int overlap, int f1, int f2)  const
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
    
    int nfield = f2-f1+1;

    const IntVec& MLGN = patch[n]->getGlobalNodeNums();
    for (int p1 = 0;p1 < nx;p1++) {
      PetscIntVec subdDofs;
      
      int min = 0;
      int i1 = std::max(p1*n1 - olow,min);
      int i2 = std::min((p1+1)*n1+ohigh,nnod);
      if ((p1 == nx-1) && (i2 < nnod)) i2 = nnod;
      
      for (int i = i1;i < i2;i++) {
	int globNode = MLGN[i];	
	int nodedof = madof[globNode-1];
	int gdof      = meqn[nodedof-1];
	int nnodedof = patch[n]->getNodalDOFs(i+1);
	int gbnodedof = (gdof-1)/nnodedof;
	for (int m = 0;m < nfield;m++) {
	  PetscInt ieq = gbnodedof*nfield+m;
	  if (ieq > -1)
	    subdDofs.push_back(ieq);
	}
      }
    
      subds.push_back(subdDofs);
    }
  }

  return true;
}


bool SAMpatchPara::getSubdomains2D(std::vector<PetscIntVec>& subds, IntVec& nxvec, IntVec& nyvec, 
				   int overlap, int f1, int f2) const 
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
    
    int nfield = f2-f1+1;

    const IntVec& MLGN = patch[n]->getGlobalNodeNums();
    for (int p2 = 0;p2 < ny;p2++) {
      int jmin = 0;
      int j1 = std::max(p2*n2-olow,jmin);
      int j2 = std::min((p2+1)*n2+ohigh,nnod2);
      if ((p2 == ny-1) && (j2 < nnod2)) j2 = nnod2;

      for (int p1 = 0;p1 < nx;p1++) {
	PetscIntVec subdDofs;
	
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

	    int gdof      = meqn[nodedof-1];
	    int gbnodedof = (gdof-1)/nnodedof;

	    for (int m = 0;m < nfield;m++) {
	      PetscInt ieq = gbnodedof*nfield+m;
	      if (ieq > -1)
		subdDofs.push_back(ieq);
	    }
	  }
      
	subds.push_back(subdDofs);
      }
    }
    
  }

  return true;
}


bool SAMpatchPara::getSubdomains3D(std::vector<PetscIntVec>& subds, IntVec& nxvec, IntVec& nyvec, IntVec& nzvec,
				   int overlap, int f1, int f2)  const  
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
    
    int nfield = f2-f2+1;

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
	  PetscIntVec subdDofs;

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

		int gdof      = meqn[nodedof-1];
		int gbnodedof = (gdof-1)/nnodedof;

		for (int m = 0;m < nfield;m++) {
		  PetscInt ieq = gbnodedof*nfield+m;
		  if (ieq > -1)
		    subdDofs.push_back(ieq);
		}
	      }
	  
	  subds.push_back(subdDofs);
	}
      }
    }
  }
  
  return true;
}
