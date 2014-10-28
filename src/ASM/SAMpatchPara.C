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
#include "LinAlgInit.h"
#include "ASMstruct.h"
#include "Vec3.h"
#include "PETScMatrix.h"
#include "ProcessAdm.h"


SAMpatchPara::SAMpatchPara (const std::map<int,int>& g2ln, const ProcessAdm& padm) : adm(padm)
{
  l2gn.resize(g2ln.size(),0);
  std::map<int,int>::const_iterator it;
  for (it = g2ln.begin(); it != g2ln.end(); it++)
    l2gn[it->second-1] = it->first;

  nProc = adm.getNoProcs();
  LinAlgInit::increfs();
}


SAMpatchPara::~SAMpatchPara ()
{
#ifdef HAS_PETSC
  if (adm.isParallel()) {
    ISDestroy(PETSCMANGLE(iglob));
    ISDestroy(PETSCMANGLE(iloc));
  }
#endif

  LinAlgInit::decrefs();
}


bool SAMpatchPara::init (const ASMVec& model, int numNod)
{
  // Get local 2 global node mapping for each patch
  patch = model;

  return SAMpatch::init(model,numNod);
}


bool SAMpatchPara::getNoDofCouplings (int ifirst, int ilast,
				      IntVec& d_nnz, IntVec& o_nnz) const
{
#ifdef HAS_PETSC
  if (adm.isParallel()) {
    int d_ldof, d_gdof, o_ldof, o_gdof;
    size_t j, k;
    
    // Find number of dof couplings for each node
    std::vector<IntSet> d_dofc(ndof), o_dofc(ndof);
    for (int iel = 1; iel <= nel; iel++) {
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
		// RUNAR: This gives a small overestimation for some nodes
		//o_gdof = meqn[o_ldof]-1;
		//if (o_gdof >= ifirst && o_gdof < ilast)
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
    PetscIntVec l2g(ndof);
    RealArray nnz(ndof);
    for (i = 0; i < ndof; i++) {
      l2g[i] = meqn[i]-1;
      nnz[i] = o_dofc[i].size();
    }

    Vec x;
    VecCreate(*adm.getCommunicator(),&x);
    VecSetSizes(x,locsize,PETSC_DETERMINE);
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
  }
  else {
    this->SAM::getNoDofCouplings(d_nnz);
    o_nnz = IntVec(ndof,0);
  }
#endif

  return true;
}



bool SAMpatchPara::getNoDofCouplings (int ifirst, int ilast,
				      const IntVec& ncomps,
				      std::vector<IntMat>& d_nnz,
				      std::vector<IntMat>& o_nnz) const
{
  // Find number of dofs per row for each submatrix
  size_t nblock = ncomps.size();
  size_t nf = 0;
  for (size_t n = 0;n < nblock;n++)
    nf += ncomps[n];
  int nlocnode = (ilast-ifirst)/nf;
  
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

  if (adm.isParallel()) {
    o_dofc.resize(nblock); 
    for (size_t i = 0;i < nblock;i++) {
      o_dofc[i].resize(nblock);
      for (size_t j = 0;j < nblock;j++) 
	o_dofc[i][j].resize(nnod*ncomps[i]);
    }
  }

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
		  int o_lbdof = (o_ldof/nf)*ncomps[j] + (o_ldof%nf)-jf1;
		  // RUNAR: This gives a small overestimation for some nodes
		  // int o_gdof = meqn[o_ldof]-1;
 		  // if (o_gdof >= ifirst && o_gdof < ilast)
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

#ifdef HAS_PETSC
  if (adm.isParallel()) {
    // Generate nnz for off-diagonal block);
    o_nnz.resize(nblock);
    if1 = 0;
    for (size_t i = 0;i < nblock;i++) {
      o_nnz[i].resize(nblock);
      for (size_t j = 0;j < nblock;j++) {
	o_nnz[i][j].resize(nlocnode*ncomps[i]);
	
	RealArray nnz(nnod*ncomps[i]);
	PetscIntVec l2g(nnod*ncomps[i]);
	
	for (int k = 0; k < nnod; k++) {
	  size_t bdof = k*ncomps[i];
	  size_t dof  = k*nf + if1;
	  for (int l = 0;l < ncomps[i];l++, bdof++, dof++) {
	    int d_gdof = meqn[dof]-1;
	    nnz[bdof] = o_dofc[i][j][bdof].size();
	    l2g[bdof] = (d_gdof/nf)*ncomps[i] + l; 
	  }
	}
	
	size_t locsize = nlocnode*ncomps[i];
	size_t nldof   = nnod*ncomps[i];

	Vec x;
	VecCreate(*adm.getCommunicator(),&x);
	VecSetSizes(x,locsize,PETSC_DETERMINE);
	VecSetFromOptions(x);
	VecSet(x,0.0);
	VecSetValues(x,nldof,&l2g[0],&(nnz[0]),ADD_VALUES);
	VecAssemblyBegin(x);
	VecAssemblyEnd(x);
	
	PetscScalar* vec;
	VecGetArray(x,&vec);
	
	o_nnz[i][j].resize(locsize);
	for (size_t k = 0; k < locsize; k++)
	  o_nnz[i][j][k] = ceil(vec[k]);
	
	VecRestoreArray(x,&vec);
	VecDestroy(PETSCMANGLE(x));
      }
      
      if1 += ncomps[i];
    }
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

  PETScVector* pvec = dynamic_cast<PETScVector*>(&sysRHS);
  if (!pvec) return false;

  PetscIntVec L2g(l2g.size());
  for (i = 0; i < l2g.size(); i++)
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


void SAMpatchPara::addToRHS (SystemVector& sysRHS, const RealArray& S) const
{
  Real* sysrhsPtr = sysRHS.getPtr();
  for (size_t i = 0; i < S.size(); i++) 
    if (mpmceq[i] == 0)
      sysrhsPtr[i] += S[i];
  sysRHS.restore(sysrhsPtr);
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

#ifdef HAS_PETSC
  if (adm.isParallel()) {
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
  }
#endif

  return this->expandVector(solVec.getRef(),dofVec,scaleSD);
}


Real SAMpatchPara::dot (const Vector& x, const Vector& y, char dofType) const
{
  Real globVal = this->SAM::dot(x,y,dofType);

#ifdef HAS_PETSC
  if (adm.isParallel()) {
    Real locVal = globVal;

    for (size_t i = 0; i < ghostNodes.size(); i++) {
      int inod = ghostNodes[i];
      if (nodeType.empty() || nodeType[inod-1] == dofType || dofType == 'A')
	for (int j = madof[inod-1]; j < madof[inod]; j++)
	  locVal -= x(j)*y(j);
    }

    globVal = adm.allReduce(locVal,MPI_SUM);
  }
#endif

  return globVal;
}


Real SAMpatchPara::normL2 (const Vector& x, char dofType) const
{
  if (adm.isParallel() && nnodGlob > 1)
    return this->norm2(x,dofType)/sqrt((madof[1]-madof[0])*nnodGlob);
  // TODO,kmo: The above is not correct for mixed methods. We need to find the
  // global number of DOFs of type dofType and use that in the denominator.
  return this->SAM::normL2(x,dofType);
}


Real SAMpatchPara::normInf (const Vector& x, size_t& comp, char dofType) const
{
  Real locmax = this->SAM::normInf(x,comp,dofType);

#ifdef HAS_PETSC
  if (adm.isParallel()) {
    int nndof = madof[1]-madof[0];
    for (size_t i = 0; i < nodeType.size(); i++)
      if (nodeType[i] == dofType)
      {
	nndof = madof[i+1]-madof[i];
	break;
      }

    // TODO,kmo: Don't think this is correct in case of mixed methods
    comp = meqn[(comp-1)*nndof]/nndof+1;

    int myRank = adm.getProcId();
    RealArray globval(2*nProc,0.0);
    globval[2*myRank]   = locmax;
    globval[2*myRank+1] = 1.0*comp;
    adm.allReduce(globval,MPI_MAX);
    
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


bool SAMpatchPara::getDirOrdering(PetscIntVec& order, int perm, int nf) const
{
  int nfield = (nf==0) ? madof[1]-madof[0] : nf;

  IntVec maxNodeId, minNodeId;
  if (!this->getMinMaxNode(minNodeId,maxNodeId))
    return false;

  int firstDof = (minNodeId.front()-1)*nfield;
  int gidx     = firstDof;
  int nlocdof  = (maxNodeId.back()-minNodeId.front()+1)*nfield;
  
  order.resize(nlocdof,0);
  // Split the patches into smaller subdomains
  for (size_t n = 0; n < patch.size(); n++) {
    int nnod[3];
    if (!static_cast<ASMstruct*>(patch[n])->getSize(nnod[0],nnod[1],nnod[2]))
      return false;
    if (!nnod[2]) nnod[2] = 1;

    const IntVec& MLGN = patch[n]->getGlobalNodeNums();
    switch (perm) {
    case(123):
    {
      for (int k = 0;k < nnod[2];k++)
	for (int j = 0;j < nnod[1];j++)
	  for (int i = 0;i < nnod[0];i++) {
	    int locNode  = k*nnod[0]*nnod[1] + j*nnod[0] + i;
	    int procNode = MLGN[locNode];
	    int globNode = l2gn[procNode-1];
	    
	    if (globNode >= minNodeId[n] && globNode <= maxNodeId[n])  
	      for (int m = 0;m < nfield;m++) {
		int globDof = (globNode-1)*nfield+m;
		int locDof  = globDof-firstDof;
		order[locDof] = gidx++;
	      }
	  }
      break;
    }
    case(132):
    {
      for (int j = 0;j < nnod[1];j++)
	for (int k = 0;k < nnod[2];k++)
	  for (int i = 0;i < nnod[0];i++) {
	    int locNode  = k*nnod[0]*nnod[1] + j*nnod[0] + i;
	    int procNode = MLGN[locNode];
	    int globNode = l2gn[procNode-1];
	    
	    if (globNode >= minNodeId[n] && globNode <= maxNodeId[n])  
	      for (int m = 0;m < nfield;m++) {
		int globDof = (globNode-1)*nfield+m;
		int locDof  = globDof-firstDof;
		order[locDof] = gidx++;
	      }
	  }
      break;
    }
    case(213):
    {
      for (int k = 0;k < nnod[2];k++)
	for (int i = 0;i < nnod[0];i++)
	  for (int j = 0;j < nnod[1];j++) {
	    int locNode  = k*nnod[0]*nnod[1] + j*nnod[0] + i;
	    int procNode = MLGN[locNode];
	    int globNode = l2gn[procNode-1];
	    
	    if (globNode >= minNodeId[n] && globNode <= maxNodeId[n])  
	      for (int m = 0;m < nfield;m++) {
		int globDof = (globNode-1)*nfield+m;
		int locDof  = globDof-firstDof;
		order[locDof] = gidx++;
	      }
	  }
      break;
    }
    case(231):
    {
      for (int i = 0;i < nnod[0];i++)
	for (int k = 0;k < nnod[2];k++)
	  for (int j = 0;j < nnod[1];j++) {
	    int locNode  = k*nnod[0]*nnod[1] + j*nnod[0] + i;
	    int procNode = MLGN[locNode];
	    int globNode = l2gn[procNode-1];
	    
	    if (globNode >= minNodeId[n] && globNode <= maxNodeId[n])  
	      for (int m = 0;m < nfield;m++) {
		int globDof = (globNode-1)*nfield+m;
		int locDof  = globDof-firstDof;
		order[locDof] = gidx++;
	      }
	  }
      break;
    }
    case(312):
    {
      for (int j = 0;j < nnod[1];j++)
	for (int i = 0;i < nnod[0];i++)
	  for (int k = 0;k < nnod[2];k++) {
	    int locNode  = k*nnod[0]*nnod[1] + j*nnod[0] + i;
	    int procNode = MLGN[locNode];
	    int globNode = l2gn[procNode-1];
	    
	    if (globNode >= minNodeId[n] && globNode <= maxNodeId[n])  
	      for (int m = 0;m < nfield;m++) {
		int globDof = (globNode-1)*nfield+m;
		int locDof  = globDof-firstDof;
		order[locDof] = gidx++;
	      }
	  }
      break;
    }
    case(321):
    {
      for (int i = 0;i < nnod[0];i++)
	for (int j = 0;j < nnod[1];j++)
	  for (int k = 0;k < nnod[2];k++) {
	    int locNode  = k*nnod[0]*nnod[1] + j*nnod[0] + i;
	    int procNode = MLGN[locNode];
	    int globNode = l2gn[procNode-1];
	    
	    if (globNode >= minNodeId[n] && globNode <= maxNodeId[n])  
	      for (int m = 0;m < nfield;m++) {
		int globDof = (globNode-1)*nfield+m;
		int locDof  = globDof-firstDof;
		order[locDof] = gidx++;
	      }
	  }
      break;
    }
    }
  }    

  return true;
}



bool SAMpatchPara::initConstraintEqs (const ASMVec& model)
{
  if (!this->SAMpatch::initConstraintEqs(model))
    return false;

  // Replace the mpmceq array by an equivalent ndof-sized array (why?)
  int* new_mpmceq = new int[ndof];
  memset(new_mpmceq,0,ndof*sizeof(int));

  if (nceq > 0)
  {
    int ip = 0;
    MPCIter cit;
    for (size_t j = 0; j < model.size(); j++)
      for (cit = model[j]->begin_MPC(); cit != model[j]->end_MPC(); ++cit, ip++)
      {
	// Slave dof ...
	int idof = madof[(*cit)->getSlave().node-1] + (*cit)->getSlave().dof-2;
	(*cit)->iceq = idof; // index into mpmceq for this MPC equation
	new_mpmceq[idof] = mpmceq[ip];
      }
  }

  // Swap the mpmceq array with the new one
  delete[] mpmceq;
  mpmceq = new_mpmceq;
  return true;
}


bool SAMpatchPara::initSystemEquations ()
{
  int i, j;

  const int nndof = madof[1]-madof[0];

  // Initialize matrix-of-equation-numbers
  meqn = new int[ndof];
  if (!l2gn.empty())
  {
    ieqmin = l2gn.front();
    ieqmax = ieqmin;

    for (i = 1; i < nnod; i++)
      if (madof[i] < madof[i+1])
	if (l2gn[i] < ieqmin)
	  ieqmin = l2gn[i];
	else if (ieqmax < l2gn[i])
	  ieqmax = l2gn[i];

#ifdef HAS_PETSC
    if (adm.isParallel()) {
      int myRank = adm.getProcId();
      if (myRank < nProc-1)
	adm.send(ieqmax,myRank+1);
      if (myRank > 0) {
	adm.receive(ieqmin,myRank-1);
	ieqmin++;
      }

      // Find number of global nodes
      nnodGlob = adm.allReduce(ieqmax,MPI_MAX);
      
      // Generate list of ghost nodes
      for (size_t k = 0; k < l2gn.size(); k++)
	if (madof[k] < madof[k+1])
	  if (l2gn[k] < ieqmin) ghostNodes.push_back(k+1);
    }
    else
      nnodGlob = ieqmax;
#endif
    
    // TODO: Fix this for mixed methods (varying DOFs per node)
    nleq = (ieqmax-ieqmin+1)*nndof;
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
  
  // Convert from node numbers to equation numbers
  ieqmin = (ieqmin-1)*nndof + 1;
  ieqmax *= nndof;
  
#ifdef HAS_PETSC
  if (adm.isParallel()) {
    // Generate 0-based local-to-global dof mapping
    PetscIntVec l2g(ndof);
    for (i = 0; i < ndof; i++)
      l2g[i] = meqn[i]-1;
    
    // Generate global and local index sets
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2
    ISCreateGeneral(*adm.getCommunicator(),ndof,&l2g[0],PETSC_COPY_VALUES,&iglob);
#else
    ISCreateGeneral(*adm.getCommunicator(),ndof,&l2g[0],&iglob);
#endif
    ISCreateStride(*adm.getCommunicator(),ndof,0,1,&iloc);
  }
#endif

  // Number of equations equals number of dofs
  neq = ndof;

  return true;
}


bool SAMpatchPara::getLocalSubdomains (PetscIntMat& locSubds,
				       int nx, int ny, int nz) const
{
  // Find min and max node for each patch on this processor
  IntVec maxNodeId, minNodeId;
  if (!this->getMinMaxNode(minNodeId,maxNodeId))
    return false;
    
  switch (patch.front()->getNoSpaceDim()) {
  case 1:
  {
    IntVec nxVec(patch.size(),nx);
    return this->getLocalSubdomains1D(nxVec,minNodeId,maxNodeId,locSubds);
  }
  case 2:
  {
    IntVec nxVec(patch.size(),nx);
    IntVec nyVec(patch.size(),ny);
    return this->getLocalSubdomains2D(nxVec,nyVec,minNodeId,maxNodeId,locSubds);
  }
  case 3:
  {
    IntVec nxVec(patch.size(),nx);
    IntVec nyVec(patch.size(),ny);
    IntVec nzVec(patch.size(),nz);
    return this->getLocalSubdomains3D(nxVec,nyVec,nzVec,minNodeId,maxNodeId,locSubds);
  }
  default:
    return false;
  }
}


bool SAMpatchPara::getLocalSubdomainsBlock (PetscIntMat& locSubds,
					    int f1, int f2,
					    int nx, int ny, int nz) const
{
  // Find min and max node for each patch on this processor
  IntVec maxNodeId, minNodeId;
  if (!this->getMinMaxNode(minNodeId,maxNodeId))
    return false;

  switch (patch.front()->getNoSpaceDim()) {
  case 1:
  {
    IntVec nxVec(patch.size(),nx);
    return this->getLocalSubdomains1D(locSubds,nxVec,minNodeId,maxNodeId,f1,f2);
  }
  case 2:
  {
    IntVec nxVec(patch.size(),nx);
    IntVec nyVec(patch.size(),ny);
    return this->getLocalSubdomains2D(locSubds,nxVec,nyVec,minNodeId,maxNodeId,f1,f2);
  }
  case 3:
  {
    IntVec nxVec(patch.size(),nx);
    IntVec nyVec(patch.size(),ny);
    IntVec nzVec(patch.size(),nz);
    return this->getLocalSubdomains3D(locSubds,nxVec,nyVec,nzVec,minNodeId,maxNodeId,f1,f2);
  }
  default:
    return false;
  }
}


bool SAMpatchPara::getSubdomains (PetscIntMat& subds, int overlap,
				  int nx, int ny, int nz) const
{
  switch (patch.front()->getNoSpaceDim()) {
  case 1:
  {
    IntVec nxVec(patch.size(),nx);
    return this->getSubdomains1D(nxVec,overlap,subds);
  }
  case 2:
  {
    IntVec nxVec(patch.size(),nx);
    IntVec nyVec(patch.size(),ny);
    return this->getSubdomains2D(nxVec,nyVec,overlap,subds);
  }
  case 3:
  {
    IntVec nxVec(patch.size(),nx);
    IntVec nyVec(patch.size(),ny);
    IntVec nzVec(patch.size(),nz);
    return this->getSubdomains3D(nxVec,nyVec,nzVec,overlap,subds);
  }
  default:
    return false;
  }
}


bool SAMpatchPara::getSubdomainsBlock (PetscIntMat& subds, int f1, int f2,
				       int overlap, int nx, int ny, int nz) const
{
  switch (patch.front()->getNoSpaceDim()) {
  case 1:
  {
    IntVec nxVec(patch.size(),nx);
    return this->getSubdomains1D(subds,nxVec,overlap,f1,f2);
  }
  case 2:
  {
    IntVec nxVec(patch.size(),nx);
    IntVec nyVec(patch.size(),ny);
    return this->getSubdomains2D(subds,nxVec,nyVec,overlap,f1,f2);
  }
  case 3:
  {
    IntVec nxVec(patch.size(),nx);
    IntVec nyVec(patch.size(),ny);
    IntVec nzVec(patch.size(),nz);
    return this->getSubdomains3D(subds,nxVec,nyVec,nzVec,overlap,f1,f2);
  }
  default:
    return false;
  }
}


bool SAMpatchPara::getLocalSubdomains1D (const IntVec& nxvec,
					 const IntVec& minNodeId,
					 const IntVec& maxNodeId,
					 PetscIntMat& locSubds) const
{
  if (nxvec.size() != patch.size())
    return false;

  // Split the patches into smaller subdomains
  for (size_t n = 0; n < patch.size(); n++) {

    int nnod, n2, n3;
    if (!static_cast<ASMstruct*>(patch[n])->getSize(nnod,n2,n3) || n2 || n3)
      return false;

    int nx = nxvec[n];
    int n1 = nnod/nx;
    
    const IntVec& MLGN = patch[n]->getGlobalNodeNums();
    for (int p1 = 0;p1 < nx;p1++) {
      std::set<PetscInt> subdDofs;
      
      int i1 = p1*n1;
      int i2 = (p1+1)*n1;
      if ((p1 == nx-1) && (i2 < nnod)) i2 = nnod;
      
      for (int i = i1;i < i2;i++) {
	int procNode = MLGN[i];
	int globNode = l2gn[procNode-1];

	if (globNode >= minNodeId[n] && globNode <= maxNodeId[n]) {
	  int nodedof = madof[procNode-1];
	  int nnodedof = patch[n]->getNodalDOFs(i+1);
	  for (int m = 0;m < nnodedof;m++) {
	    PetscInt ieq = meqn[nodedof+m-1];
	    if (ieq > 0)
	      subdDofs.insert(ieq-1);
	  }
	}
      }

      PetscIntVec sdofs;
      for (std::set<PetscInt>::iterator it = subdDofs.begin(); it != subdDofs.end();it++)
	sdofs.push_back(*it);
      locSubds.push_back(sdofs);
    }
  }
  
  return true;
}


bool SAMpatchPara::getLocalSubdomains2D (const IntVec& nxvec,
					 const IntVec& nyvec,
					 const IntVec& minNodeId,
					 const IntVec& maxNodeId,
					 PetscIntMat& locSubds) const
{
  if (nxvec.size() != patch.size() || nyvec.size() != patch.size())
    return false;

  // Split the patches into smaller subdomains
  for (size_t n = 0; n < patch.size(); n++) {

    int nnod1, nnod2, n3;
    if (!static_cast<ASMstruct*>(patch[n])->getSize(nnod1,nnod2,n3) || n3)
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
	std::set<PetscInt> subdDofs;
	
	int i1 = p1*n1;
	int i2 = (p1+1)*n1;
	if ((p1 == nx-1) && (i2 < nnod1)) i2 = nnod1;

	for (int j = j1;j < j2;j++)
	  for (int i = i1;i < i2;i++) {
	    int locNode = j*nnod1 + i;
	    int procNode = MLGN[locNode];
	    int globNode = l2gn[procNode-1];

	    if (globNode >= minNodeId[n] && globNode <= maxNodeId[n]) {
	      int nodedof = madof[procNode-1];
	      int nnodedof = patch[n]->getNodalDOFs(locNode+1);
	      for (int m = 0;m < nnodedof;m++) {
		PetscInt ieq = meqn[nodedof+m-1];
		if (ieq > 0)
		  subdDofs.insert(ieq-1);
	      }
	    }
	  }

	PetscIntVec sdofs;
	for (std::set<PetscInt>::iterator it = subdDofs.begin(); it != subdDofs.end();it++)
	  sdofs.push_back(*it);
	locSubds.push_back(sdofs);
      }
    }
  }

  return true;
}


bool SAMpatchPara::getLocalSubdomains3D (const IntVec& nxvec,
					 const IntVec& nyvec,
					 const IntVec& nzvec,
					 const IntVec& minNodeId,
					 const IntVec& maxNodeId,
					 PetscIntMat& locSubds) const
{
  if (nxvec.size() != patch.size() ||
      nyvec.size() != patch.size() ||
      nzvec.size() != patch.size())
    return false;

  // Split the patches into smaller subdomains
  for (size_t n = 0; n < patch.size(); n++) {

    int nnod1, nnod2, nnod3;
    if (!static_cast<ASMstruct*>(patch[n])->getSize(nnod1,nnod2,nnod3))
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
	  std::set<PetscInt> subdDofs;

	  int i1 = p1*n1;
	  int i2 = (p1+1)*n1;
	  if ((p1 == nx-1) && (i2 < nnod1)) i2 = nnod1;

	  for (int k = k1;k < k2;k++)
	    for (int j = j1;j < j2;j++)
	      for (int i = i1;i < i2;i++) {
		int locNode = k*nnod1*nnod2 + j*nnod1 + i;
		int procNode = MLGN[locNode];
		int globNode = l2gn[procNode-1];

		if (globNode >= minNodeId[n] && globNode <= maxNodeId[n]) {
		  int nodedof = madof[procNode-1];
		  int nnodedof = patch[n]->getNodalDOFs(locNode+1);
		  for (int m = 0;m < nnodedof;m++) {
		    PetscInt ieq = meqn[nodedof+m-1];
		    if (ieq > 0)
		      subdDofs.insert(ieq-1);
		  }
		}
	      }

	  PetscIntVec sdofs;
	  for (std::set<PetscInt>::iterator it = subdDofs.begin(); it != subdDofs.end();it++)
	    sdofs.push_back(*it);
	  locSubds.push_back(sdofs);
	}
      }
    }
  }
  
  return true;
}


bool SAMpatchPara::getSubdomains1D (const IntVec& nxvec, int overlap,
				    PetscIntMat& subds) const
{
  if (nxvec.size() != patch.size())
    return false;

  // Overlap
  int olow  = overlap/2;
  int ohigh = overlap/2 + overlap%2; 
    
  // Split the patches into smaller subdomains
  for (size_t n = 0; n < patch.size(); n++) {

    int nnod, n2, n3;
    if (!static_cast<ASMstruct*>(patch[n])->getSize(nnod,n2,n3) || n2 || n3)
      return false;

    int nx = nxvec[n];
    int n1 = nnod/nx;
    
    const IntVec& MLGN = patch[n]->getGlobalNodeNums();
    for (int p1 = 0;p1 < nx;p1++) {
      std::set<PetscInt> subdDofs;
      
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
	    subdDofs.insert(ieq-1);
	}
      }
    
      PetscIntVec sdofs;
      for (std::set<PetscInt>::iterator it = subdDofs.begin(); it != subdDofs.end();it++)
	sdofs.push_back(*it);
      subds.push_back(sdofs);
    }
  }

  return true;
}


bool SAMpatchPara::getSubdomains2D (const IntVec& nxvec,
				    const IntVec& nyvec, int overlap,
				    PetscIntMat& subds) const
{
  if (nxvec.size() != patch.size() || nyvec.size() != patch.size())
    return false;

  // Overlap
  int olow  = overlap/2;
  int ohigh = overlap/2 + overlap%2; 

  // Split the patches into smaller subdomains
  for (size_t n = 0; n < patch.size(); n++) {

    int nnod1, nnod2, n3;
    if (!static_cast<ASMstruct*>(patch[n])->getSize(nnod1,nnod2,n3) || n3)
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
	std::set<PetscInt> subdDofs;
	
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
		subdDofs.insert(ieq-1);
	    }
	  }

	PetscIntVec sdofs;
	for (std::set<PetscInt>::iterator it = subdDofs.begin(); it != subdDofs.end();it++)
	  sdofs.push_back(*it);
	subds.push_back(sdofs);
      }
    }
  }
   
  return true;
}


bool SAMpatchPara::getSubdomains3D (const IntVec& nxvec,
				    const IntVec& nyvec,
				    const IntVec& nzvec, int overlap,
				    PetscIntMat& subds) const
{
  if (nxvec.size() != patch.size() ||
      nyvec.size() != patch.size() ||
      nzvec.size() != patch.size())
    return false;

  // Overlap
  int olow  = overlap/2;
  int ohigh = overlap/2 + overlap%2; 

  // Split the patches into smaller subdomains
  for (size_t n = 0; n < patch.size(); n++) {

    int nnod1, nnod2, nnod3;
    if (!static_cast<ASMstruct*>(patch[n])->getSize(nnod1,nnod2,nnod3))
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
	  std::set<PetscInt> subdDofs;

	  int imin = 0;
	  int i1 = std::max(p1*n1-olow,imin);
	  int i2 = std::min((p1+1)*n1+ohigh,nnod1);
	  if ((p1 == nx-1) && (i2 < nnod1)) i2 = nnod1;

	  for (int k = k1;k < k2;k++)
	    for (int j = j1;j < j2;j++)
	      for (int i = i1;i < i2;i++) {
		int locNode = k*nnod1*nnod2 + j*nnod1 + i;
		int globNode = MLGN[locNode];
		
		int nodedof = madof[globNode-1];
		int nnodedof = patch[n]->getNodalDOFs(locNode+1);
		for (int m = 0;m < nnodedof;m++) {
		  PetscInt ieq = meqn[nodedof+m-1];
		  if (ieq > 0)
		    subdDofs.insert(ieq-1);
		}
	      }
	  
	  PetscIntVec sdofs;
	  for (std::set<PetscInt>::iterator it = subdDofs.begin(); it != subdDofs.end();it++)
	    sdofs.push_back(*it);
	  subds.push_back(sdofs);
	}
      }
    }
  }
  
  return true;
}


bool SAMpatchPara::getLocalSubdomains1D (PetscIntMat& locSubds,
					 const IntVec& nxvec,
					 const IntVec& minNodeId,
					 const IntVec& maxNodeId,
					 int f1, int f2) const
{
  if (nxvec.size() != patch.size())
    return false;

  // Split the patches into smaller subdomains
  for (size_t n = 0; n < patch.size(); n++) {

    int nnod, n2, n3;
    if (!static_cast<ASMstruct*>(patch[n])->getSize(nnod,n2,n3) || n2 || n3)
      return false;

    int nx = nxvec[n];
    int n1 = nnod/nx;
    
    int nfield   = f2-f1+1;

    const IntVec& MLGN = patch[n]->getGlobalNodeNums();
    for (int p1 = 0;p1 < nx;p1++) {
      std::set<PetscInt> subdDofs;
      
      int i1 = p1*n1;
      int i2 = (p1+1)*n1;
      if ((p1 == nx-1) && (i2 < nnod)) i2 = nnod;
      
      for (int i = i1;i < i2;i++) {
	int procNode = MLGN[i];
	int globNode = l2gn[procNode-1];

	if (globNode >= minNodeId[n] && globNode <= maxNodeId[n]) {
	  int nodedof   = madof[procNode-1];
	  int gdof      = meqn[nodedof-1];
	  int nnodedof  = patch[n]->getNodalDOFs(i+1);
	  int gbnodedof = (gdof-1)/nnodedof;
	  for (int m = 0;m < nfield;m++) {
	    PetscInt ieq = gbnodedof*nfield+m;
	    if (ieq > -1)
	      subdDofs.insert(ieq);
	  }
	}
      }

      PetscIntVec sdofs;
      for (std::set<PetscInt>::iterator it = subdDofs.begin(); it != subdDofs.end();it++)
	sdofs.push_back(*it);
      locSubds.push_back(sdofs);
    }
  }
  
  return true;
}


bool SAMpatchPara::getLocalSubdomains2D (PetscIntMat& locSubds,
					 const IntVec& nxvec,
					 const IntVec& nyvec,
					 const IntVec& minNodeId,
					 const IntVec& maxNodeId,
					 int f1, int f2) const
{
  if (nxvec.size() != patch.size() || nyvec.size() != patch.size())
    return false;

  // Split the patches into smaller subdomains
  for (size_t n = 0; n < patch.size(); n++) {

    int nnod1, nnod2, n3;
    if (!static_cast<ASMstruct*>(patch[n])->getSize(nnod1,nnod2,n3) || n3)
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
	std::set<PetscInt> subdDofs;
	
	int i1 = p1*n1;
	int i2 = (p1+1)*n1;
	if ((p1 == nx-1) && (i2 < nnod1)) i2 = nnod1;

	for (int j = j1;j < j2;j++)
	  for (int i = i1;i < i2;i++) {
	    int locNode = j*nnod1 + i;
	    int procNode = MLGN[locNode];
	    int globNode = l2gn[procNode-1];

	    if (globNode >= minNodeId[n] && globNode <= maxNodeId[n]) {
	      int nodedof  = madof[procNode-1];
	      int gdof      = meqn[nodedof-1];
	      int nnodedof  = patch[n]->getNodalDOFs(locNode+1);
	      int gbnodedof = (gdof-1)/nnodedof;
	      for (int m = 0;m < nfield;m++) {
		PetscInt ieq = gbnodedof*nfield+m;
		if (ieq > -1)
		  subdDofs.insert(ieq);
	      }
	    }
	  }

	PetscIntVec sdofs;
	for (std::set<PetscInt>::iterator it = subdDofs.begin(); it != subdDofs.end();it++)
	  sdofs.push_back(*it);
	locSubds.push_back(sdofs);

      }
    }
  }

  return true;
}


bool SAMpatchPara::getLocalSubdomains3D (PetscIntMat& locSubds,
					 const IntVec& nxvec,
					 const IntVec& nyvec,
					 const IntVec& nzvec,
					 const IntVec& minNodeId,
					 const IntVec& maxNodeId,
					 int f1, int f2) const
{
  if (nxvec.size() != patch.size() ||
      nyvec.size() != patch.size() ||
      nzvec.size() != patch.size())
    return false;

  // Split the patches into smaller subdomains
  for (size_t n = 0; n < patch.size(); n++) {

    int nnod1, nnod2, nnod3;
    if (!static_cast<ASMstruct*>(patch[n])->getSize(nnod1,nnod2,nnod3))
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
	  std::set<PetscInt> subdDofs;

	  int i1 = p1*n1;
	  int i2 = (p1+1)*n1;
	  if ((p1 == nx-1) && (i2 < nnod1)) i2 = nnod1;

	  for (int k = k1;k < k2;k++)
	    for (int j = j1;j < j2;j++)
	      for (int i = i1;i < i2;i++) {
		int locNode = k*nnod1*nnod2 + j*nnod1 + i;
		int procNode = MLGN[locNode];
		int globNode = l2gn[procNode-1];

		if (globNode >= minNodeId[n] && globNode <= maxNodeId[n]) {
		  int nodedof = madof[procNode-1];
		  int gdof      = meqn[nodedof-1];
		  int nnodedof = patch[n]->getNodalDOFs(locNode+1);
		  int gbnodedof = (gdof-1)/nnodedof;
		  for (int m = 0;m < nfield;m++) {
		    PetscInt ieq = gbnodedof*nfield+m;
		    if (ieq > -1)
		      subdDofs.insert(ieq);
		  }
		}
	      }

	  PetscIntVec sdofs;
	  for (std::set<PetscInt>::iterator it = subdDofs.begin();it != subdDofs.end();it++)
	    sdofs.push_back(*it);
	  locSubds.push_back(sdofs);
	}
      }
    }
  }

  return true;
}


bool SAMpatchPara::getSubdomains1D (PetscIntMat& subds,
				    const IntVec& nxvec,
				    int overlap, int f1, int f2) const
{
  if (nxvec.size() != patch.size())
    return false;

  // Overlap
  int olow  = overlap/2;
  int ohigh = overlap/2 + overlap%2; 

  // Split the patches into smaller subdomains
  for (size_t n = 0; n < patch.size(); n++) {

    int nnod, n2, n3;
    if (!static_cast<ASMstruct*>(patch[n])->getSize(nnod,n2,n3) || n2 || n3)
      return false;

    int nx = nxvec[n];
    int n1 = nnod/nx;
    
    int nfield = f2-f1+1;

    const IntVec& MLGN = patch[n]->getGlobalNodeNums();
    for (int p1 = 0;p1 < nx;p1++) {
      std::set<PetscInt> subdDofs;
      
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
	    subdDofs.insert(ieq);
	}
      }

      PetscIntVec sdofs;
      for (std::set<PetscInt>::iterator it = subdDofs.begin();it != subdDofs.end();it++)
	sdofs.push_back(*it);
      subds.push_back(sdofs);
    }
  }

  return true;
}


bool SAMpatchPara::getSubdomains2D (PetscIntMat& subds,
				    const IntVec& nxvec,
				    const IntVec& nyvec,
				    int overlap, int f1, int f2) const
{
  if (nxvec.size() != patch.size() || nyvec.size() != patch.size())
    return false;

  // Overlap
  int olow  = overlap/2;
  int ohigh = overlap/2 + overlap%2; 

  // Split the patches into smaller subdomains
  for (size_t n = 0; n < patch.size(); n++) {

    int nnod1, nnod2, n3;
    if (!static_cast<ASMstruct*>(patch[n])->getSize(nnod1,nnod2,n3) || n3)
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
	std::set<PetscInt> subdDofs;
	
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
		subdDofs.insert(ieq);
	    }
	  }

	PetscIntVec sdofs;
	for (std::set<PetscInt>::iterator it = subdDofs.begin();it != subdDofs.end();it++)
	  sdofs.push_back(*it);
	subds.push_back(sdofs);
      }
    }
    
  }

  return true;
}


bool SAMpatchPara::getSubdomains3D (PetscIntMat& subds,
				    const IntVec& nxvec,
				    const IntVec& nyvec,
				    const IntVec& nzvec,
				    int overlap, int f1, int f2) const
{
  if (nxvec.size() != patch.size() ||
      nyvec.size() != patch.size() ||
      nzvec.size() != patch.size())
    return false;

  // Overlap
  int olow  = overlap/2;
  int ohigh = overlap/2 + overlap%2; 

  // Split the patches into smaller subdomains
  for (size_t n = 0; n < patch.size(); n++) {

    int nnod1, nnod2, nnod3;
    if (!static_cast<ASMstruct*>(patch[n])->getSize(nnod1,nnod2,nnod3))
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
	  std::set<PetscInt> subdDofs;

	  int imin = 0;
	  int i1 = std::max(p1*n1-olow,imin);
	  int i2 = std::min((p1+1)*n1+ohigh,nnod1);
	  if ((p1 == nx-1) && (i2 < nnod1)) i2 = nnod1;

	  for (int k = k1;k < k2;k++)
	    for (int j = j1;j < j2;j++)
	      for (int i = i1;i < i2;i++) {
		int locNode = k*nnod1*nnod2 + j*nnod1 + i;
		int globNode = MLGN[locNode];
		
		int nodedof = madof[globNode-1];
		int nnodedof = patch[n]->getNodalDOFs(locNode+1);

		int gdof      = meqn[nodedof-1];
		int gbnodedof = (gdof-1)/nnodedof;

		for (int m = 0;m < nfield;m++) {
		  PetscInt ieq = gbnodedof*nfield+m;
		  if (ieq > -1)
		    subdDofs.insert(ieq);
		}
	      }

	  PetscIntVec sdofs;
	  for (std::set<PetscInt>::iterator it = subdDofs.begin();it != subdDofs.end();it++)
	    sdofs.push_back(*it);
	  subds.push_back(sdofs);
	}
      }
    }
  }

  return true;
}


int SAMpatchPara::getLocalNodeCoordinates (PetscRealVec& coords) const
{
  // Find min and max node for each patch on this processor
  IntVec maxNodeId, minNodeId;
  if (!this->getMinMaxNode(minNodeId,maxNodeId))
    return 0;

  int min     = minNodeId.front();
  int nlocnod = maxNodeId.back()-min+1;
  int nsd     = patch.front()->getNoSpaceDim();
  coords.resize(nlocnod*nsd,0.0);

  for (size_t n = 0; n < patch.size(); n++) {
    const IntVec& MLGN = patch[n]->getGlobalNodeNums();
    for (size_t i = 0;i < MLGN.size();i++) {
      int gnod = l2gn[MLGN[i]-1];
      if ((gnod >= minNodeId[n]) && (gnod <= maxNodeId[n]))
      {
	Vec3 X = patch[n]->getCoord(i+1);
	size_t locnode = gnod-min+1;
	for (int k = 0;k < nsd;k++)
	  coords[(locnode-1)*nsd + k] = X[k];
      }
    }
  }

  return nsd;
}


bool SAMpatchPara::getMinMaxNode(IntVec& minNodeId, IntVec& maxNodeId) const
{
  // Find min and max node for each patch on this processor
  size_t i, n, npatch = patch.size();
  maxNodeId.resize(npatch,true);
  minNodeId.resize(npatch,true);

  for (n = 0; n < npatch; n++) {
    const IntVec& MLGN = patch[n]->getGlobalNodeNums();
    int min = 1e9;
    int max = 0;
    int ieq;
    for (i = 0; i < MLGN.size(); i++) {
      ieq = l2gn[MLGN[i]-1];
      if (ieq > max)
	max = ieq;
      if (ieq < min)
	min = ieq;
    }
    minNodeId[n] = min;
    maxNodeId[n] = max;
  }

  minNodeId[0] = 1;
  for (n = 1; n < npatch; n++)
    minNodeId[n] = maxNodeId[n-1] + 1;

#ifdef HAS_PETSC
  if (adm.isParallel()) {
    int myRank = adm.getProcId();
    if (myRank < nProc-1)
      adm.send(maxNodeId[npatch-1],myRank+1);
    if (myRank > 0) {
      adm.receive(minNodeId[0],myRank-1);
      minNodeId[0]++;
    }
  }
#endif

  return true;
}


bool SAMpatchPara::isDirichlet(int inod, int dof) const
{
  std::vector<int> idofs;

  int nnodedof = madof[inod]-madof[inod-1];

  int idof = dof;
  int rdof;
  if (rdof = idof/100) {
    if (rdof <= nnodedof)
      idofs.push_back(rdof);
    idof = idof%100;
  }
  if (rdof = idof/10) {
    if (rdof <= nnodedof)
      idofs.push_back(rdof);
    idof =idof%10;
  }
  if (idof <= nnodedof)
    idofs.push_back(idof);

  bool fixed = false;
  for (size_t i = 0;i < idofs.size();i++) {
    idof = idofs[i];
    if (mpmceq[madof[inod-1]+idof-2] != 0)
      fixed = true;
  }

  return fixed;
}
