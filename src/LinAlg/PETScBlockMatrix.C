// $Id$
//==============================================================================
//!
//! \file PETScMatrix.C
//!
//! \date Oct 20 2012
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Representation of the system matrix in PETSc block format.
//!
//==============================================================================

#include "PETScBlockMatrix.h"
#include "ProcessAdm.h"
#ifdef HAS_PETSC
#include "LinSolParams.h"
#include "LinAlgInit.h"
#include "SAMpatchPara.h"
#include "PCPerm.h"
#include "PCScale.h"
#include <sstream>

#ifdef USE_OPENMP
#include <omp.h>
#endif


PETScBlockMatrix::PETScBlockMatrix (const ProcessAdm& padm, const LinSolParams& par) : PETScMatrix(padm,par)
{
  ncomps.resize(par.ncomps.size());
  nblocks = ncomps.size();
  for (size_t i = 0;i < nblocks;i++)
    ncomps[i] = par.ncomps[i];
  
  isvec = (IS*) malloc(sizeof(IS)*nblocks);
  matvec = (Mat*) malloc(sizeof(Mat)*nblocks*nblocks);
  for (size_t m = 0;m < nblocks*nblocks;m++) {
    MatCreate(*adm.getCommunicator(),&matvec[m]);
    MatSetOption(matvec[m],MAT_STRUCTURALLY_SYMMETRIC,PETSC_TRUE);
    MatSetFromOptions(matvec[m]);
  }
  
  S  = NULL;
  Sp = NULL;
  Fp = NULL;

  LinAlgInit::increfs();
}


PETScBlockMatrix::PETScBlockMatrix (const ProcessAdm& padm, const IntVec& nc, const LinSolParams& par)
  : PETScMatrix(padm,par)
{
  nblocks = nc.size();
  ncomps.resize(nblocks);
  for (size_t i = 0;i < nblocks;i++)
    ncomps[i]  = nc[i];

  isvec = (IS*) malloc(sizeof(IS)*nblocks);
  matvec = (Mat*) malloc(sizeof(Mat)*nblocks*nblocks);
  for (size_t m = 0;m < nblocks*nblocks;m++) {
    MatCreate(*adm.getCommunicator(),&matvec[m]);
    MatSetFromOptions(matvec[m]);
  }

  LinAlgInit::increfs();
}


PETScBlockMatrix::~PETScBlockMatrix ()
{
  // Deallocation of null space
  if (solParams.getNullSpace() == CONSTANT)
    MatNullSpaceDestroy(PETSCMANGLE(nsp));

  // Deallocation of linear solver object.
  KSPDestroy(PETSCMANGLE(ksp));

  for (size_t m = 0;m < nblocks*nblocks;m++)
    MatDestroy(&matvec[m]);
  free(matvec);
  for (size_t m = 0;m < nblocks;m++)
    ISDestroy(&isvec[m]);
  free(isvec);

  // Deallocation of matrix object.
  MatDestroy(PETSCMANGLE(A));
  LinAlgInit::decrefs();

  // Deallocation of Schur matrix
  MatDestroy(PETSCMANGLE(Sp));

  for (PetscInt i = 0; i < ISsize; i++)
    ISDestroy(PETSCMANGLE(elmIS[i]));

  delete elmIS;
}


void PETScBlockMatrix::assemPETSc (const Matrix& eM, PETScVector& SV,
				   const std::vector<int>& meen, const int* meqn,
				   const int* mpmceq, const int* mmceq,
				   const Real* ttcc)
{
  Real   c0;
  size_t i, j;
  int    jp, jceq;

  // Number of degrees of freedom for element
  size_t nedof = meen.size();

  // Convert meen to 0-based C array
  PetscIntVec l2g(nedof);
  for (i = 0; i < nedof; i++)
    l2g[i] = meqn[meen[i]-1]-1;

  // Cast to non-constant Matrix to modify for Dirichlet BCs
  Matrix& Amat = const_cast<Matrix&>(eM);

  std::vector<Real> uc(nedof), bc(nedof);

  bool rhsMod = false;
  for (j = 1;j <= nedof;j++) {
    uc[j-1] = bc[j-1] = 0.0;
    jceq = mpmceq[meen[j-1]-1];
    if (jceq < 1) continue;

    jp = jceq-1;
    c0 = ttcc[jp];

    if (c0 != 0.0) {
      uc[j-1] = -c0;
      rhsMod = true;
    }
  }

  if (rhsMod)
    Amat.multiply(uc,bc);

  // Eliminate constrained dofs from element matrix
  for (j = 1;j <= nedof;j++) {
    jceq = mpmceq[meen[j-1]-1];
    if (jceq < 1) continue;

    bc[j-1] = -uc[j-1];
    for (i = 1;i <= nedof;i++)
      Amat(i,j) = Amat(j,i) = 0.0;
    Amat(j,j) = 1.0;
  }

  // Add contributions to SV (righthand side)
  VecSetValues(SV.getVector(),nedof,&(l2g[0]),&(bc.front()),ADD_VALUES);

  std::vector< std::vector<Matrix> > eMb;
  PetscIntMat l2gb;
  this->getBlockElmMatData(eM,l2g,eMb,l2gb);

  int m = 0;
  for (i = 0;i < nblocks;i++)
    for (j = 0;j < nblocks;j++, m++) {
      (eMb[i][j]).transpose();
      MatSetValues(matvec[m],l2gb[i].size(),&(l2gb[i][0]),l2gb[j].size(),&(l2gb[j][0]),(eMb[i][j]).ptr(),ADD_VALUES);
    }
}


void PETScBlockMatrix::getBlockElmMatData(const Matrix& Amat,
					  const PetscIntVec& l2g,
					  std::vector< std::vector<Matrix> >& Ab,
					  PetscIntMat& l2gb) const
{
  int nf = 0;
  for (size_t i = 0;i < nblocks;i++)
    nf += ncomps[i];

  int nenod = Amat.rows()/nf;
  l2gb.resize(nblocks);
  int stride = 0;
  for (size_t i = 0;i < nblocks;i++) {
    l2gb[i].resize(nenod*ncomps[i]);

    for (int j = 0;j < nenod;j++)
      for (int k = 0;k < ncomps[i];k++) {
	int dof = l2g[j*nf + stride + k];
	int node = dof/nf;
	l2gb[i][j*ncomps[i] + k] = node*ncomps[i] + k;
      }

    stride += ncomps[i];
  }
  
  Ab.resize(nblocks);
  for (size_t i = 0;i < nblocks;i++)
    Ab[i].resize(nblocks);
  
  int istride = 0;
  for (size_t i = 0;i < nblocks;i++) {
    int jstride = 0;
    for (size_t j = 0;j < nblocks;j++) {
      Ab[i][j].resize(nenod*ncomps[i],nenod*ncomps[j]);
      for (int m = 1;m <= nenod;m++) 
	for (int n = 1;n <= nenod;n++) 
	  for (int r = 1;r <= ncomps[i];r++) 
	    for (int s = 1;s <= ncomps[j];s++) 
	      (Ab[i][j])((m-1)*ncomps[i]+r,(n-1)*ncomps[j]+s) = Amat((m-1)*nf + istride + r,(n-1)*nf + jstride + s);
      jstride += ncomps[j];
    }
    istride += ncomps[i];
  }
}


void PETScBlockMatrix::initAssembly (const SAM& sam, bool)
{
  const SAMpatchPara* sampch = dynamic_cast<const SAMpatchPara*>(&sam);

  bool getcoords = false;
  if (!strncasecmp(solParams.prec.c_str(),"gamg",4) || !strncasecmp(solParams.prec.c_str(),"ml",2))
    getcoords = true;
  else
    for (size_t i = 0;i < solParams.subprec.size();i++)
      if (!strncasecmp(solParams.subprec[i].c_str(),"gamg",4) || !strncasecmp(solParams.subprec[i].c_str(),"ml",2))
	getcoords = true;
  if (getcoords)
    sampch->getLocalNodeCoordinates(coords);

  if (!solParams.dirOrder.empty()) {
    dirIndexSet.resize(solParams.dirOrder.size());
    for (size_t i = 0;i < solParams.dirOrder.size();i++)
      for (size_t j = 0;j < solParams.dirOrder[i].size();j++) {
	IS permIndex;
	if (solParams.dirOrder[i][j] != 123) {
	  PetscIntVec perm;
	  sampch->getDirOrdering(perm,solParams.dirOrder[i][j],solParams.ncomps[i]);
	  ISCreateGeneral(*adm.getCommunicator(),perm.size(),&(perm[0]),PETSC_COPY_VALUES,&permIndex);
	  ISSetPermutation(permIndex);
	}
	else {
	  ISCreate(*adm.getCommunicator(),&permIndex);
	  ISSetIdentity(permIndex);
	}
	dirIndexSet[i].push_back(permIndex);
      }
  }

  int nx = solParams.getLocalPartitioning(0);
  int ny = solParams.getLocalPartitioning(1);
  int nz = solParams.getLocalPartitioning(2);
  int overlap = solParams.getOverlap();

  if (nx+ny+nz > 0) {
    locSubdDofsBlock.resize(nblocks);
    subdDofsBlock.resize(nblocks);
    int f1 = 0;
    int f2 = -1;
    for (size_t k = 0;k < nblocks;k++) {
      f2 += ncomps[k];
      sampch->getLocalSubdomainsBlock(locSubdDofsBlock[k],f1,f2,nx,ny,nz);
      sampch->getSubdomainsBlock(subdDofsBlock[k],f1,f2,overlap,nx,ny,nz);
      f1 = f2+1;
    }
  }

  // Get number of equations in linear system
  const PetscInt neq = sam.getNoEquations();

  // Get number of fields and number of nodes
  size_t nf = 0;
  size_t ncomp = ncomps.size();
  for (size_t i = 0;i < ncomp;i++)
    nf += ncomps[i];
  size_t nnod = neq/nf;

  // Set size of matrices
   int k = 0;
  for (size_t m = 0;m < nblocks;m++)
    for (size_t n = 0;n < nblocks;n++, k++) {
      MatSetSizes(matvec[k],nnod*ncomps[m],nnod*ncomps[n],PETSC_DETERMINE,PETSC_DETERMINE);
      //MatSetBlockSize(matvec[k],ncomps[n]);
      MatSetFromOptions(matvec[k]);
    }

  // Equation numbers 
  int ifirst = sampch->getMinEqNumber();
  ifirst--;
  int ilast  = sampch->getMaxEqNumber();
  
  // Allocation of sparsity pattern
  std::vector<std::vector<IntVec> > d_nnz, o_nnz;  
  if (adm.isParallel()) {
    int myRank = adm.getProcId(); 
    int nProc  = adm.getNoProcs();
    if (myRank < nProc-1)
      adm.send(ilast,myRank+1);
    if (myRank > 0) {
      adm.receive(ifirst,myRank-1);
    }
    
    if (sampch->getNoDofCouplings(ifirst,ilast,ncomps,d_nnz,o_nnz)) {
      PetscIntVec d_Nnz;
      PetscIntVec o_Nnz;
      int id = 0;
      for (size_t m = 0;m < nblocks;m++)
	for (size_t n = 0;n < nblocks;n++) {
	  size_t dsize = d_nnz[m][n].size();
	  d_Nnz.resize(dsize);	
	  for (size_t i = 0; i < dsize; i++)
	    d_Nnz[i] = d_nnz[m][n][i];
	  
	  size_t osize = o_nnz[m][n].size();
	  o_Nnz.resize(osize);
	  for (size_t i = 0; i < osize; i++)
	    o_Nnz[i] = o_nnz[m][n][i];
	  
	  MatMPIAIJSetPreallocation(matvec[id++],PETSC_DEFAULT,&(d_Nnz[0]),PETSC_DEFAULT,&(o_Nnz[0]));
	}
    }
  }
  else {
    if (sampch->getNoDofCouplings(ifirst,ilast,ncomps,d_nnz,o_nnz)) {
      PetscIntVec d_Nnz;
      int id = 0;
      for (size_t m = 0;m < nblocks;m++) {
	for (size_t n = 0;n < nblocks;n++) {
	  size_t dsize = d_nnz[m][n].size();
	  d_Nnz.resize(dsize);	
	  for (size_t i = 0; i < dsize; i++)
	    d_Nnz[i] = d_nnz[m][n][i];
	  
	  MatSeqAIJSetPreallocation(matvec[id++],PETSC_DEFAULT,&(d_Nnz[0]));
	}	  
      }
    }
  }

  // Initialize the block matrix
  PetscIntVec  idx;

  int gidx = ifirst;
  for (size_t i = 0;i < ncomps.size();i++) {
    PetscInt ncomp  = ncomps[i];
    PetscInt nentry = ncomp*nnod;

    PetscInt lidx = 0; 
    idx.resize(nentry);
    for (size_t n = 0;n < nnod;n++) {
      for (int k = 0;k < ncomp;k++)
	idx[lidx++] = gidx++; 
    }

    ISCreateGeneral(*adm.getCommunicator(),nentry,&(idx[0]),PETSC_COPY_VALUES,&(isvec[i]));
  }

  MatCreateNest(*adm.getCommunicator(),nblocks,isvec,nblocks,isvec,matvec,&A);

#ifndef SP_DEBUG
  // Do not abort program for allocation error in release mode
  for (size_t i = 0;i < nblocks*nblocks;i++)
    MatSetOption(matvec[i],MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
#endif

  if (solParams.getNullSpace(1) == CONSTANT) {
     Vec const_mode;
     VecCreate(*adm.getCommunicator(), &const_mode);
     VecSetSizes(const_mode, ilast-ifirst, PETSC_DECIDE);
     VecSetFromOptions(const_mode);
     std::vector<PetscInt> ix((ilast-ifirst)/4);
     std::vector<PetscScalar> y(ix.size(), 1.0);
     for (size_t i=0;i<ix.size();++i)
       ix[i] = ifirst+i*4+3;
     VecSetValues(const_mode, ix.size(), &ix[0], &y[0], INSERT_VALUES);
     VecAssemblyBegin(const_mode);
     VecAssemblyEnd(const_mode);
     Vec x;
     VecDuplicate(const_mode,&x);
     this->renumberRHS(const_mode,x,true);
     VecNormalize(x, NULL);
     MatNullSpaceCreate(*adm.getCommunicator(),PETSC_FALSE,1,&x,&nsp);
     KSPSetNullSpace(ksp,nsp);
     VecDestroy(&x);
     VecDestroy(&const_mode);
  }
}


void PETScBlockMatrix::init ()
{
  // Set all matrix elements to zero for each submatrix
  for (size_t m = 0;m < nblocks*nblocks;m++) 
    MatZeroEntries(matvec[m]);

  // Set all matrix elements to zero
  MatZeroEntries(A);
}


bool PETScBlockMatrix::assemble (const Matrix& eM, const SAM& sam,
                                 SystemVector& B, int e)
{
  // Get mapping "meen" between local degrees of freedom in element e
  // and global degrees of freedom.
  std::vector<int> meen;
  if (!sam.getElmEqns(meen,e,eM.rows()))
    return false;
  
  PETScVector* Bptr = dynamic_cast<PETScVector*>(&B);
  if (!Bptr)
    return false;

  // Assemble local stiffness matrix into global system.
  this->assemPETSc(eM,*Bptr,meen,sam.meqn,sam.mpmceq,sam.mmceq,sam.ttcc);

  return true;
}


bool PETScBlockMatrix::solve (SystemVector& B, bool newLHS)
{
  // Reset linear solver
  if (nLinSolves && solParams.nResetSolver)
    if (nLinSolves%solParams.nResetSolver == 0) {
      KSPDestroy(&ksp);
      KSPCreate(*adm.getCommunicator(),&ksp);
      if (solParams.schur)
	MatDestroy(&Sp);
      setParams = true;
    }

  PETScVector* Bptr = dynamic_cast<PETScVector*>(&B);
  if (!Bptr)
    return false;

  Vec x;
  VecDuplicate(Bptr->getVector(),&x);
  VecCopy(Bptr->getVector(),x);
  this->renumberRHS(Bptr->getVector(),x,true);

  // Has lefthand side changed?
  static int firstIt = true;
  if (firstIt)
    //if (newLHS)
#if PETSC_VERSION_MINOR < 5
    KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);
#else
    KSPSetOperators(ksp,A,A);
#endif
  else
#if PETSC_VERSION_MINOR < 5
    KSPSetOperators(ksp,A,A,SAME_PRECONDITIONER);
#else
    KSPSetOperators(ksp,A,A);
#endif
  firstIt = false;

  if (setParams) {
    this->setParameters();
    setParams = false;
  }
  KSPSetInitialGuessKnoll(ksp,PETSC_TRUE);
  KSPSolve(ksp,x,Bptr->getVector());

  // Renumber back to usual numbering
  this->renumberRHS(Bptr->getVector(),Bptr->getVector(),false);

  PetscInt its;
  KSPGetIterationNumber(ksp,&its);
  PetscPrintf(PETSC_COMM_WORLD,"\n Iterations for %s = %D\n",solParams.getMethod(),its);
  VecDestroy(PETSCMANGLE(x));

  nIts += its;
  nLinSolves++;

  return true;
}


bool PETScBlockMatrix::solve (const SystemVector& b, SystemVector& x, bool newLHS)
{
  // Reset linear solver
  if (solParams.nResetSolver)
    if (nLinSolves%solParams.nResetSolver == 0) {
      KSPDestroy(&ksp);
      KSPCreate(*adm.getCommunicator(),&ksp);
      setParams = true;
    }
  
  SystemVector* Bp = const_cast<SystemVector*>(&b);
  PETScVector* Bptr = dynamic_cast<PETScVector*>(Bp);
  if (!Bptr)
    return false;
  PETScVector* Xptr = dynamic_cast<PETScVector*>(&x);
  if (!Xptr)
    return false;

  // Renumber RHS blockwise
  this->renumberRHS(Bptr->getVector(),Bptr->getVector(),true);

  // Has lefthand side changed?
  if (newLHS)
#if PETSC_VERSION_MINOR < 5
    KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);
#else
    KSPSetOperators(ksp,A,A);
#endif
  else
#if PETSC_VERSION_MINOR < 5
    KSPSetOperators(ksp,A,A,SAME_PRECONDITIONER);
#else
    KSPSetOperators(ksp,A,A);
#endif

  if (setParams) {
    this->setParameters();
    setParams = false;
  }

  KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
  KSPSolve(ksp,Bptr->getVector(),Xptr->getVector());

  // Renumber back 
  this->renumberRHS(Xptr->getVector(),Xptr->getVector(),false);

  PetscInt its;
  KSPGetIterationNumber(ksp,&its);
  PetscPrintf(PETSC_COMM_WORLD,"\n Iterations for %s = %D\n",solParams.getMethod(),its);

  nIts += its;
  nLinSolves++;

  return true;
}


bool PETScBlockMatrix::solve (SystemVector& B, SystemMatrix& P, SystemVector& Pb, bool newLHS)
{
  // Reset linear solver
  if (solParams.nResetSolver)
    if (nLinSolves%solParams.nResetSolver == 0) {
      KSPDestroy(&ksp);
      KSPCreate(*adm.getCommunicator(),&ksp);
      setParams = true;
    }

  PETScVector* Bptr = dynamic_cast<PETScVector*>(&B);
  if (!Bptr)
    return false;

  PETScBlockMatrix* Pptr = dynamic_cast<PETScBlockMatrix*>(&P);
  if (!Pptr)
    return false;

  PETScVector* Pbptr = dynamic_cast<PETScVector*>(&Pb);
  if (!Pbptr)
    return false;

  Vec x;
  VecDuplicate(Bptr->getVector(),&x);
  VecCopy(Bptr->getVector(),x);
  this->renumberRHS(Bptr->getVector(),x,true);

  // Has lefthand side changed?
  static int firstIt = true;
  //if (newLHS)
  if (firstIt)
#if PETSC_VERSION_MINOR < 5
    KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);
#else
    KSPSetOperators(ksp,A,A);
#endif
  //else
  //  KSPSetOperators(ksp,A,A,SAME_PRECONDITIONER);
  firstIt = false;

  if (setParams) {
    if (!this->setParameters(Pptr,Pbptr))
      return false;
    setParams = false;
  }
  KSPSetInitialGuessKnoll(ksp,PETSC_TRUE);
  KSPSolve(ksp,x,Bptr->getVector());

  // Renumber back to usual numbering
  this->renumberRHS(Bptr->getVector(),Bptr->getVector(),false);

  PetscInt its;
  KSPGetIterationNumber(ksp,&its);
  PetscPrintf(PETSC_COMM_WORLD,"\n Iterations for %s = %D\n",solParams.getMethod(),its);
  VecDestroy(PETSCMANGLE(x));

  nIts += its;
  nLinSolves++;

  return true;
}

void PETScBlockMatrix::renumberRHS(const Vec& b, Vec& bnew, bool renum2block)
{
  PetscInt low, high;
  VecGetOwnershipRange(b,&low,&high);


  PetscInt nldof;
  VecGetLocalSize(b,&nldof);
  PetscScalar* bp;
  VecGetArray(b,&bp);

  int nf = 0;
  for (size_t i = 0;i < nblocks;i++)
    nf += ncomps[i];
  int nnod = nldof/nf;

  std::vector<PetscScalar> y;
  PetscIntVec l2g;
  y.resize(nldof);
  l2g.resize(nldof);
  int id = 0;
  int stride = 0;
  int gidx = low;
  if (renum2block)
    for (size_t i = 0;i < nblocks;i++) {
      for (int j = 0;j < nnod;j++)
	for (int k = 0;k < ncomps[i];k++,id++) {
	  y[id] = bp[j*nf + stride + k];
	  l2g[id] = gidx++;
	}
      
      stride += ncomps[i];
    }
  else 
    for (size_t i = 0;i < nblocks;i++) {
      for (int j = 0;j < nnod;j++)
	for (int k = 0;k < ncomps[i];k++,id++) {
	  y[j*nf + stride + k] = bp[id];
	  l2g[id] = gidx++;
	}
      
      stride += ncomps[i];
    }
  
  VecRestoreArray(b,&bp);
  VecSetValues(bnew,nldof,&(l2g[0]),&(y[0]),INSERT_VALUES);
  VecAssemblyBegin(bnew);
  VecAssemblyEnd(bnew);
} 


bool PETScBlockMatrix::setParameters(PETScBlockMatrix *P, PETScVector *Pb)
{
  // Set linear solver method
  KSPSetType(ksp,solParams.method.c_str());
  KSPSetTolerances(ksp,solParams.rtol,solParams.atol,solParams.dtol,solParams.maxIts);

  // Set preconditioner
  PC pc;
  KSPGetPC(ksp,&pc);
  PCSetType(pc,solParams.prec.c_str());
  if (!strncasecmp(solParams.prec.c_str(),"gamg",4) || !strncasecmp(solParams.prec.c_str(),"ml",2)) {
    PetscInt nloc = coords.size()/solParams.nsd;
    PCSetCoordinates(pc,solParams.nsd,nloc,&coords[0]);
  }
  //PCFactorSetMatSolverPackage(pc,solParams.package[0].c_str());
  if (!strncasecmp(solParams.prec.c_str(),"fieldsplit",10)) {
    PetscInt m1, m2, n1, n2, nr, nc, nsplit;
    KSP  *subksp;
    PC   subpc[2];
    Vec diagA00;
    MatGetLocalSize(matvec[0],&m1,&n1);
    MatGetLocalSize(matvec[3],&m2,&n2);
    if (adm.isParallel())
      VecCreateMPI(*adm.getCommunicator(),m1,PETSC_DETERMINE,&diagA00);
    else
      VecCreateSeq(PETSC_COMM_SELF,m1,&diagA00);

    if (!P || (solParams.schurPrec == SIMPLE)) {
      MatGetDiagonal(matvec[0],diagA00);
    }
    else { 
      MatGetDiagonal(P->getMatrixBlock(0,0),diagA00);
    }

    VecReciprocal(diagA00);
    VecScale(diagA00,-1.0);
    MatDiagonalScale(matvec[1],diagA00,PETSC_NULL);
    MatMatMult(matvec[2],matvec[1],MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Sp);
    VecReciprocal(diagA00);
    MatDiagonalScale(matvec[1],diagA00,PETSC_NULL);
    VecDestroy(&diagA00);
    MatGetSize(Sp,&nr,&nc);
    MatAXPY(Sp,1.0,matvec[3],DIFFERENT_NONZERO_PATTERN);
    //MatCreateSchurComplement(matvec[0],matvec[0],matvec[1],matvec[2],matvec[3],&S);
    PCFieldSplitSetIS(pc,"u",isvec[0]);
    PCFieldSplitSetIS(pc,"p",isvec[1]);
    PCFieldSplitSetType(pc,PC_COMPOSITE_SCHUR);
    PCFieldSplitSetSchurFactType(pc,PC_FIELDSPLIT_SCHUR_FACT_UPPER);
    //PCFieldSplitSchurPrecondition(pc,PC_FIELDSPLIT_SCHUR_PRE_USER,Sp);
    PCSetFromOptions(pc);
    PCSetUp(pc);
    PCFieldSplitGetSubKSP(pc,&nsplit,&subksp);
#if PETSC_VERSION_MINOR < 5
    KSPSetOperators(subksp[1],Sp,Sp,SAME_PRECONDITIONER);
#else
    KSPSetOperators(subksp[1],Sp,Sp);
#endif

    if (P && Pb && (solParams.schurPrec == PCD)) {
      KSPSetType(subksp[1],"preonly");
      KSPGetPC(subksp[1],&subpc[1]);
      PCSetType(subpc[1],PCSHELL);
      PCCreate(*adm.getCommunicator(),&S);
      PCCreate(*adm.getCommunicator(),&Fp);

      if (!strncasecmp(solParams.subprec[1].c_str(),"compositedir",12)) {
	if (!solParams.addDirSmoother(S,Sp,1,dirIndexSet))
	  return false;
      }
      else
	PCSetType(S,solParams.subprec[1].c_str());
      PCSetType(Fp,PCMAT);
      Vec Pbr;
      VecDuplicate(Pb->getVector(),&Pbr);
      VecCopy(Pb->getVector(),Pbr);
      this->renumberRHS(Pb->getVector(),Pbr,true);
      Vec QpLvec;
      VecGetSubVector(Pbr,isvec[1],&QpLvec);
      VecDuplicate(QpLvec,&QpL);
      VecCopy(QpLvec,QpL);
      VecReciprocal(QpL);
      //MatDiagonalScale(P->getMatrixBlock(1,1),PETSC_NULL,QpL);
      VecRestoreSubVector(Pbr,isvec[1],&QpLvec);
      VecDestroy(&Pbr);
#if PETSC_VERSION_MINOR < 5
      PCSetOperators(S,Sp,Sp,SAME_PRECONDITIONER);
#else
      PCSetOperators(S,Sp,Sp);
#endif
      
      PCProdCreate(&pcprod);
      PCShellSetApply(subpc[1],PCProdApply);
      PCShellSetContext(subpc[1],pcprod);
      PCShellSetDestroy(subpc[1],PCProdDestroy);
      PCShellSetName(subpc[1],"PCD");

      KSPSetType(subksp[1],"preonly");
      if (!strncasecmp(solParams.subprec[1].c_str(),"compositedir",12)) {
	if (!solParams.addDirSmoother(S,Sp,1,dirIndexSet))
	  return false;
      }
      else
	PCSetType(S,solParams.subprec[1].c_str());
      PCFactorSetLevels(S,solParams.levels[1]);
      if (!strncasecmp(solParams.subprec[1].c_str(),"asm",3)) {
	PCASMSetType(S,PC_ASM_BASIC);
	PCASMSetOverlap(S,solParams.overlap[1]);
	if (!locSubdDofsBlock.empty() && !subdDofsBlock.empty()) {
	  const size_t nsubds = subdDofsBlock[1].size();
	  
	  IS isLocSubdDofs[nsubds], isSubdDofs[nsubds];
	  for (size_t i = 0;i < nsubds;i++) {
	    ISCreateGeneral(PETSC_COMM_SELF,locSubdDofsBlock[1][i].size(),&(locSubdDofsBlock[1][i][0]),PETSC_USE_POINTER,&(isLocSubdDofs[i]));
	    ISCreateGeneral(PETSC_COMM_SELF,subdDofsBlock[1][i].size(),&(subdDofsBlock[1][i][0]),PETSC_USE_POINTER,&(isSubdDofs[i]));
	  }
	  PCASMSetLocalSubdomains(S,nsubds,isSubdDofs,isLocSubdDofs);
	}
	
	if (solParams.asmlu[1]) {
	  KSP* subksp;
	  PC   subsubpc;
	  PetscInt first, nlocal;
	  PCSetUp(S);
	  PCASMGetSubKSP(S,&nlocal,&first,&subksp);
	  
	  for (int i = 0; i < nlocal; i++) {
	    KSPGetPC(subksp[i],&subsubpc);
	    PCSetType(subsubpc,PCLU);
	    KSPSetType(subksp[i],KSPPREONLY);
	  }
	}
      }
      else if (!strncasecmp(solParams.subprec[1].c_str(),"ml",2)) {
	PetscInt n;
        solParams.setMLOptions("fieldsplit_p", 1);

	PCSetFromOptions(S);
	PCSetUp(S);

	// Set coarse solver package
	if (!solParams.MLCoarsePackage.empty()) {
	  KSP cksp;
	  PC  cpc;
	  PCMGGetCoarseSolve(S,&cksp);
	  KSPGetPC(cksp,&cpc);
	  PCFactorSetMatSolverPackage(cpc,solParams.MLCoarsePackage[1].c_str());
	}
	
	PCMGGetLevels(S,&n);
	// Presmoother
	for (int i = 1;i < n;i++) {
	  KSP preksp;
	  PC  prepc;
	  
	  // Set smoother
	  std::string smoother;
	  PetscInt noSmooth;

	  PCMGGetSmoother(S,i,&preksp);
	  KSPSetType(preksp,"richardson");

          if ((i == n-1) && (!solParams.finesmoother.empty())) {
            smoother = solParams.finesmoother[1];
	    noSmooth = solParams.noFineSmooth[1];
	  }
          else {
	    smoother = solParams.presmoother[1];
	    noSmooth = solParams.noPreSmooth[1];
	  }

	  KSPSetTolerances(preksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,noSmooth);
	  KSPGetPC(preksp,&prepc); 

	  if (smoother == "asm" || smoother == "asmlu" ) {
	    PCSetType(prepc,"asm");
	    PCASMSetType(prepc,PC_ASM_BASIC);
	    PCASMSetOverlap(prepc,solParams.overlap[1]);
	    
	    if (!locSubdDofs.empty() && !subdDofs.empty()) {
	      const size_t nsubds = subdDofs.size();
	      
	      IS isLocSubdDofs[nsubds], isSubdDofs[nsubds];
	      for (size_t k = 0;k < nsubds;k++) {
		ISCreateGeneral(PETSC_COMM_SELF,locSubdDofs[k].size(),&(locSubdDofs[k][0]),PETSC_USE_POINTER,&(isLocSubdDofs[k]));
		ISCreateGeneral(PETSC_COMM_SELF,subdDofs[k].size(),&(subdDofs[i][0]),PETSC_USE_POINTER,&(isSubdDofs[k]));
	      }
	      PCASMSetLocalSubdomains(prepc,nsubds,isSubdDofs,isLocSubdDofs);
	    }
	    
	    // If LU factorization is used on each subdomain
	    if (smoother == "asmlu") {
	      KSP* subksp;
	      PC   subpc;
	      PetscInt first, nlocal;
	      PCSetFromOptions(prepc);
	      PCSetUp(prepc);
	      PCASMGetSubKSP(prepc,&nlocal,&first,&subksp);
	      
	      for (int k = 0; k < nlocal; k++) {
		KSPGetPC(subksp[i],&subpc);
		PCSetType(subpc,PCLU);
		KSPSetType(subksp[k],KSPPREONLY);
	      }
	    }
	  }
	  else if (!strncasecmp(smoother.c_str(),"compositeDir",12) && (i==n-1)) {
	    if (!solParams.addDirSmoother(prepc,Sp,1,dirIndexSet))
	      return false;
	  }
	  else
	    PCSetType(prepc,smoother.c_str());
	  
          PCFactorSetLevels(prepc,solParams.levels[1]);
          KSPSetUp(preksp);
	}
      }
      else if (!strncasecmp(solParams.subprec[1].c_str(),"gamg",4) || !strncasecmp(solParams.subprec[1].c_str(),"ml",2) ) {
        PetscInt nloc = coords.size()/solParams.nsd;
        PCSetCoordinates(pc,solParams.nsd,nloc,&coords[0]);
      }

      PCProdSetUp(subpc[1],&QpL,&Fp,&S);
      nsplit--;
    }
      
    // Preconditioner for blocks
    for (PetscInt m = 0; m < nsplit; m++) {
      std::string prefix;
      if (m == 0)
	prefix = "fieldsplit_u";
      else
	prefix = "fieldsplit_p";

      KSPSetType(subksp[m],"preonly");
      KSPGetPC(subksp[m],&subpc[m]);
      if (!strncasecmp(solParams.subprec[m].c_str(),"compositedir",12)) {
	Mat mat;
	Mat Pmat;
#if PETSC_VERSION_MINOR < 5
	MatStructure flag;
	PCGetOperators(subpc[m],&mat,&Pmat,&flag);
#else
	PCGetOperators(subpc[m],&mat,&Pmat);
#endif
	
	if (!solParams.addDirSmoother(subpc[m],Pmat,m,dirIndexSet))
	  return false;
      }
      else
	PCSetType(subpc[m],solParams.subprec[m].c_str());
      
      if (!strncasecmp(solParams.subprec[m].c_str(),"gamg",4) || !strncasecmp(solParams.subprec[m].c_str(),"ml",2) ) {
	PetscInt nloc = coords.size()/solParams.nsd;
	PCSetCoordinates(subpc[m],solParams.nsd,nloc,&coords[0]);
        PCGAMGSetType(subpc[m], PCGAMGAGG); // TODO?
      }
      
      PCFactorSetLevels(subpc[m],solParams.levels[m]);
      if (!strncasecmp(solParams.subprec[m].c_str(),"asm",3)) {
	PCASMSetType(subpc[m],PC_ASM_BASIC);
	PCASMSetOverlap(subpc[m],solParams.overlap[m]);
	if (!locSubdDofsBlock.empty() && !subdDofsBlock.empty()) {
	  const size_t nsubds = subdDofsBlock[m].size();
	  
	  IS isLocSubdDofs[nsubds], isSubdDofs[nsubds];
	  for (size_t i = 0;i < nsubds;i++) {
	    ISCreateGeneral(PETSC_COMM_SELF,locSubdDofsBlock[m][i].size(),&(locSubdDofsBlock[m][i][0]),PETSC_USE_POINTER,&(isLocSubdDofs[i]));
	    ISCreateGeneral(PETSC_COMM_SELF,subdDofsBlock[m][i].size(),&(subdDofsBlock[m][i][0]),PETSC_USE_POINTER,&(isSubdDofs[i]));
	  }
	  PCASMSetLocalSubdomains(subpc[m],nsubds,isSubdDofs,isLocSubdDofs);
	}
	
	if (solParams.asmlu[m]) {
	  KSP* subksp;
	  PC   subsubpc;
	  PetscInt first, nlocal;
	  PCSetUp(subpc[m]);
	  PCASMGetSubKSP(subpc[m],&nlocal,&first,&subksp);
	  
	  for (int i = 0; i < nlocal; i++) {
	    KSPGetPC(subksp[i],&subsubpc);
	    PCSetType(subsubpc,PCLU);
	    KSPSetType(subksp[i],KSPPREONLY);
	  }
	}
      }
      else if (!strncasecmp(solParams.subprec[m].c_str(),"ml",2)) {
	PetscInt n;
        solParams.setMLOptions(prefix, m);

	PCSetFromOptions(subpc[m]);
	PCSetUp(subpc[m]);

	// Settings for coarse solver
	if (!solParams.MLCoarseSolver.empty())
          solParams.setupCoarseSolver(subpc[m], prefix, m);
	
	PCMGGetLevels(subpc[m],&n);
	// Presmoother
	for (int i = 1;i < n;i++) {
	  KSP preksp;
	  PC  prepc;
	  PCMGGetSmoother(subpc[m],i,&preksp);
	  KSPSetType(preksp,"richardson");
	  KSPSetTolerances(preksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,solParams.noPreSmooth[m]);
	  KSPGetPC(preksp,&prepc);

	  // Set smoother
	  std::string smoother;
          if ((i == n-1) && (!solParams.finesmoother.empty())) 	    
            smoother = solParams.finesmoother[m];
          else
	    smoother = solParams.presmoother[m];

	  if (smoother == "asm" || smoother == "asmlu" ) {
	    PCSetType(prepc,"asm");
	    PCASMSetOverlap(prepc,solParams.overlap[m]);
	
	    if (!locSubdDofsBlock.empty() && !subdDofsBlock.empty() && (i==n-1)) {
	      const size_t nsubds = subdDofsBlock[m].size();
	      
	      IS isLocSubdDofs[nsubds], isSubdDofs[nsubds];
	      for (size_t k = 0;k < nsubds;k++) {
		ISCreateGeneral(PETSC_COMM_SELF,locSubdDofsBlock[m][k].size(),&(locSubdDofsBlock[m][k][0]),PETSC_USE_POINTER,&(isLocSubdDofs[k]));
		ISCreateGeneral(PETSC_COMM_SELF,subdDofsBlock[m][k].size(),&(subdDofsBlock[m][k][0]),PETSC_USE_POINTER,&(isSubdDofs[k]));
	      }

	      PCASMSetLocalSubdomains(prepc,nsubds,isSubdDofs,isLocSubdDofs);
	    }

	    PCSetFromOptions(prepc);
	    PCSetUp(prepc);

	    // If LU factorization is used on each subdomain
	    if (smoother == "asmlu") {
	      KSP* subdksp;
	      PC   subdpc;
	      PetscInt first, nlocal;
	      PCASMGetSubKSP(prepc,&nlocal,&first,&subdksp);

	      for (int k = 0; k < nlocal; k++) {
		KSPGetPC(subdksp[k],&subdpc);
		PCSetType(subdpc,PCLU);
		KSPSetType(subdksp[k],KSPPREONLY);
	      }
	    }

	    PCSetFromOptions(prepc);
	    PCSetUp(prepc);
	  }
	  else if (smoother == "compositedir" && (i==n-1)) {
	    if (!solParams.addDirSmoother(prepc,Sp,m,dirIndexSet))
	      return false;
	  } 
	  else 
	    PCSetType(prepc,smoother.c_str());

          PCFactorSetLevels(prepc,solParams.levels[m]);
          KSPSetUp(preksp);
	}
      }
      else if (!strncasecmp(solParams.subprec[m].c_str(),"hypre",5))
        solParams.setHypreOptions(prefix,m);
    }
  }
  else {
    PCSetFromOptions(pc);
    PCSetUp(pc);
  }

  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);  

  PCView(pc,PETSC_VIEWER_STDOUT_WORLD); 

  return true;
}
#endif // HAS_PETSC
