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
#include "SIMenums.h"
#include "PETScPCPerm.h"
#include "PETScPCScale.h"
#include <sstream>
#ifdef USE_OPENMP
#include <omp.h>
#endif


PETScBlockMatrix::PETScBlockMatrix (const ProcessAdm& padm, const LinSolParams& par) :
  PETScMatrix(padm,par,LinAlg::GENERAL_MATRIX)
{
  ncomps.resize(par.getNoBlocks());
  nblocks = par.getNoBlocks();
  for (size_t i = 0;i < nblocks;i++)
    ncomps[i] = par.getBlock(i).comps/10 + 1;
  
  isvec = (IS*) malloc(sizeof(IS)*nblocks);
  matvec = (Mat*) malloc(sizeof(Mat)*nblocks*nblocks);
  for (size_t m = 0;m < nblocks*nblocks;m++) {
    MatCreate(*adm.getCommunicator(),&matvec[m]);
    MatSetFromOptions(matvec[m]);
  }
  
  S  = nullptr;
  Sp = nullptr;
  Fp = nullptr;

  LinAlgInit::increfs();
}


PETScBlockMatrix::~PETScBlockMatrix ()
{
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

  if (!strncasecmp(solParams.getPreconditioner(),"gamg",4) ||
      !strncasecmp(solParams.getPreconditioner(),"ml",2))
    sampch->getLocalNodeCoordinates(coords);
  else for (size_t i = 0; i < solParams.getNoBlocks(); i++)
    if (!strncasecmp(solParams.getBlock(i).prec.c_str(),"gamg",4) ||
        !strncasecmp(solParams.getBlock(i).prec.c_str(),"ml",2)) {
      sampch->getLocalNodeCoordinates(coords);
      break;
    }

  bool dirsmoothers=false;
  for (size_t i=0; i < solParams.getNoBlocks();++i)
    if (!solParams.getBlock(i).dirSmoother.empty())
      dirsmoothers = true;
  if (dirsmoothers) {
    dirIndexSet.resize(solParams.getNoBlocks());
    for (size_t i = 0;i < solParams.getNoBlocks();i++)
      for (size_t j = 0;j < solParams.getBlock(i).dirSmoother.size();j++) {
	IS permIndex;
	if (solParams.getBlock(i).dirSmoother[j].order != 123) {
	  PetscIntVec perm;
	  sampch->getDirOrdering(perm,solParams.getBlock(i).dirSmoother[j].order,
                                 solParams.getBlock(i).comps);
	  ISCreateGeneral(*adm.getCommunicator(),perm.size(),&(perm[0]),
                          PETSC_COPY_VALUES,&permIndex);
	  ISSetPermutation(permIndex);
	}
	else {
	  ISCreate(*adm.getCommunicator(),&permIndex);
	  ISSetIdentity(permIndex);
	}
	dirIndexSet[i].push_back(permIndex);
      }
  }

  int nx = solParams.getBlock(0).subdomains[0];
  int ny = solParams.getBlock(0).subdomains[1];
  int nz = solParams.getBlock(0).subdomains[2];
  int overlap = solParams.getBlock(1).overlap;

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

  if (solParams.getBlock(1).nullspace == CONSTANT) {
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
     VecNormalize(x, nullptr);
     nsp = new MatNullSpace;
     MatNullSpaceCreate(*adm.getCommunicator(),PETSC_FALSE,1,&x,nsp);
#if PETSC_VERSION_MAJOR > 6
     KSPSetNullSpace(ksp,*nsp);
#else
     MatSetNullSpace(matvec[1],*nsp);
#endif
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


bool PETScBlockMatrix::solve (SystemVector& B, bool newLHS, Real*)
{
  newLHS = nLinSolves==0; // don't reset preconditioner unless requested

  // Reset linear solver
  if (nLinSolves && solParams.getResetSolver())
    if (nLinSolves%solParams.getResetSolver()== 0) {
      if (solParams.getNoBlocks() > 1)
	MatDestroy(&Sp);
      newLHS = true;
    }

  PETScVector* Bptr = dynamic_cast<PETScVector*>(&B);
  if (!Bptr)
    return false;

  Vec x;
  VecDuplicate(Bptr->getVector(),&x);
  this->renumberRHS(Bptr->getVector(),x,true);

  bool result = this->PETScMatrix::solve(x,Bptr->getVector(),newLHS,true);

  // Renumber back to usual numbering
  if (result)
    this->renumberRHS(Bptr->getVector(),Bptr->getVector(),false);

  VecDestroy(&x);

  return result;
}


bool PETScBlockMatrix::solve (const SystemVector& B, SystemVector& X, bool newLHS)
{
  newLHS = nLinSolves==0; // don't reset preconditioner unless requested

  // Reset linear solver
  if (nLinSolves && solParams.getResetSolver())
    if (nLinSolves%solParams.getResetSolver() == 0) {
      if (solParams.getNoBlocks() > 1)
	MatDestroy(&Sp);
      newLHS = true;
    }

  const PETScVector* Bcptr = dynamic_cast<const PETScVector*>(&B);
  PETScVector* Xptr = dynamic_cast<PETScVector*>(&X);
  if (!Bcptr || !Xptr)
    return false;

  // cast off constness
  PETScVector* Bptr = const_cast<PETScVector*>(Bcptr);

  this->renumberRHS(Bptr->getVector(),Bptr->getVector(),true);
  this->renumberRHS(Xptr->getVector(),Xptr->getVector(),true);

  bool result = PETScMatrix::solve(Bptr->getVector(),Xptr->getVector(),newLHS,false);

  // Renumber back to usual numbering
  if (result)
    this->renumberRHS(Xptr->getVector(),Xptr->getVector(),false);

  return result;
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


bool PETScBlockMatrix::setParameters(PETScMatrix *P2, PETScVector *Pb)
{
  PETScBlockMatrix* P = dynamic_cast<PETScBlockMatrix*>(P2);
  if (!P && P2)
    return false;

  // Set linear solver method
  KSPSetType(ksp,solParams.getMethod().c_str());
  KSPSetTolerances(ksp,solParams.getRelTolerance(),
                       solParams.getAbsTolerance(),
                       solParams.getDivTolerance(),
                       solParams.getMaxIterations());

  // Set preconditioner
  PC pc;
  KSPGetPC(ksp,&pc);
  PCSetType(pc,solParams.getPreconditioner());
  if (!strncasecmp(solParams.getPreconditioner(),"gamg",4) ||
      !strncasecmp(solParams.getPreconditioner(),"ml",2)) {
    PetscInt nloc = coords.size()/solParams.getDimension();
    PCSetCoordinates(pc,solParams.getDimension(),nloc,&coords[0]);
  }
  if (!strncasecmp(solParams.getPreconditioner(),"fieldsplit",10)) {
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

    if (!P || (solParams.getSchurType() == SIMPLE))
      MatGetDiagonal(matvec[0],diagA00);
    else
      MatGetDiagonal(P->getMatrixBlock(0,0),diagA00);

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
    KSPSetReusePreconditioner(subksp[1], PETSC_TRUE);
#endif

    if (P && Pb && (solParams.getSchurType() == PCD)) {
      KSPSetType(subksp[1],"preonly");
      KSPGetPC(subksp[1],&subpc[1]);
      PCSetType(subpc[1],PCSHELL);
      PCCreate(*adm.getCommunicator(),&S);
      PCCreate(*adm.getCommunicator(),&Fp);

      const char* prec1 = solParams.getBlock(1).prec.c_str();
      if (strncasecmp(prec1,"compositedir",12))
        PCSetType(S,prec1);
      else if (!solParams.addDirSmoother(S,Sp,1,dirIndexSet))
        return false;

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
      PCSetReusePreconditioner(S, PETSC_TRUE);
#endif
      
      PCProdCreate(&pcprod);
      PCShellSetApply(subpc[1],PCProdApply);
      PCShellSetContext(subpc[1],pcprod);
      PCShellSetDestroy(subpc[1],PCProdDestroy);
      PCShellSetName(subpc[1],"PCD");

      KSPSetType(subksp[1],"preonly");
      if (strncasecmp(prec1,"compositedir",12))
        PCSetType(S,prec1);
      else if (!solParams.addDirSmoother(S,Sp,1,dirIndexSet))
        return false;

      PCFactorSetLevels(S,solParams.getBlock(1).ilu_fill_level);
      if (!strncasecmp(solParams.getBlock(1).prec.c_str(),"asm",3)) {
        solParams.setupAdditiveSchwarz(S, solParams.getBlock(1).overlap,
                                       !strncasecmp(solParams.getBlock(1).prec.c_str(),"asmlu",5),
                                       locSubdDofsBlock.empty()?PetscIntMat():locSubdDofsBlock[1],
                                       subdDofsBlock.empty()?PetscIntMat():subdDofsBlock[1],false);
      }
      else if (!strncasecmp(solParams.getBlock(1).prec.c_str(),"ml",2)) {
	PetscInt n;
        solParams.setMLOptions("fieldsplit_p", 1);

	PCSetFromOptions(S);
	PCSetUp(S);

	// Set coarse solver package
	if (!solParams.getBlock(1).ml.coarsePackage.empty()) {
	  KSP cksp;
	  PC  cpc;
	  PCMGGetCoarseSolve(S,&cksp);
	  KSPGetPC(cksp,&cpc);
	  PCFactorSetMatSolverPackage(cpc,solParams.getBlock(1).ml.coarsePackage.c_str());
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

          if ((i == n-1) && (!solParams.getBlock(1).finesmoother.empty())) {
            smoother = solParams.getBlock(1).finesmoother;
	    noSmooth = solParams.getBlock(1).noFineSmooth;
	  }
          else {
	    smoother = solParams.getBlock(1).presmoother;
	    noSmooth = solParams.getBlock(1).noPreSmooth;
	  }

	  KSPSetTolerances(preksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,noSmooth);
	  KSPGetPC(preksp,&prepc); 

	  if (smoother == "asm" || smoother == "asmlu" ) {
            solParams.setupAdditiveSchwarz(prepc, solParams.getBlock(1).overlap,
                                           smoother == "asmlu",
                                           locSubdDofs, subdDofs, true);
	  }
	  else if (!strncasecmp(smoother.c_str(),"compositeDir",12) && (i==n-1)) {
	    if (!solParams.addDirSmoother(prepc,Sp,1,dirIndexSet))
	      return false;
	  }
	  else
	    PCSetType(prepc,smoother.c_str());
	  
          PCFactorSetLevels(prepc,solParams.getBlock(1).ilu_fill_level);
          KSPSetUp(preksp);
	}
      }
      else if (!strncasecmp(solParams.getBlock(1).prec.c_str(),"gamg",4) ||
               !strncasecmp(solParams.getBlock(1).prec.c_str(),"ml",2) ) {
        PetscInt nloc = coords.size()/solParams.getDimension();
        PCSetCoordinates(pc,solParams.getDimension(),nloc,&coords[0]);
      }

      PCProdSetUp(subpc[1],&QpL,&Fp,&S);
      nsplit--;
    }
      
    // Preconditioner for blocks
    char pchar='1';
    for (PetscInt m = 0; m < nsplit; m++, pchar++) {
      std::string prefix;
      if (nsplit == 2) {
        if (m == 0)
          prefix = "fieldsplit_u";
        else
          prefix = "fieldsplit_p";
      } else
        prefix = std::string("fieldsplit_b")+pchar;

      KSPSetType(subksp[m],"preonly");
      KSPGetPC(subksp[m],&subpc[m]);
      if (!strncasecmp(solParams.getBlock(m).prec.c_str(),"compositedir",12)) {
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
	PCSetType(subpc[m],solParams.getBlock(m).prec.c_str());
      
      if (!strncasecmp(solParams.getBlock(m).prec.c_str(),"gamg",4) ||
          !strncasecmp(solParams.getBlock(m).prec.c_str(),"ml",2) ) {
	PetscInt nloc = coords.size()/solParams.getDimension();
	PCSetCoordinates(subpc[m],solParams.getDimension(),nloc,&coords[0]);
        PCGAMGSetType(subpc[m], PCGAMGAGG); // TODO?
      }
      
      PCFactorSetLevels(subpc[m],solParams.getBlock(m).ilu_fill_level);
      if (!strncasecmp(solParams.getBlock(m).prec.c_str(),"asm",3)) {
        solParams.setupAdditiveSchwarz(subpc[m], solParams.getBlock(m).overlap,
                                       !strncasecmp(solParams.getBlock(m).prec.c_str(),"asmlu",5),
                                       locSubdDofsBlock.empty()?PetscIntMat():locSubdDofsBlock[m],
                                       subdDofsBlock.empty()?PetscIntMat():subdDofsBlock[m],false);
      } else if (!strncasecmp(solParams.getBlock(m).prec.c_str(),"ml",2)) {
        solParams.setMLOptions(prefix, m);

	PCSetFromOptions(subpc[m]);
	PCSetUp(subpc[m]);

	// Settings for coarse solver
	if (!solParams.getBlock(m).ml.coarseSolver.empty())
          solParams.setupCoarseSolver(subpc[m], prefix, m);
	
        solParams.setupSmoothers(subpc[m], m, dirIndexSet,
                                 locSubdDofsBlock.empty()?PetscIntMat():locSubdDofsBlock[m],
                                 subdDofsBlock.empty()?PetscIntMat():subdDofsBlock[m]);
      }
      else if (!strncasecmp(solParams.getBlock(m).prec.c_str(),"hypre",5))
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
