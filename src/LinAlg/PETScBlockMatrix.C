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
#ifdef HAS_PETSC
#include "LinSolParams.h"
#include "SAMpatchPara.h"
#include "petscversion.h"
#include "petscis.h"
#include "petscsys.h"

#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2
#include "petscpcmg.h"
#define PETSCMANGLE(x) &x
#else
#include "petscmg.h"
#define PETSCMANGLE(x) x
#endif

#ifdef HAS_SLEPC
#include "slepceps.h"
#endif

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "LinAlgInit.h"


PETScBlockVector::PETScBlockVector() {}


PETScBlockVector::PETScBlockVector(size_t n, IntVec nc)
{
  nblocks = nc.size();
  ncomps  = nc;

  is = (IS*) malloc(sizeof(IS)*nblocks);
  bvecs = (Vec*) malloc(sizeof(Vec)*nblocks);
  for (size_t m = 0;m < nblocks;m++) {
    VecCreate(PETSC_COMM_WORLD,&bvecs[m]);
    VecSetFromOptions(bvecs[m]);
  }

  LinAlgInit::increfs();
}


PETScBlockVector::PETScBlockVector(const PETScBlockVector& vec)
{
  // VecDuplicate(vec.x,&x);
  // VecCopy(vec.x,x);
  // LinAlgInit::increfs();
}


PETScBlockVector::~PETScBlockVector()
{
  // Deallocate vector blocks
  for (size_t m = 0;m < nblocks; m++) 
    VecDestroy(PETSCMANGLE(bvecs[m]));
  free(bvecs);

  // Deallocate index sets
  for (size_t m = 0;m < nblocks;m++)
    ISDestroy(&is[m]);
  free(is);


  LinAlgInit::decrefs();
}


size_t PETScBlockVector::dim() const
{
  PetscInt size;

  VecGetLocalSize(x,&size);
  return size;
}


size_t PETScBlockVector::dim(size_t n) const
{
  PetscInt size;

  VecGetLocalSize(bvecs[n],&size);
  return size;
}


void PETScBlockVector::redim(size_t n)
{
  size_t nf = 0;
  for (size_t i = 0;i < ncomps.size();i++)
    nf += ncomps[i];
  size_t nnod = n/nf;
  
  for (size_t i = 0;i < nblocks;i++) {
    VecCreate(PETSC_COMM_WORLD,&(bvecs[i]));
    VecSetSizes(bvecs[i],nnod*ncomps[i],PETSC_DECIDE);
    VecSetFromOptions(bvecs[i]);
  }
}


void PETScBlockVector::redim(size_t n, IntVec nc)
{
  ncomps = nc;

  if (nc.size() != nblocks) {
    for (size_t i = 0;i < nblocks;i++)
      VecDestroy(PETSCMANGLE(bvecs[i]));
    free(bvecs);
 
    bvecs = (Vec*) malloc(sizeof(Vec)*nblocks);
    nblocks = nc.size();
  }

  size_t nf = 0;
  for (size_t i = 0;i < ncomps.size();i++)
    nf += ncomps[i];
  size_t nnod = n/nf;
  
  for (size_t i = 0;i < ncomps.size();i++)
    VecSetSizes(bvecs[i],nnod*ncomps[i],PETSC_DECIDE);
}


Real* PETScBlockVector::getPtr()
{
  Real* ptr = 0;

  VecGetArray(x,&ptr);
  return ptr;
}


Real* PETScBlockVector::getPtr(size_t i)
{
  Real* ptr = 0;

  VecGetArray(bvecs[i],&ptr);
  return ptr;
}


const Real* PETScBlockVector::getRef() const
{
  return const_cast<PETScBlockVector*>(this)->getPtr();
}


const Real* PETScBlockVector::getRef(size_t i) const
{
  return const_cast<PETScBlockVector*>(this)->getPtr(i);
}


void PETScBlockVector::restore(const Real* ptr)
{
  PetscScalar* pptr = (PetscScalar*) ptr;
  VecRestoreArray(x,&pptr);
}


void PETScBlockVector::restore(const Real* ptr, size_t i)
{
  PetscScalar* pptr = (PetscScalar*) ptr;
  VecRestoreArray(bvecs[i],&pptr);
}



void PETScBlockVector::init(Real value)
{
  VecSet(x,value);
}

void PETScBlockVector::initAssembly(const SAM& sam)
{
  // Get number of equations in linear system
  const PetscInt neq = sam.getNoEquations();

  // Get number of fields and number of nodes
  size_t nf = 0;
  size_t ncomp = ncomps.size();
  for (size_t i = 0;i < ncomp;i++)
    nf += ncomps[i];
  size_t nnod = neq/nf;

  // Initialize the block vector
  PetscIntVec  idx;

  const int* meqn = sam.getMEQN();

  int stride = 0;
  for (size_t i = 0;i < ncomps.size();i++) {
    PetscInt ncomp  = ncomps[i];
    PetscInt nentry = ncomp*nnod;

    PetscInt lidx = 0; 
    PetscInt gidx = stride;
    idx.resize(nentry);
    for (size_t n = 0;n < nnod;n++) {
      for (int k = 0;k < ncomp;k++)
	idx[lidx++] = meqn[gidx+k]-1; 
      
      gidx += nf;
    }

    stride += ncomps[i];
    ISCreateGeneral(PETSC_COMM_WORLD,nentry,&(idx[0]),PETSC_COPY_VALUES,&(is[i]));
  }

  VecCreateNest(PETSC_COMM_WORLD,nblocks,&(is[0]),&(bvecs[0]),&x);
  VecSetUp(x);
}

bool PETScBlockVector::beginAssembly()
{
  //for (size_t m = 0;m < nblocks;m++)
  //   VecAssemblyBegin(bvecs[m]);
  
  VecAssemblyBegin(x);
  return true;
}


bool PETScBlockVector::endAssembly()
{
  // for (size_t m = 0;m < nblocks;m++)
  //   VecAssemblyEnd(bvecs[m]);

  VecAssemblyEnd(x);
  return true;
}


void PETScBlockVector::mult(Real alpha)
{
  VecScale(x,alpha);
}


Real PETScBlockVector::L1norm() const
{
  PetscReal val;

  VecNorm(x,NORM_1,&val);
  return val;
}


Real PETScBlockVector::L2norm() const
{
  PetscReal val;

  VecNorm(x,NORM_2,&val);
  return val;
}


Real PETScBlockVector::Linfnorm() const
{
  PetscReal val;

  VecNorm(x,NORM_INFINITY,&val);
  return val;
}


PETScBlockMatrix::PETScBlockMatrix(const LinSolParams& spar) : PETScMatrix(spar) 
{
  ncomps.resize(spar.ncomps.size());
  nblocks = ncomps.size();
  for (size_t i = 0;i < nblocks;i++)
    ncomps[i]  = spar.ncomps[i];
  
  isvec = (IS*) malloc(sizeof(IS)*nblocks);
  matvec = (Mat*) malloc(sizeof(Mat)*nblocks*nblocks);
  for (size_t m = 0;m < nblocks*nblocks;m++) {
    MatCreate(PETSC_COMM_WORLD,&matvec[m]);
    MatSetFromOptions(matvec[m]);
  }
  
  LinAlgInit::increfs();
}


PETScBlockMatrix::PETScBlockMatrix(IntVec nc, const LinSolParams& spar) : PETScMatrix(spar)
{
  nblocks = nc.size();
  ncomps.resize(nblocks);
  for (size_t i = 0;i < nblocks;i++)
    ncomps[i]  = nc[i];

  isvec = (IS*) malloc(sizeof(IS)*nblocks);
  matvec = (Mat*) malloc(sizeof(Mat)*nblocks*nblocks);
  for (size_t m = 0;m < nblocks*nblocks;m++) {
    MatCreate(PETSC_COMM_WORLD,&matvec[m]);
    MatSetFromOptions(matvec[m]);
  }

  LinAlgInit::increfs();
}


// PETScBlockMatrix::PETScBlockMatrix (const PETScBlockMatrix& B)
// {
//   MatDuplicate(B.A,MAT_COPY_VALUES,&A);
//
  // // Create linear solver object.
  // KSPCreate(PETSC_COMM_WORLD,&ksp);

  // // Create null space, if any
  // if (solParams.getNullSpace() == CONSTANT) {
  //   MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nsp);
  //   KSPSetNullSpace(ksp,nsp);
  // }
  // LinAlgInit::increfs();

  // elmIS = 0;
  // ISsize = 0;
//}


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

  for (int i = 0;i < ISsize;i++)
    ISDestroy(PETSCMANGLE(elmIS[i]));
  delete elmIS;
}

void PETScBlockMatrix::assemPETScBlock (const Matrix& eM, Mat SM, PETScVector& SV,
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
  
  std::vector<std::vector<Matrix> > eMb;
  std::vector<PetscIntVec> l2gb;
  this->getBlockElmMatData(eM,l2g,eMb,l2gb);

  int m = 0;
  for (i = 0;i < nblocks;i++)
    for (j = 0;j < nblocks;j++, m++) {
      (eMb[i][j]).transpose();
      MatSetValues(matvec[m],l2gb[i].size(),&(l2gb[i][0]),l2gb[j].size(),&(l2gb[j][0]),(eMb[i][j]).ptr(),ADD_VALUES);
    }
}


void PETScBlockMatrix::getBlockElmMatData(const Matrix& Amat, const PetscIntVec& l2g, 
					  std::vector<std::vector<Matrix> >& Ab, 
					  std::vector<PetscIntVec>& l2gb) const
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
//#endif


void PETScBlockMatrix::initAssembly (const SAM& sam, bool)
{
  const SAMpatchPara* sampch = dynamic_cast<const SAMpatchPara*>(&sam);
  
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
#ifdef PARALLEL_PETSC
  int myRank, nProc;
  MPI_Status status;
  MPI_Comm_rank(PETSC_COMM_WORLD,&myRank);
  MPI_Comm_size(PETSC_COMM_WORLD,&nProc);
  if (myRank < nProc-1)
    MPI_Send(&ilast,1,MPI_INT,myRank+1,101,PETSC_COMM_WORLD);
  if (myRank > 0) {
    MPI_Recv(&ifirst,1,MPI_INT,myRank-1,101,PETSC_COMM_WORLD,&status);
  }
  
  if (sampch->getNoDofCouplings(ifirst,ilast,ncomps,d_nnz,o_nnz)) {
    std::vector<PetscInt> d_Nnz;
    std::vector<PetscInt> o_Nnz;
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
#else
  if (sampch->getNoDofCouplings(ifirst,ilast,ncomps,d_nnz,o_nnz)) {
    std::vector<PetscInt> d_Nnz;
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
#endif
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

    ISCreateGeneral(PETSC_COMM_WORLD,nentry,&(idx[0]),PETSC_COPY_VALUES,&(isvec[i]));
  }

  MatCreateNest(PETSC_COMM_WORLD,nblocks,isvec,nblocks,isvec,matvec,&A);

#ifndef SP_DEBUG 
  // Do not abort program for allocation error in release mode
  for (size_t i = 0;i < nblocks*nblocks;i++)
    MatSetOption(matvec[i],MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
#endif

}


bool PETScBlockMatrix::beginAssembly()
{
  // Starts parallel assembly for each submatrix
  //for (size_t m = 0;m < nblocks*nblocks;m++) 
  //  MatAssemblyBegin(matvec[m],MAT_FINAL_ASSEMBLY);
  
  // Starts parallel assembly process
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);

  return true;
}


bool PETScBlockMatrix::endAssembly()
{
  // Finalizes parallel assembly for each submatrix
  //for (size_t m = 0;m < nblocks*nblocks;m++)
  //  MatAssemblyEnd(matvec[m],MAT_FINAL_ASSEMBLY);

  // Finalizes parallel assembly process
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  return true;
}


void PETScBlockMatrix::init ()
{
  // Set all matrix elements to zero for each submatrix
  for (size_t m = 0;m < nblocks*nblocks;m++) 
    MatZeroEntries(matvec[m]);

  // Set all matrix elements to zero
  MatZeroEntries(A);
}


bool PETScBlockMatrix::assemble (const Matrix& eM, const SAM& sam, int e)
{
  return PETScMatrix::assemble(eM,sam,e);
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
  assemPETScBlock(eM,A,*Bptr,meen,sam.meqn,sam.mpmceq,sam.mmceq,sam.ttcc);

  return true;
}


bool PETScBlockMatrix::multiply (const SystemVector& B, SystemVector& C)
{
  const PETScVector* Bptr = dynamic_cast<const PETScVector*>(&B);
        PETScVector* Cptr = dynamic_cast<PETScVector*>(&C);

  if ((!Bptr) || (!Cptr))
    return false;

  MatMult(A,Bptr->getVector(),Cptr->getVector());
  return true;
}


bool PETScBlockMatrix::solve (SystemVector& B, bool newLHS)
{
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
    KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);
  else
    KSPSetOperators(ksp,A,A,SAME_PRECONDITIONER);
  firstIt = false;

  if (setParams) 
    this->setParameters();
  KSPSetInitialGuessKnoll(ksp,PETSC_TRUE);
  KSPSolve(ksp,x,Bptr->getVector());

  // Renumber back to usual numbering
  this->renumberRHS(Bptr->getVector(),Bptr->getVector(),false);

  PetscInt its;
  KSPGetIterationNumber(ksp,&its);
  PetscPrintf(PETSC_COMM_WORLD,"\n Iterations for %s = %D\n",solParams.getMethod(),its);
  VecDestroy(PETSCMANGLE(x));

  setParams = false;

  return true;
}


bool PETScBlockMatrix::solve (const SystemVector& b, SystemVector& x, bool newLHS)
{
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
    KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);
  else
    KSPSetOperators(ksp,A,A,SAME_PRECONDITIONER);

  if (setParams)
    this->setParameters();
  KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
  KSPSolve(ksp,Bptr->getVector(),Xptr->getVector());

  // Renumber back 
  this->renumberRHS(Xptr->getVector(),Xptr->getVector(),false);

  PetscInt its;
  KSPGetIterationNumber(ksp,&its);
  PetscPrintf(PETSC_COMM_WORLD,"\n Iterations for %s = %D\n",solParams.getMethod(),its);

  setParams = false;

  return true;
}


bool PETScBlockMatrix::solve (SystemVector& B, SystemMatrix& P, SystemVector& Pb, bool newLHS)
{
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
    KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);
  //else
  //  KSPSetOperators(ksp,A,A,SAME_PRECONDITIONER);
  firstIt = false;

  if (setParams)
    if (!this->setParameters(Pptr,Pbptr))
      return false;
  KSPSetInitialGuessKnoll(ksp,PETSC_TRUE);
  KSPSolve(ksp,x,Bptr->getVector());

  // Renumber back to usual numbering
  this->renumberRHS(Bptr->getVector(),Bptr->getVector(),false);

  PetscInt its;
  KSPGetIterationNumber(ksp,&its);
  PetscPrintf(PETSC_COMM_WORLD,"\n Iterations for %s = %D\n",solParams.getMethod(),its);
  VecDestroy(PETSCMANGLE(x));

  setParams = false;

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


bool PETScBlockMatrix::solveEig (PETScBlockMatrix& B, RealArray& val,
			    Matrix& vec, int nv, Real shift, int iop)
{
#ifdef HAS_SLEPC
  ST          st;
  PetscInt    m, n, nconv;
  PetscScalar kr, ki;
  PetscScalar *xrarr;
  Vec         xr, xi;

  EPS eps;
  EPSCreate(PETSC_COMM_WORLD,&eps);

  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(B.A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B.A,MAT_FINAL_ASSEMBLY);

  EPSSetOperators(eps,A,B.A);
  EPSSetProblemType(eps,EPS_GHEP);
  EPSSetType(eps,EPSKRYLOVSCHUR);
  EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE);
  EPSGetST(eps,&st);
  STSetShift(st,shift);
  EPSSetDimensions(eps,nv,4*nv,PETSC_NULL);
  EPSSetFromOptions(eps);
  EPSSolve(eps);
  EPSGetConverged(eps,&nconv);

  MatGetSize(A,&m,&n);
  if (m != n) return false;

  VecCreate(PETSC_COMM_WORLD,&xr);
  VecSetSizes(xr,n,PETSC_DECIDE);
  VecSetFromOptions(xr);
  VecDuplicate(xr,&xi);

  val.resize(nv);
  vec.resize(n,nv);
  for (int i = 0; i < nv; i++) {
    EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);
    VecGetArray(xr,&xrarr);
    val[i] = kr;
    vec.fillColumn(i+1,xrarr);
    VecRestoreArray(xr,&xrarr);
  }

  VecDestroy(PETSCMANGLE(xi));
  VecDestroy(PETSCMANGLE(xr));

  EPSDestroy(PETSCMANGLE(eps));

  return true;
#endif
  return false;
}


Real PETScBlockMatrix::Linfnorm () const
{
  PetscReal norm;
  MatNorm(A,NORM_INFINITY,&norm);
  return norm;
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
  //PCFactorSetMatSolverPackage(pc,solParams.package[0].c_str());
  if (!strncasecmp(solParams.prec.c_str(),"fieldsplit",10)) {
    PetscInt m1, m2, n1, n2, nr, nc, nsplit;
    KSP  *subksp;
    PC   subpc[2];
    Vec diagA00;
    MatGetLocalSize(matvec[0],&m1,&n1);
    MatGetLocalSize(matvec[3],&m2,&n2);
#ifdef PARALLEL_PETSC
    VecCreateMPI(PETSC_COMM_WORLD,m1,PETSC_DETERMINE,&diagA00);
#else
    VecCreateSeq(PETSC_COMM_SELF,m1,&diagA00);
#endif

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
    KSPSetOperators(subksp[1],Sp,Sp,SAME_PRECONDITIONER);

    if (P && Pb && (solParams.schurPrec == PCD)) {
      Vec QpL;
      
      KSPSetType(subksp[1],"preonly");
      KSPGetPC(subksp[1],&subpc[1]);
      PCSetType(subpc[1],PCSHELL);
#ifdef PARALLEL_PETSC
      PCCreate(PETSC_COMM_WORLD,&S);
      PCCreate(PETSC_COMM_WORLD,&Fp);
#else
      PCCreate(PETSC_COMM_SELF,&S);
      PCCreate(PETSC_COMM_SELF,&Fp);
#endif
      PCSetType(S,PCLU);
      PCSetType(Fp,PCMAT);
      Vec Pbr;
      VecDuplicate(Pb->getVector(),&Pbr);
      VecCopy(Pb->getVector(),Pbr);
      this->renumberRHS(Pb->getVector(),Pbr,true);
      VecGetSubVector(Pbr,isvec[1],&QpL);
      VecReciprocal(QpL);
      MatDiagonalScale(P->getMatrixBlock(1,1),PETSC_NULL,QpL);
      VecRestoreSubVector(Pbr,isvec[1],&QpL);
      VecDestroy(&Pbr);
      PCSetOperators(S,Sp,Sp,SAME_PRECONDITIONER);
      PCSetOperators(Fp,P->getMatrixBlock(1,1),P->getMatrixBlock(1,1),SAME_NONZERO_PATTERN);
      
      PCProdCreate(&pcprod);
      PCShellSetApply(subpc[1],PCProdApply);
      PCShellSetContext(subpc[1],pcprod);
      PCShellSetDestroy(subpc[1],PCProdDestroy);
      PCShellSetName(subpc[1],"PCD");

      PCProdSetUp(subpc[1],&Fp,&S);
      nsplit--;
    }

    // Preconditioner for blocks
    for (int m = 0;m < nsplit;m++) {
      KSPSetType(subksp[m],"preonly");
      KSPGetPC(subksp[m],&subpc[m]);
      PCSetType(subpc[m],solParams.subprec[m].c_str());
      PCFactorSetLevels(subpc[m],solParams.levels[m]);
      if (!strncasecmp(solParams.subprec[m].c_str(),"asm",3)) {
	PCASMSetType(subpc[m],PC_ASM_BASIC);
	PCASMSetOverlap(subpc[m],solParams.overlap[m]);
	if (!locSubdDofsBlock.empty() && !subdDofsBlock.empty()) {
	  const size_t nsubds = subdDofsBlock[m].size();
	  
	  IS isLocSubdDofs[nsubds], isSubdDofs[nsubds];
	  for (size_t i = 0;i < nsubds;i++) {
	    ISCreateGeneral(PETSC_COMM_WORLD,locSubdDofsBlock[m][i].size(),&(locSubdDofsBlock[m][i][0]),PETSC_USE_POINTER,&(isLocSubdDofs[i]));
	    ISCreateGeneral(PETSC_COMM_WORLD,subdDofsBlock[m][i].size(),&(subdDofsBlock[m][i][0]),PETSC_USE_POINTER,&(isSubdDofs[i]));
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
	
	//PCMGSetNumberSmoothDown(subpc[m],noPreSmooth[m]);
	//PCMGSetNumberSmoothUp(subpc[m],noPostSmooth[m]);
	//PCGAMGSetNlevels(subpc[m],solParams.mglevels[m]);
	//PCGAMGSetCoarseEqLim(subpc[m],solParams.maxCoarseSize[m]);
	std::stringstream maxLevel;
	maxLevel << solParams.mglevels[m];
        PetscOptionsSetValue("-pc_ml_maxNLevels",maxLevel.str().c_str());
	std::stringstream maxCoarseDof;
	maxCoarseDof << solParams.maxCoarseSize[m];
        PetscOptionsSetValue("-pc_ml_maxCoarseSize",maxCoarseDof.str().c_str());
	if (!solParams.MLCoarsenScheme.empty()) {
	  std::stringstream coarsenScheme;
	  coarsenScheme << solParams.MLCoarsenScheme[m];
	  PetscOptionsSetValue("-pc_ml_CoarsenScheme",coarsenScheme.str().c_str());
	}
	if (!solParams.MLThreshold.empty()) {
	  std::stringstream threshold;
	  threshold << solParams.MLThreshold[m];
	  PetscOptionsSetValue("-pc_ml_Threshold",threshold.str().c_str());
	}
	if (!MLDampingFactor.empty()) {
	  std::stringstream damping;
	  damping << MLDampingFactor[m];
	  PetscOptionsSetValue("-pc_ml_DampingFactor",damping.str().c_str());
	}

	PCSetFromOptions(pc);
	PCSetUp(subpc[m]);
	//KSPSetUp(subksp[m]);

	PCMGGetLevels(subpc[m],&n);
	// Presmoother
	for (int i = 1;i < n;i++) {
	  KSP preksp;
	  PC  prepc;
	  // Not working for some reason
	  //PCMGGetSmootherDown(subpc[m],i,&preksp);
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
	    PCASMSetType(prepc,PC_ASM_BASIC);
	    PCASMSetOverlap(prepc,solParams.overlap[m]);
	    
	    if (!locSubdDofs.empty() && !subdDofs.empty()) {
	      const size_t nsubds = subdDofs.size();
	      
	      IS isLocSubdDofs[nsubds], isSubdDofs[nsubds];
	      for (size_t i = 0;i < nsubds;i++) {
		ISCreateGeneral(PETSC_COMM_WORLD,locSubdDofs[i].size(),&(locSubdDofs[i][0]),PETSC_USE_POINTER,&(isLocSubdDofs[i]));
		ISCreateGeneral(PETSC_COMM_WORLD,subdDofs[i].size(),&(subdDofs[i][0]),PETSC_USE_POINTER,&(isSubdDofs[i]));
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
	      
	      for (int i = 0; i < nlocal; i++) {
		KSPGetPC(subksp[i],&subpc);
		PCSetType(subpc,PCLU);
		KSPSetType(subksp[i],KSPPREONLY);
	      }
	    }
	  }
	  else 
	    PCSetType(prepc,smoother.c_str());

          PCFactorSetLevels(prepc,solParams.levels[m]);
          KSPSetUp(preksp);
	}
	
	// Postsmoother
	// for (int i = 1;i < n;i++) {
	//   KSP postksp;
	//   PC  postpc;
	//   // Not working for some reason
	//   //PCMGGetSmootherUp(subpc[m],i,&postksp);
	//   PCMGGetSmoother(subpc[m],i,&postksp);
	//   KSPSetType(postksp,"richardson");
	//   KSPSetTolerances(postksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,solParams.noPostSmooth[m]);
	//   KSPGetPC(postksp,&postpc);

	//   // Set smoother
	//   std::string smoother;

        //   if ((i == n-1) && (!solParams.finesmoother.empty())) 	    
        //     smoother = solParams.finesmoother[m];
        //   else
	//     smoother = solParams.postsmoother[m];

	//   if (smoother == "asm" || smoother == "asmlu" ) {
	//     PCSetType(postpc,"asm");
	//     PCASMSetType(postpc,PC_ASM_BASIC);
	//     PCASMSetOverlap(postpc,solParams.overlap[m]);
	    
	//     if (!locSubdDofs.empty() && !subdDofs.empty()) {
	//       const size_t nsubds = subdDofs.size();
	      
	//       IS isLocSubdDofs[nsubds], isSubdDofs[nsubds];
	//       for (size_t i = 0;i < nsubds;i++) {
	// 	ISCreateGeneral(PETSC_COMM_WORLD,locSubdDofs[i].size(),&(locSubdDofs[i][0]),PETSC_USE_POINTER,&(isLocSubdDofs[i]));
	// 	ISCreateGeneral(PETSC_COMM_WORLD,subdDofs[i].size(),&(subdDofs[i][0]),PETSC_USE_POINTER,&(isSubdDofs[i]));
	//       }
	//       PCASMSetLocalSubdomains(postpc,nsubds,isSubdDofs,isLocSubdDofs);
	//     }

	//     // If LU factorization is used on each subdomain
	//     if (smoother == "asmlu") {
	//       KSP* subksp;
	//       PC   subpc;
	//       PetscInt first, nlocal;
	//       PCSetFromOptions(postpc);
	//       PCSetUp(postpc);
	//       PCASMGetSubKSP(postpc,&nlocal,&first,&subksp);
	      
	//       for (int i = 0; i < nlocal; i++) {
	// 	KSPGetPC(subksp[i],&subpc);
	// 	PCSetType(subpc,PCLU);
	// 	KSPSetType(subksp[i],KSPPREONLY);
	//       }
	//     }
	//   }
	//   else 
	//     PCSetType(postpc,smoother.c_str());

        //   PCFactorSetLevels(postpc,solParams.levels[m]);
        //   KSPSetUp(postksp);
	// }
      }      
    }
  }
  else {
    PCSetFromOptions(pc);
    PCSetUp(pc);
  }

  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);  

  return true;
}


#endif // HAS_PETSC
