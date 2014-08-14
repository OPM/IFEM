// $Id$
//==============================================================================
//!
//! \file PETScMatrix.C
//!
//! \date Jan 15 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Representation of the system matrix in PETSc format.
//!
//==============================================================================

#include "PETScMatrix.h"
#ifdef HAS_PETSC
#include "LinSolParams.h"
#include "LinAlgInit.h"
#include "SAMpatchPara.h"
#include "ProcessAdm.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif


PETScVector::PETScVector(const ProcessAdm& padm) : adm(padm)
{
  VecCreate(*padm.getCommunicator(),&x);
  VecSetFromOptions(x);
  LinAlgInit::increfs();
}


PETScVector::PETScVector(const ProcessAdm& padm, size_t n) : adm(padm)
{
  VecCreate(*adm.getCommunicator(),&x);
  VecSetSizes(x,n,PETSC_DECIDE);
  VecSetFromOptions(x);
  LinAlgInit::increfs();
}


PETScVector::PETScVector(const ProcessAdm& padm, const Real* values, size_t n) : adm(padm)
{
  VecCreate(*adm.getCommunicator(),&x);
  VecSetSizes(x,n,PETSC_DECIDE);
  VecSetFromOptions(x);
  this->restore(values);
  LinAlgInit::increfs();
}


PETScVector::PETScVector(const PETScVector& vec) : adm(vec.adm)
{
  VecDuplicate(vec.x,&x);
  VecCopy(vec.x,x);
  LinAlgInit::increfs();
}


PETScVector::~PETScVector()
{
  VecDestroy(&x);
  LinAlgInit::decrefs();
}


size_t PETScVector::dim() const
{
  PetscInt size;

  VecGetLocalSize(x,&size);
  return size;
}


void PETScVector::redim(size_t n)
{
  VecDestroy(&x);
  VecCreate(*adm.getCommunicator(),&x);
  VecSetSizes(x,n,PETSC_DECIDE);
  VecSetFromOptions(x);
}


Real* PETScVector::getPtr()
{
  Real* ptr = 0;

  VecGetArray(x,&ptr);
  return ptr;
}


const Real* PETScVector::getRef() const
{
  return const_cast<PETScVector*>(this)->getPtr();
}


void PETScVector::restore(const Real* ptr)
{
  PetscScalar* pptr = (PetscScalar*)ptr;

  VecRestoreArray(x,&pptr);
}


void PETScVector::init(Real value)
{
  VecSet(x,value);
}


bool PETScVector::beginAssembly()
{
  VecAssemblyBegin(x);
  return true;
}


bool PETScVector::endAssembly()
{
  VecAssemblyEnd(x);
  return true;
}


void PETScVector::mult(Real alpha)
{
  VecScale(x,alpha);
}


Real PETScVector::L1norm() const
{
  PetscReal val;

  VecNorm(x,NORM_1,&val);
  return val;
}


Real PETScVector::L2norm() const
{
  PetscReal val;

  VecNorm(x,NORM_2,&val);
  return val;
}


Real PETScVector::Linfnorm() const
{
  PetscReal val;

  VecNorm(x,NORM_INFINITY,&val);
  return val;
}


PETScMatrix::PETScMatrix (const ProcessAdm& padm, const LinSolParams& spar) : adm(padm), solParams(spar)
{
  // Create matrix object, by default the matrix type is AIJ
  MatCreate(*adm.getCommunicator(),&A);

  // Create linear solver object
  KSPCreate(*adm.getCommunicator(),&ksp);

  // Create null space if any
  if (solParams.getNullSpace() == CONSTANT) {
    MatNullSpaceCreate(*adm.getCommunicator(),PETSC_TRUE,0,0,&nsp);
    KSPSetNullSpace(ksp,nsp);
  }

  LinAlgInit::increfs();

  setParams = true;
  elmIS = 0;
  ISsize = 0;
  nIts = 0;
  nLinSolves = 0;
}


PETScMatrix::PETScMatrix (const PETScMatrix& B) : adm(B.adm), solParams(B.solParams)
{
  MatDuplicate(B.A,MAT_COPY_VALUES,&A);

  // Create linear solver object.
  KSPCreate(*adm.getCommunicator(),&ksp);

  // Create null space, if any
  if (solParams.getNullSpace() == CONSTANT) {
    MatNullSpaceCreate(*adm.getCommunicator(),PETSC_TRUE,0,0,&nsp);
    KSPSetNullSpace(ksp,nsp);
  }
  LinAlgInit::increfs();

  setParams = true;
  elmIS = 0;
  ISsize = 0;
}


PETScMatrix::~PETScMatrix ()
{
  // Deallocation of null space
  if (solParams.getNullSpace() == CONSTANT)
    MatNullSpaceDestroy(&nsp);

  // Deallocation of linear solver object.
  KSPDestroy(&ksp);

  // Deallocation of matrix object.
  MatDestroy(&A);
  LinAlgInit::decrefs();

  for (PetscInt i = 0; i < ISsize; i++)
    ISDestroy(&elmIS[i]);

  // Deallocation of index sets
  if (!dirIndexSet.empty()) {
    for (size_t i = 0;i < dirIndexSet.size();i++)
      for (size_t j = 0;j < dirIndexSet[i].size();j++)
	ISDestroy(&dirIndexSet[i][j]);
  }

  delete elmIS;
}


/*!
  \brief This is a C++ version of the F77 subroutine ADDEM2 (SAM library).
  \details It performs exactly the same tasks, except that \a NRHS always is 1.
  and that the system matrix \a SM here is an object of the PETScMatrix class.
  \note This function does not account for multi-point constraint equations.
*/

static bool assemPETSc (const Matrix& eM, Mat SM, PETScVector* SV,
                        const std::vector<int>& meen, const int* meqn,
                        const int* mpmceq, const Real* ttcc)
{
  size_t i, j;

  // Number of degrees of freedom for element
  size_t nedof = meen.size();

  // Convert meen to 0-based C array
  PetscInt* l2g = new PetscInt[nedof];
  for (i = 0; i < nedof; i++)
    l2g[i] = meqn[meen[i]-1] - 1;

  // Make a (transposed) copy of eM to modify for Dirichlet BCs
  Matrix A(eM,true);

  // Get vector of prescribed degrees of freedom for this element
  bool rhsMod = false;
  Vector uc(SV ? nedof : 0), bc;
  for (j = 1; j <= uc.size(); j++) {
    int jceq = mpmceq[meen[j-1]-1];

    if ((jceq > 0) && ((uc(j) = -ttcc[jceq-1]) != Real(0)))
      rhsMod = true;
  }

  // Multiply element matrix with prescribed values to get the RHS-contributions
  if (rhsMod)
    eM.multiply(uc,bc);

  // Eliminate constrained degrees of freedom from element matrix
  for (j = 1; j <= nedof; j++)
    if (mpmceq[meen[j-1]-1] > 0) {
      for (i = 1; i <= nedof; i++)
        A(i,j) = A(j,i) = Real(0);
      A(j,j) = Real(1);

      if (rhsMod)
        bc(j) = -uc(j);
    }

  if (rhsMod) // Add contributions to right-hand side vector (SV)
    VecSetValues(SV->getVector(),nedof,l2g,bc.ptr(),ADD_VALUES);

  // Add element stiffness matrix to global matrix
  MatSetValues(SM,nedof,l2g,nedof,l2g,A.ptr(),ADD_VALUES);

  delete[] l2g;
  return true;
}


void PETScMatrix::initAssembly (const SAM& sam, bool)
{
  const SAMpatchPara* sampch = dynamic_cast<const SAMpatchPara*>(&sam);

  if (!strncasecmp(solParams.getPreconditioner(),"gamg",4) ||
      !strncasecmp(solParams.getPreconditioner(),"ml",2)   ||
      solParams.getNullSpace(1) == RIGID_BODY)
    sampch->getLocalNodeCoordinates(coords);

  if (solParams.getNullSpace(1) == RIGID_BODY) {
#ifdef PARALLEL_PETSC
    std::cerr << "WARNING: Rigid body null space not implemented in parallel, ignoring" << std::endl;
#else
    Vec coordVec;
    VecCreate(PETSC_COMM_SELF, &coordVec);
    VecSetBlockSize(coordVec, sampch->getNoSpaceDim());
    VecSetSizes(coordVec, coords.size(), PETSC_DECIDE);
    VecSetFromOptions(coordVec);
    for (size_t i=0;i<coords.size();++i)
      VecSetValue(coordVec, i, coords[i], INSERT_VALUES);
    MatNullSpaceCreateRigidBody(coordVec, &nsp);
    KSPSetNullSpace(ksp,nsp);
#endif
  }

  if (!solParams.dirOrder.empty()) {
    dirIndexSet.resize(1);
    for (size_t j = 0;j < solParams.dirOrder[0].size();j++) {
      IS permIndex;
      if (solParams.dirOrder[0][j] != 123) {
	PetscIntVec perm;
	sampch->getDirOrdering(perm,solParams.dirOrder[0][j]);
	ISCreateGeneral(*adm.getCommunicator(),perm.size(),&(perm[0]),PETSC_COPY_VALUES,&permIndex);
	ISSetPermutation(permIndex);
      }
      else {
	ISCreate(*adm.getCommunicator(),&permIndex);
	ISSetIdentity(permIndex);
      }
      dirIndexSet[0].push_back(permIndex);
    }
  }

  int nx = solParams.getLocalPartitioning(0);
  int ny = solParams.getLocalPartitioning(1);
  int nz = solParams.getLocalPartitioning(2);
  int overlap = solParams.getOverlap();

  if (nx+ny+nz > 0) {
    sampch->getLocalSubdomains(locSubdDofs,nx,ny,nz);
    sampch->getSubdomains(subdDofs,overlap,nx,ny,nz);
  }

  // Get number of equations in linear system
  const PetscInt neq   = sam.getNoEquations();
  const PetscInt nnod  = sam.getNoNodes();
  const PetscInt bsize = sam.getNoDOFs()/nnod;

  // Set correct number of rows and columns for matrix.
  MatSetSizes(A,neq,neq,PETSC_DECIDE,PETSC_DECIDE);
  MatSetBlockSize(A,bsize);
  MPI_Barrier(*adm.getCommunicator());
  MatSetFromOptions(A);

  // Allocation of sparsity pattern
  if (adm.isParallel()) {
    PetscInt ifirst, ilast;
    std::vector<int> d_nnz, o_nnz;

    // Determine rows owned by this process
    ifirst = sampch->getMinEqNumber();
    ifirst--;
    ilast  = sampch->getMaxEqNumber();

    if (sam.getNoDofCouplings(ifirst,ilast,d_nnz,o_nnz)) {
      size_t i;
      PetscIntVec d_Nnz(d_nnz.size());
      PetscIntVec o_Nnz(o_nnz.size());
      for (i = 0; i < d_nnz.size(); i++)
	d_Nnz[i] = d_nnz[i];
      for (i = 0; i < o_nnz.size(); i++)
	o_Nnz[i] = o_nnz[i];
      MatMPIAIJSetPreallocation(A,PETSC_DEFAULT,&(d_Nnz[0]),PETSC_DEFAULT,&(o_Nnz[0]));
    }
    else {
      const PetscInt maxdofc = sam.getMaxDofCouplings();
      MatMPIAIJSetPreallocation(A,maxdofc,PETSC_NULL,maxdofc,PETSC_NULL);
    }
  }
  else {
    // RUNAR
    std::vector<int> nnz;
    if (sam.getNoDofCouplings(nnz)) {
      PetscIntVec Nnz(nnz.size());
      for (size_t i = 0; i < nnz.size(); i++)
	Nnz[i] = nnz[i];
      MatSeqAIJSetPreallocation(A,PETSC_DEFAULT,&(Nnz[0]));
    }
    else {
      const PetscInt maxdofc = sam.getMaxDofCouplings();
      MatSeqAIJSetPreallocation(A,maxdofc,PETSC_NULL);
    }
#ifdef USE_OPENMP
    // dummy assembly loop to avoid matrix resizes during assembly
    if (omp_get_max_threads() > 1) {
      std::vector<int> irow;
      std::vector<int> jcol;
    PetscIntVec col;
    sam.getDofCouplings(irow,jcol);
    for (size_t i=0;i<jcol.size();++i)
      col.push_back(jcol[i]-1);
    MatSeqAIJSetColumnIndices(A,&col[0]);
    MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
    MatSetOption(A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    }
#endif
  }

#ifndef SP_DEBUG
  // Do not abort program for allocation error in release mode
  MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
#endif
}


bool PETScMatrix::beginAssembly()
{
  // Starts parallel assembly process
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  return true;
}


bool PETScMatrix::endAssembly()
{
  // Finalizes parallel assembly process
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  return true;
}


void PETScMatrix::init ()
{
  // Set all matrix elements to zero
  MatZeroEntries(A);
}


bool PETScMatrix::assemble (const Matrix& eM, const SAM& sam, int e)
{
  // Get mapping "meen" between local degrees of freedom in element e
  // and global degrees of freedom.
  std::vector<int> meen;
  if (!sam.getElmEqns(meen,e,eM.rows()))
    return false;

  // Assemble local stiffness matrix into global system.
  return assemPETSc(eM,A,NULL,meen,sam.meqn,sam.mpmceq,sam.ttcc);
}


bool PETScMatrix::assemble (const Matrix& eM, const SAM& sam,
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
  return assemPETSc(eM,A,Bptr,meen,sam.meqn,sam.mpmceq,sam.ttcc);
}


bool PETScMatrix::multiply (const SystemVector& B, SystemVector& C)
{
  const PETScVector* Bptr = dynamic_cast<const PETScVector*>(&B);
        PETScVector* Cptr = dynamic_cast<PETScVector*>(&C);

  if ((!Bptr) || (!Cptr))
    return false;

  MatMult(A,Bptr->getVector(),Cptr->getVector());
  return true;
}


bool PETScMatrix::solve (SystemVector& B, bool newLHS)
{
  // Reset linear solver
  if (nLinSolves && solParams.nResetSolver)
    if (nLinSolves%solParams.nResetSolver == 0) {
      KSPDestroy(&ksp);
      KSPCreate(*adm.getCommunicator(),&ksp);
      setParams = true;
    }

  const PETScVector* Bptr = dynamic_cast<PETScVector*>(&B);
  if (!Bptr)
    return false;

  Vec x;
  VecDuplicate(Bptr->getVector(),&x);
  VecCopy(Bptr->getVector(),x);

  if (setParams) {
    KSPSetOperators(ksp,A,A, newLHS ? SAME_NONZERO_PATTERN:SAME_PRECONDITIONER);
    solParams.setParams(ksp,locSubdDofs,subdDofs,coords,dirIndexSet);
    setParams = false;
  }
  KSPSetInitialGuessKnoll(ksp,PETSC_TRUE);
  KSPSolve(ksp,x,Bptr->getVector());
  VecDestroy(&x);

  PetscInt its;
  KSPGetIterationNumber(ksp,&its);
  PetscPrintf(PETSC_COMM_WORLD,"\n Iterations for %s = %D\n",
              solParams.getMethod(),its);
  nIts += its;
  nLinSolves++;

  return true;
}


bool PETScMatrix::solve (const SystemVector& b, SystemVector& x, bool newLHS)
{
  // Reset linear solver
  if (nLinSolves && solParams.nResetSolver)
    if (nLinSolves%solParams.nResetSolver == 0) {
      KSPDestroy(&ksp);
      KSPCreate(*adm.getCommunicator(),&ksp);
      setParams = true;
    }

  const PETScVector* Bptr = dynamic_cast<const PETScVector*>(&b);
  if (!Bptr)
    return false;

  PETScVector* Xptr = dynamic_cast<PETScVector*>(&x);
  if (!Xptr)
    return false;

  if (setParams) {
    KSPSetOperators(ksp,A,A, newLHS ? SAME_NONZERO_PATTERN : SAME_PRECONDITIONER);
    solParams.setParams(ksp,locSubdDofs,subdDofs,coords,dirIndexSet);
    setParams = false;
  }
  KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
  KSPSolve(ksp,Bptr->getVector(),Xptr->getVector());

  PetscInt its;
  KSPGetIterationNumber(ksp,&its);
  PetscPrintf(PETSC_COMM_WORLD,"\n Iterations for %s = %D\n",
	      solParams.getMethod(),its);
  nIts += its;
  nLinSolves++;

  return true;
}


bool PETScMatrix::solve (SystemVector& B, SystemMatrix& P, bool newLHS)
{
  // Reset linear solver
  if (nLinSolves && solParams.nResetSolver)
    if (nLinSolves%solParams.nResetSolver == 0) {
      KSPDestroy(&ksp);
      KSPCreate(*adm.getCommunicator(),&ksp);
      setParams = true;
    }

  const PETScVector* Bptr = dynamic_cast<PETScVector*>(&B);
  if (!Bptr)
    return false;

  const PETScMatrix* Pptr = dynamic_cast<PETScMatrix*>(&P);
  if (!Pptr)
    return false;

  Vec x;
  VecDuplicate(Bptr->getVector(),&x);
  VecCopy(Bptr->getVector(),x);

  if (setParams) {
    KSPSetOperators(ksp,A,Pptr->getMatrix(),
                  newLHS ? SAME_NONZERO_PATTERN : SAME_PRECONDITIONER);
    solParams.setParams(ksp,locSubdDofs,subdDofs,coords,dirIndexSet);
    setParams = false;
  }
  KSPSetInitialGuessKnoll(ksp,PETSC_TRUE);
  KSPSolve(ksp,x,Bptr->getVector());
  VecDestroy(&x);

  PetscInt its;
  KSPGetIterationNumber(ksp,&its);
  PetscPrintf(*adm.getCommunicator(),"\n Iterations for %s = %D\n",
		solParams.getMethod(),its);
  nIts += its;
  nLinSolves++;

  return true;
}


bool PETScMatrix::solveEig (PETScMatrix& B, RealArray& val,
			    Matrix& vec, int nv, Real shift, int iop)
{
#ifdef HAS_SLEPC
  ST          st;
  PetscInt    m, n, nconv;
  PetscScalar kr, ki;
  PetscScalar *xrarr;
  Vec         xr, xi;

  EPS eps;
  EPSCreate(*adm.getCommunicator(),&eps);

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

  VecCreate(*adm.getCommunicator(),&xr);
  VecSetSizes(xr,n,PETSC_DETERMINE);
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

  VecDestroy(&xi);
  VecDestroy(&xr);

  EPSDestroy(&eps);

  return true;
#endif
  return false;
}


Real PETScMatrix::Linfnorm () const
{
  PetscReal norm;
  MatNorm(A,NORM_INFINITY,&norm);
  return norm;
}


bool PETScMatrix::makeElementIS(const SAM& sam)
{
  std::vector<int> meen;

  ISsize = sam.getNoElms();

  elmIS = new IS[ISsize];
  for (PetscInt e = 1; e <= ISsize; e++) {
    if (!sam.getElmEqns(meen,e))
      return false;
    std::sort(meen.begin(),meen.end());

    for (size_t i = 0;i < meen.size();i++)
      meen[i]--;

    int ndof = 0;
    for (size_t i = 0;i < meen.size();i++)
      if (meen[i] >= 0)
	ndof++;

    PetscInt l2g[ndof];
    ndof = 0;
    for (size_t i = 0;i < meen.size();i++)
      if (meen[i] >= 0)
	l2g[ndof++] = meen[i];

    ISCreateGeneral(PETSC_COMM_SELF,ndof,l2g,PETSC_COPY_VALUES,&(elmIS[e-1]));
  }

  return true;
}


bool PETScMatrix::makeEBEpreconditioner(const Mat A, Mat* AeI)
{
  PetscInt        nedof;
  const PetscInt* indx;
  Vector          vals;
  PetscIntVec     locidx;

  if (!elmIS)
    return false;

  Mat* subMats;
  MatGetSubMatrices(A,ISsize,elmIS,elmIS,MAT_INITIAL_MATRIX,&subMats);

  MatDuplicate(A,MAT_DO_NOT_COPY_VALUES,AeI);
  MatZeroEntries(*AeI);

  Matrix elmMat;
  for (PetscInt e = 0; e < ISsize; e++) {
    ISGetSize(elmIS[e],&nedof);

    locidx.resize(nedof);
    for (int i = 0;i < nedof;i++)
      locidx[i] = i;

    vals.resize(nedof*nedof,true);
    MatGetValues(subMats[e],nedof,&(locidx.front()),nedof,&(locidx.front()),&(vals.front()));

    int dof = 0;
    elmMat.resize(nedof,nedof);
    for (int i = 1;i <= nedof;i++)
      for (int j = 1;j <= nedof;j++)
	elmMat(i,j) = vals[dof++];

    utl::invert(elmMat);
    elmMat.transpose();

    ISGetIndices(elmIS[e],&indx);
    MatSetValues(*AeI,nedof,indx,nedof,indx,elmMat.ptr(),ADD_VALUES);
    ISRestoreIndices(elmIS[e],&indx);
  }

   MatDestroyMatrices(ISsize,&subMats);

   MatAssemblyBegin(*AeI,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(*AeI,MAT_FINAL_ASSEMBLY);

   return true;
}

#endif // HAS_PETSC
