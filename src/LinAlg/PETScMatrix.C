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
#include "SAM.h"
#include "petscversion.h"

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


PETScVector::PETScVector()
{
  VecCreate(PETSC_COMM_WORLD,&x);
  VecSetFromOptions(x);
  LinAlgInit::increfs();
}


PETScVector::PETScVector(size_t n)
{
  VecCreate(PETSC_COMM_WORLD,&x);
  VecSetSizes(x,n,PETSC_DECIDE);
  VecSetFromOptions(x);
  LinAlgInit::increfs();
}


PETScVector::PETScVector(const real* values, size_t n)
{
  PetscScalar *x_array;

  VecCreate(PETSC_COMM_WORLD,&x);
  VecSetSizes(x,n,PETSC_DECIDE);
  VecSetFromOptions(x);
  VecGetArray(x,&x_array);
  *x_array = *values;
  VecRestoreArray(x,&x_array);    
  LinAlgInit::increfs();
}


PETScVector::PETScVector(const PETScVector& vec)
{
  VecDuplicate(vec.x,&x);
  VecCopy(vec.x,x);
  LinAlgInit::increfs();
}


PETScVector::~PETScVector()
{
  // Deallocation of vector
  VecDestroy(PETSCMANGLE(x));
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
  VecSetSizes(x,n,PETSC_DECIDE);
}


real* PETScVector::getPtr()
{
  real* ptr = 0;

  VecGetArray(x,&ptr); 
  return ptr;
}


const real* PETScVector::getRef() const
{
  return const_cast<PETScVector*>(this)->getPtr();
}


void PETScVector::restore(const real* ptr)
{
  PetscScalar* pptr = (PetscScalar*) ptr;
  VecRestoreArray(x,&pptr);
}


void PETScVector::init(real value)
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


void PETScVector::mult(real alpha)
{
  VecScale(x,alpha);
}


real PETScVector::L1norm() const
{
  PetscReal val;

  VecNorm(x,NORM_1,&val);
  return val;
}


real PETScVector::L2norm() const
{
  PetscReal val;

  VecNorm(x,NORM_2,&val);
  return val;
}


real PETScVector::Linfnorm() const
{
  PetscReal val;

  VecNorm(x,NORM_INFINITY,&val);
  return val;
}


PETScMatrix::PETScMatrix(const LinSolParams& spar) : solParams(spar)
{
  // Create matrix object, by default the matrix type is AIJ
  MatCreate(PETSC_COMM_WORLD,&A);

  // Create linear solver object
  KSPCreate(PETSC_COMM_WORLD,&ksp);

  // Create null space if any
  if (solParams.getNullSpace() == CONSTANT) {
    MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nsp);
    KSPSetNullSpace(ksp,nsp);
  }
  LinAlgInit::increfs();
}


PETScMatrix::PETScMatrix (const PETScMatrix& B) : solParams(B.solParams)
{
  // Duplicate matrix. 
  //A = MatDuplicate(B.A);
  MatDuplicate(B.A,MAT_COPY_VALUES,&A);

  // Create linear solver object.
  KSPCreate(PETSC_COMM_WORLD,&ksp);

  // Create null space, if any
  if (solParams.getNullSpace() == CONSTANT) {
    MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nsp);
    KSPSetNullSpace(ksp,nsp);
  }
  LinAlgInit::increfs();
}


PETScMatrix::~PETScMatrix ()
{
  // Deallocation of null space
  if (solParams.getNullSpace() == CONSTANT) 
    MatNullSpaceDestroy(PETSCMANGLE(nsp));

  // Deallocation of linear solver object.
  KSPDestroy(PETSCMANGLE(ksp));

  // Deallocation of matrix object.
  MatDestroy(PETSCMANGLE(A));   
  LinAlgInit::decrefs();
}


#ifdef PARALLEL_PETSC
static void assemPETSc (const Matrix& eM, Mat SM, PETScVector& SV,
                        const std::vector<int>& meen, const int* meqn,
                        const int* mpmceq, const int* mmceq,
                        const real* ttcc)
{
  real   c0;
  size_t i, j;
  int    jp, jceq;

  // Number of degrees of freedom for element
  size_t nedof = meen.size();

  // Convert meen to 0-based C array
  PetscInt* l2g = new PetscInt[nedof];
  for (i = 0; i < nedof; i++)
    l2g[i] = meqn[meen[i]-1]-1;

  // Cast to non-constant Matrix to modify for Dirichlet BCs
  Matrix& A = const_cast<Matrix&>(eM);

  std::vector<real> uc, bc;
  uc.resize(nedof);
  bc.resize(nedof);

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
    A.multiply(uc,bc);

  // Eliminate constrained dofs from element matrix
  for (j = 1;j <= nedof;j++) {
    jceq = mpmceq[meen[j-1]-1];
    if (jceq < 1) continue;

    bc[j-1] = -uc[j-1];
    for (i = 1;i <= nedof;i++)
      A(i,j) = A(j,i) = 0.0;
    
    
    A(j,j) = 1.0;
  }

  // Add contributions to SV (righthand side)
  VecSetValues(SV.getVector(),nedof,l2g,&(bc.front()),ADD_VALUES);

  // Add element stiffness matrix to global matrix
  A.transpose();
  MatSetValues(SM,nedof,l2g,nedof,l2g,A.ptr(),ADD_VALUES);
  A.transpose();

  delete[] l2g;
}
#else


static void assemPETSc (const Matrix& eM, Mat SM, PETScVector& SV,
			const std::vector<int>& meen, const int* meqn,
			const int* mpmceq, const int* mmceq, const real* ttcc)
{
  real   c0;
  size_t i, j;
  int    ieq, jeq, ip, jp, iceq, jceq;

  // Get C array
  PetscScalar* svec;
  VecGetArray(SV.getVector(),&svec); 
  
  // Number of degrees of freedom for element
  size_t nedof = meen.size();

  // Convert meen to 0-based C array
  PetscInt* l2g = new PetscInt[nedof];
  for (i = 0; i < nedof; i++)
    l2g[i] = meen[i]-1;

  // Cast to non-constant Matrix 
  Matrix& A = const_cast<Matrix&>(eM);

  // Add element stiffness matrix to global matrix
  A.transpose();
  MatSetValues(SM,nedof,l2g,nedof,l2g,eM.ptr(),ADD_VALUES);
  A.transpose();

  // Add (appropriately weighted) elements corresponding to constrained
  // (dependent and prescribed) dofs in eM into SM and/or SV
  for (j = 1; j <= nedof; j++) {
    jceq = -meen[j-1];
    if (jceq < 1) continue;

    jp = mpmceq[jceq-1];
    c0 = ttcc[jp-1];

    // Add contributions to SV (righthand side)
    if (SV.dim() > 0) 
      for (i = 1; i <= nedof; i++) {
	ieq  = meen[i-1];
	iceq = -ieq;
	if (ieq > 0) 
	  svec[ieq-1] -= c0*eM(i,j);
	else if (iceq > 0)
	  for (ip = mpmceq[iceq-1]; ip < mpmceq[iceq]-1; ip++)
	    if (mmceq[ip] > 0) {
	      ieq = meqn[mmceq[ip]-1];
	      svec[ieq-1] -= c0*ttcc[ip]*eM(i,j);
	    }
      }

    // Add contributions to SM
    for (jp = mpmceq[jceq-1]; jp < mpmceq[jceq]-1; jp++) {
      std::cout << "jp = " << jp << std::endl;
      if (mmceq[jp] > 0) {
	jeq = meqn[mmceq[jp]-1];
	for (i = 1; i <= nedof; i++) {
	  ieq = meen[i-1];
	  iceq = -ieq;
	  if (ieq > 0)
	    if (ieq == jeq) 
	      MatSetValue(SM,ieq-1,jeq-1,(ttcc[jp]+ttcc[jp])*eM(i,j),ADD_VALUES);
	  
	    else {
	      MatSetValue(SM,ieq-1,jeq-1,ttcc[jp]*eM(i,j),ADD_VALUES);
	      MatSetValue(SM,jeq-1,ieq-1,ttcc[jp]*eM(j,i),ADD_VALUES);
	    }
	  else if (iceq > 0)
	    for (ip = mpmceq[iceq-1]; ip < mpmceq[iceq]-1; ip++)
	      if (mmceq[ip] > 0) {
		ieq = meqn[mmceq[ip]-1];
		MatSetValue(SM,ieq-1,jeq-1,ttcc[ip]*ttcc[jp]*eM(i,j),ADD_VALUES);
	      }
	}
      }
    }
    
  }
   
  // Restore vector values
  VecRestoreArray(SV.getVector(),&svec);

  // Deallocate memory for l2g mapping
  delete [] l2g;
}
#endif


static void assemPETSc (const Matrix& eM, Mat SM, const std::vector<int>& meen,
			const int* meqn, const int* mpmceq, const int* mmceq,
			const real* ttcc)
{
  size_t i, j;
  int    ieq, jeq, ip, jp, iceq, jceq;

  // Number of degrees of freedom for element
  size_t nedof = meen.size();

  // Convert meen to 0-based C array
  PetscInt* l2g = new PetscInt[nedof];
  for (i = 0; i < nedof; i++)
    l2g[i] = meen[i]-1;

  // Cast to non-constant Matrix 
  Matrix& A = const_cast<Matrix&>(eM);

  // Add element stiffness matrix to global matrix
  A.transpose();
  MatSetValues(SM,nedof,l2g,nedof,l2g,eM.ptr(),ADD_VALUES);
  A.transpose();

  // Add (appropriately weighted) elements corresponding to constrained
  // (dependent and prescribed) dofs in eM into SM 
  for (j = 1; j <= nedof; j++) {
    jceq = -meen[j-1];
    if (jceq < 1) continue;

    jp = mpmceq[jceq-1];

    // Add contributions to SM
    for (jp = mpmceq[jceq-1]; jp < mpmceq[jceq]-1; jp++) {
      std::cout << "jp = " << jp << std::endl;
      if (mmceq[jp] > 0) {
	jeq = meqn[mmceq[jp]-1];
	for (i = 1; i <= nedof; i++) {
	  ieq = meen[i-1];
	  iceq = -ieq;
	  if (ieq > 0) {
	      MatSetValue(SM,ieq-1,jeq-1,ttcc[jp]*eM(i,j),ADD_VALUES);
	      MatSetValue(SM,jeq-1,ieq-1,ttcc[jp]*eM(j,i),ADD_VALUES);
	  }
	  else if (iceq > 0)
	    for (ip = mpmceq[iceq-1]; ip < mpmceq[iceq]-1; ip++)
	      if (mmceq[ip] > 0) {
		ieq = meqn[mmceq[ip]-1];
		MatSetValue(SM,ieq-1,jeq-1,ttcc[ip]*ttcc[jp]*eM(i,j),ADD_VALUES);
	      }
	}
      }
    }
    
  }

  // Deallocate memory for l2g mapping
  delete [] l2g;
}


void PETScMatrix::initAssembly (const SAM& sam)
{
  // Get number of equations in linear system
  const PetscInt neq = sam.getNoEquations();

  // Set correct number of rows and columns for matrix.
  MatSetSizes(A,neq,neq,PETSC_DECIDE,PETSC_DECIDE);
  MPI_Barrier(PETSC_COMM_WORLD);
  MatSetFromOptions(A);

  // Allocation of sparse pattern.
#ifdef PARALLEL_PETSC
  PetscInt ifirst, ilast;
  std::vector<int> d_nnz, o_nnz;

  MatGetOwnershipRange(A,&ifirst,&ilast);
  if (sam.getNoDofCouplings(ifirst,ilast,d_nnz,o_nnz))
  {
    std::vector<PetscInt> d_Nnz;
    std::vector<PetscInt> o_Nnz;
    d_Nnz.resize(d_nnz.size());
    for (size_t i=0;i<d_nnz.size();++i) 
      d_Nnz[i] = d_nnz[i];
    o_Nnz.resize(o_nnz.size());
    for (size_t i=0;i<o_nnz.size();++i) 
      o_Nnz[i] = o_nnz[i];
    MatMPIAIJSetPreallocation(A,PETSC_DEFAULT,&(d_Nnz[0]),PETSC_DEFAULT,&(o_Nnz[0]));
  }
  else {
    const PetscInt maxdofc = sam.getMaxDofCouplings();
    MatMPIAIJSetPreallocation(A,maxdofc,PETSC_NULL,maxdofc,PETSC_NULL);
  }
#else
  std::vector<int> nnz;

  if (sam.getNoDofCouplings(nnz)) {
    std::vector<PetscInt> Nnz;
    Nnz.resize(nnz.size());
    for (size_t i=0;i<nnz.size();++i) 
      Nnz[i] = nnz[i];
    PetscInt* nnzPtr = &(Nnz[0]);
    MatSeqAIJSetPreallocation(A,PETSC_DEFAULT,nnzPtr);
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
    std::vector<PetscInt> col;
    sam.getDofCouplings(irow,jcol);
    for (size_t i=0;i<jcol.size();++i)
      col.push_back(jcol[i]-1);
    MatSeqAIJSetColumnIndices(A,&col[0]);
    MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
    MatSetOption(A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
  }
#endif
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
  assemPETSc(eM,A,meen,sam.meqn,sam.mpmceq,sam.mmceq,sam.ttcc);

  return true;
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
  assemPETSc(eM,A,*Bptr,meen,sam.meqn,sam.mpmceq,sam.mmceq,sam.ttcc);
  return true;
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
  const PETScVector* Bptr = dynamic_cast<PETScVector*>(&B);
  if (!Bptr)
    return false;

  Vec x;
  VecDuplicate(Bptr->getVector(),&x);
  VecCopy(Bptr->getVector(),x);
  
  // Has lefthand side changed?
  if (newLHS)
    KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);
  else
    KSPSetOperators(ksp,A,A,SAME_PRECONDITIONER);
  
  solParams.setParams(ksp);
  KSPSetInitialGuessKnoll(ksp,PETSC_TRUE);
  KSPSolve(ksp,x,Bptr->getVector());

  PetscInt its;
  KSPGetIterationNumber(ksp,&its);
  PetscPrintf(PETSC_COMM_WORLD,"\n Iterations for %s = %D\n",solParams.getMethod(),its);
  VecDestroy(PETSCMANGLE(x));

  return true;
}


bool PETScMatrix::solve (const SystemVector& b, SystemVector& x, bool newLHS)
{
  SystemVector* Bp = const_cast<SystemVector*>(&b);
  PETScVector* Bptr = dynamic_cast<PETScVector*>(Bp);
  if (!Bptr)
    return false;
  PETScVector* Xptr = dynamic_cast<PETScVector*>(&x);
  if (!Xptr)
    return false;

  // Has lefthand side changed?
  if (newLHS)
    KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);
  else
    KSPSetOperators(ksp,A,A,SAME_PRECONDITIONER);
  
  solParams.setParams(ksp);
  KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
  KSPSolve(ksp,Bptr->getVector(),Xptr->getVector());

  PetscInt its;
  KSPGetIterationNumber(ksp,&its);
  PetscPrintf(PETSC_COMM_WORLD,"\n Iterations for %s = %D\n",solParams.getMethod(),its);

  return true;
}


bool PETScMatrix::solve (SystemVector& B, SystemMatrix& P, bool newLHS)
{
  const PETScVector* Bptr = dynamic_cast<PETScVector*>(&B);
  if (!Bptr)
    return false;

  const PETScMatrix* Pptr = dynamic_cast<PETScMatrix*>(&P);
  if (!Pptr)
    return false;

  Vec x;
  VecDuplicate(Bptr->getVector(),&x);
  VecCopy(Bptr->getVector(),x);

  // Has lefthand side changed?
  if (newLHS)
    KSPSetOperators(ksp,A,Pptr->getMatrix(),SAME_NONZERO_PATTERN);
  else
    KSPSetOperators(ksp,A,Pptr->getMatrix(),SAME_PRECONDITIONER);
  
  solParams.setParams(ksp);
  KSPSetInitialGuessKnoll(ksp,PETSC_TRUE);
  KSPSolve(ksp,x,Bptr->getVector());

  PetscInt its;
  KSPGetIterationNumber(ksp,&its);
  PetscPrintf(PETSC_COMM_WORLD,"\n Iterations for %s = %D\n",solParams.getMethod(),its);
  VecDestroy(PETSCMANGLE(x));

  return true;
}

bool PETScMatrix::solveEig (PETScMatrix& B, RealArray& val,
			    Matrix& vec, int nv, real shift, int iop)
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
  
  VecDestroy(xi);
  VecDestroy(xr);

  EPSDestroy(eps);

  return true;
#endif
  return false;
}


real PETScMatrix::Linfnorm () const
{
  PetscReal norm;
  MatNorm(A,NORM_INFINITY,&norm);
  return norm;
}

#endif // HAS_PETSC
