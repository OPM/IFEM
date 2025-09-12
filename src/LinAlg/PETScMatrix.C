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
#include "PETScSchurPC.h"
#include "ProcessAdm.h"
#include "LinAlgInit.h"
#include "SAM.h"

#include <algorithm>
#include <numeric>


namespace {

/*!
  \brief This is a C++ version of the F77 subroutine ADDEM2 (SAM library).
  \details It performs exactly the same tasks, except that \a NRHS always is 1
           and that the system matrix \a SM here is an instance of the
           PETScMatrix class.
*/
void assemPETSc (const Matrix& eM, PETScMatrix& SM, StdVector* SV,
                 const DomainDecomposition& dd,
                 const std::vector<std::array<int,2>>& glb2Blk,
                 const IntVec& meen, const int* meqn,
                 const int* mpmceq, const int* mmceq, const Real* ttcc)
{
  // Get block for equation
  auto getBlk = [&glb2Blk, nBlock = dd.getNoBlocks()](const int ieq, const int jeq)
  {
    if (nBlock > 1)
      return glb2Blk[ieq-1][0] * nBlock + glb2Blk[jeq-1][0];
    else
      return size_t{0};
  };

  // Get equation number in block
  auto getEq = [&glb2Blk, &dd](const int ieq)
  {
    if (dd.getNoBlocks() < 2)
      return dd.getGlobalEq(ieq) - 1;
    else
      return glb2Blk[ieq-1][1] - 1;
  };

  int nedof = meen.size();
  auto A = SM.getBlockMatrices();
  if (A.empty())
    A.push_back(SM.getMatrix());
  for (int j = 1; j <= nedof; ++j) {
    int jeq = meen[j-1];
    if (jeq < 1)
      continue;

    MatSetValue(A[getBlk(jeq, jeq)], getEq(jeq), getEq(jeq), eM(j,j), ADD_VALUES);

    for (int i = 1; i < j; ++i) {
      int ieq = meen[i-1];
      if (ieq < 1)
        continue;

      MatSetValue(A[getBlk(ieq, jeq)], getEq(ieq), getEq(jeq), eM(i,j), ADD_VALUES);
      MatSetValue(A[getBlk(jeq, ieq)], getEq(jeq), getEq(ieq), eM(j,i), ADD_VALUES);
    }
  }

  // Add (appropriately weighted) elements corresponding to constrained
  // (dependent and prescribed) dofs in eM into SM and/or SV
  for (int j = 1; j <= nedof; ++j) {
    int jceq = -meen[j-1];
    if (jceq < 1)
      continue;

    int jp = mpmceq[jceq-1];
    Real c0 = ttcc[jp-1];

    // Add contributions to SV (right-hand-side)
    if (SV)
      for (int i = 1; i <= nedof; ++i) {
        int ieq = meen[i-1];
        int iceq = -ieq;
        if (ieq > 0)
          (*SV)(ieq) -= c0*eM(i,j);
        else if (iceq > 0)
          for (int ip = mpmceq[iceq-1]; ip < mpmceq[iceq]-1; ++ip)
            if (mmceq[ip] > 0) {
              ieq = meqn[mmceq[ip]-1];
              (*SV)(ieq) -= c0*ttcc[ip]*eM(i,j);
            }
      }

    // Add contributions to SM
    for (jp = mpmceq[jceq-1]; jp < mpmceq[jceq]-1; ++jp)
      if (mmceq[jp] > 0) {
        int jeq = meqn[mmceq[jp]-1];
        for (int i = 1; i <= nedof; ++i) {
          int ieq = meen[i-1];
          int iceq = -ieq;
          if (ieq > 0) {
            MatSetValue(A[getBlk(ieq, jeq)], getEq(ieq), getEq(jeq),
                        ttcc[jp]*eM(i,j), ADD_VALUES);
            MatSetValue(A[getBlk(jeq, ieq)], getEq(jeq), getEq(ieq),
                        ttcc[jp]*eM(j,i), ADD_VALUES);
          }
          else if (iceq > 0)
            for (int ip = mpmceq[iceq-1]; ip < mpmceq[iceq]-1; ++ip)
              if (mmceq[ip] > 0) {
                ieq = meqn[mmceq[ip]-1];
                MatSetValue(A[getBlk(ieq, jeq)], getEq(ieq), getEq(jeq),
                            ttcc[ip]*ttcc[jp]*eM(i,j), ADD_VALUES);
              }
        }
      }
  }
}

}


PETScVector::PETScVector (const ProcessAdm& padm) : adm(padm)
{
  VecCreate(*padm.getCommunicator(),&x);
  VecSetFromOptions(x);
  LinAlgInit::increfs();
}


PETScVector::PETScVector (const ProcessAdm& padm, size_t n)
  : StdVector(n), adm(padm)
{
  if (adm.isParallel())
    n = adm.dd.getMaxEq() - adm.dd.getMinEq() + 1;

  VecCreate(*adm.getCommunicator(),&x);
  VecSetSizes(x,n,PETSC_DECIDE);
  VecSetFromOptions(x);
  LinAlgInit::increfs();
}


PETScVector::PETScVector (const ProcessAdm& padm, const Real* values, size_t n)
  : StdVector(values,n), adm(padm)
{
  if (adm.isParallel())
    n = adm.dd.getMaxEq() - adm.dd.getMinEq() + 1;

  VecCreate(*adm.getCommunicator(),&x);
  VecSetSizes(x,n,PETSC_DECIDE);
  VecSetFromOptions(x);
  LinAlgInit::increfs();
}


PETScVector::PETScVector (const PETScVector& vec) :
  StdVector(vec), adm(vec.adm)
{
  VecDuplicate(vec.x,&x);
  VecCopy(vec.x,x);
  LinAlgInit::increfs();
}


PETScVector::~PETScVector ()
{
  VecDestroy(&x);
  LinAlgInit::decrefs();
}


void PETScVector::init (Real value)
{
  StdVector::init(value);
  VecSet(x,value);
}


void PETScVector::redim (size_t n)
{
  VecDestroy(&x);
  VecCreate(*adm.getCommunicator(),&x);
  VecSetSizes(x,adm.dd.getMaxEq()-adm.dd.getMinEq() + 1,PETSC_DECIDE);
  VecSetFromOptions(x);
  StdVector::redim(n);
}


bool PETScVector::endAssembly ()
{
  for (size_t i = 0; i < this->size(); ++i)
    VecSetValue(x, adm.dd.getGlobalEq(i+1)-1, (*this)[i], ADD_VALUES);

  VecAssemblyBegin(x);
  VecAssemblyEnd(x);

  return true;
}


Real PETScVector::L1norm () const
{
  PetscReal val;

  VecNorm(x,NORM_1,&val);
  return val;
}


Real PETScVector::L2norm () const
{
  PetscReal val;

  VecNorm(x,NORM_2,&val);
  return val;
}


Real PETScVector::Linfnorm () const
{
  PetscReal val;

  VecNorm(x,NORM_INFINITY,&val);
  return val;
}


PETScVectors::PETScVectors (const PETScMatrix& A, int size)
  : myA(A)
{
  PetscInt r, c;
  MatGetSize(A.getMatrix(), &r, &c);
  myDim = c;
  vectors.resize(size);
  for (int i = 0; i < size; ++i)
    MatCreateVecs(A.getMatrix(), nullptr, &vectors[i]);
}


PETScVectors::~PETScVectors ()
{
  for (Vec& v : vectors)
    VecDestroy(&v);
}


void PETScVectors::assemble (const Vectors& vecs,
                             const IntVec& meqn, int)
{
  const DomainDecomposition& dd = myA.getDD();
#pragma omp critical
  for (size_t i = 0; i < meqn.size(); ++i)
    for (size_t r = 0; r < vecs.size(); r++)
      VecSetValue(vectors[r], dd.getGlobalEq(meqn[i]+1)-1,
                  vecs[r](1+i), ADD_VALUES);
}


PETScMatrix::PETScMatrix (const ProcessAdm& padm, const LinSolParams& spar)
  : nrow(0), ncol(0), nsp(nullptr), adm(padm), solParams(spar, adm)
{
  // Create matrix object, by default the matrix type is AIJ
  MatCreate(*adm.getCommunicator(),&pA);

  // Create linear solver object
  KSPCreate(*adm.getCommunicator(),&ksp);

  LinAlgInit::increfs();

  if (spar.getNoBlocks() > 1) {
    matvec.resize(spar.getNoBlocks()*spar.getNoBlocks());
    for (Mat& m : matvec) {
      MatCreate(*adm.getCommunicator(), &m);
      MatSetFromOptions(m);
    }
  }

  setParams = true;
  ISsize = 0;
  nLinSolves = 0;
  assembled = false;
  factored = false;
}


PETScMatrix::PETScMatrix (const ProcessAdm& padm, const PETScSolParams& spar)
  : nrow(0), ncol(0), nsp(nullptr), adm(padm), solParams(spar)
{
  // Create linear solver object
  KSPCreate(*adm.getCommunicator(),&ksp);

  LinAlgInit::increfs();

  setParams = true;
  ISsize = 0;
  nLinSolves = 0;
  assembled = false;
  factored = false;
}


PETScMatrix::~PETScMatrix ()
{
  // Deallocation of linear solver object.
  KSPDestroy(&ksp);

  // Deallocation of matrix object.
  MatDestroy(&pA);
  LinAlgInit::decrefs();
  for (Mat& m : matvec)
    MatDestroy(&m);

  for (IS& v : isvec)
    ISDestroy(&v);

  matvec.clear();
}


SystemMatrix* PETScMatrix::copy () const
{
  PETScMatrix* result = new PETScMatrix(this->adm, this->solParams);
  if (m_dd)
    result->m_dd = std::make_unique<DomainDecomposition>(*m_dd);
  MatDuplicate(this->pA, assembled ? MAT_COPY_VALUES : MAT_DO_NOT_COPY_VALUES, &result->pA);
  result->assembled = assembled;
  result->nrow = nrow;
  result->ncol = ncol;
  return result;
}


size_t PETScMatrix::dim (int idim) const
{
  switch (idim) {
  case 1: return nrow;
  case 2: return ncol;
  case 3: return nrow*ncol;
  default: return 0;
  }
}


void PETScMatrix::preAssemble (const std::vector<IntVec>& MMNPC, size_t nel)
{
  int neq = nrow;
  if (m_dd && m_dd->isPartitioned())
    neq = m_dd->getMaxEq() - m_dd->getMinEq() + 1;

  Mat prealloc = preAllocator(neq);

  std::swap(pA, prealloc);
  // Compute the nodal sparsity pattern
  int inod, jnod;
  if (m_dd && m_dd->isPartitioned()) {
    for (int iel : m_dd->getElms())
      for (size_t j = 0; iel > -1 && j < MMNPC[iel].size(); j++)
        if ((jnod = MMNPC[iel][j]+1) > 0)
        {
          const int gjnod = m_dd->getGlobalEq(jnod) - 1;
          MatSetValue(pA, gjnod, gjnod, 0.0, INSERT_VALUES);
          for (size_t i = 0; i < j; i++)
            if ((inod = MMNPC[iel][i]+1) > 0) {
              const int ginod = m_dd->getGlobalEq(inod) - 1;
              MatSetValue(pA, ginod, gjnod, 0.0, INSERT_VALUES);
              MatSetValue(pA, gjnod, ginod, 0.0, INSERT_VALUES);
            }
        }
    this->endAssembly();
  } else {
    for (size_t iel = 0; iel < nel; iel++)
      for (size_t j = 0; j < MMNPC[iel].size(); j++)
        if ((jnod = MMNPC[iel][j]+1) > 0)
        {
          MatSetValue(pA, jnod-1, jnod-1, 0.0 ,INSERT_VALUES);
          for (size_t i = 0; i < j; i++)
            if ((inod = MMNPC[iel][i]+1) > 0) {
              MatSetValue(pA, inod-1, jnod-1, 0.0, INSERT_VALUES);
              MatSetValue(pA, jnod-1, inod-1, 0.0, INSERT_VALUES);
            }
        }
  }
  std::swap(pA, prealloc);

  MatPreallocatorPreallocate(prealloc, PETSC_TRUE, pA);

  MatDestroy(&prealloc);
  MatSetOption(pA, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  MatSetOption(pA, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
}


void PETScMatrix::initAssembly (const SAM& sam, char)
{
  // Get number of local equations in linear system
  PetscInt neq = adm.dd.getMaxEq() - adm.dd.getMinEq() + 1;
  // Set correct number of rows and columns for matrix.
  MatSetSizes(pA,neq,neq,PETSC_DETERMINE,PETSC_DETERMINE);

  // Allocate sparsity pattern
  if (matvec.empty()) {
    MatSetFromOptions(pA);

    // Allocate sparsity pattern
    if (adm.dd.isPartitioned())
      this->setupSparsity(adm.dd.getElms(), sam);
    else {
      std::vector<int> elms(sam.nel);
      std::iota(elms.begin(), elms.end(), 0);
      this->setupSparsity(elms, sam);
    }

    MatSetUp(pA);

    switch (solParams.getLinSysType()) {
    case LinAlg::SPD:
      MatSetOption(pA, MAT_SPD, PETSC_TRUE);
      break;
    case LinAlg::SYMMETRIC:
      MatSetOption(pA, MAT_SYMMETRIC, PETSC_TRUE);
      break;
    default:
      break;
    }

#ifndef SP_DEBUG
    // Do not abort program for allocation error in release mode
    MatSetOption(pA,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
#endif
  } else {
    isvec.resize(adm.dd.getNoBlocks());
    // index sets
    for (size_t i = 0; i < isvec.size(); ++i) {
      IntVec blockEq;
      blockEq.reserve(adm.dd.getMaxEq(i+1) - adm.dd.getMinEq(i+1) + 1);
      for (int leq : adm.dd.getBlockEqs(i)) {
        int eq = adm.dd.getGlobalEq(leq);
        if (eq >= adm.dd.getMinEq() && eq <= adm.dd.getMaxEq())
          blockEq.push_back(eq-1);
      }
      if (adm.dd.isPartitioned())
        std::sort(blockEq.begin(), blockEq.end());

      ISCreateGeneral(*adm.getCommunicator(),blockEq.size(),
                      blockEq.data(),PETSC_COPY_VALUES,&isvec[i]);
    }

    if (adm.dd.isPartitioned())
      this->setupBlockSparsity(adm.dd.getElms(), sam);
    else {
      std::vector<int> elms(sam.nel);
      std::iota(elms.begin(), elms.end(), 0);
      this->setupBlockSparsity(elms, sam);
    }

    MatCreateNest(*adm.getCommunicator(),solParams.getNoBlocks(),isvec.data(),
                  solParams.getNoBlocks(),isvec.data(),matvec.data(),&pA);

 #ifndef SP_DEBUG
    // Do not abort program for allocation error in release mode
    for (Mat& m : matvec)
      MatSetOption(m,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
 #endif
  }

  MatGetSize(pA, &nrow, &ncol);
  assembled = false;
}


bool PETScMatrix::init (int maxEq,
                        const IntMat* elms,
                        const IntMat* neighs,
                        const IntVec* part)
{
  int neq = maxEq;
  m_dd = std::make_unique<DomainDecomposition>();
  if (adm.dd.isPartitioned()) {
    if (part)
      m_dd->setElms(*part, "");
    m_dd->setup(adm,*neighs,*elms);
    neq = m_dd->getMaxEq() - m_dd->getMinEq() + 1;
  } else
    m_dd->setup(maxEq);

  MatSetFromOptions(pA);
  MatSetSizes(pA, neq, neq, PETSC_DETERMINE, PETSC_DETERMINE);
  MatSetUp(pA);
  MatGetSize(pA, &nrow, &ncol);

  return true;
}


Mat PETScMatrix::preAllocator (const int nrows, const int ncols) const
{
  Mat prealloc;
  MatCreate(*adm.getCommunicator(), &prealloc);
  MatSetType(prealloc, MATPREALLOCATOR);
  MatSetSizes(prealloc, nrows, ncols > 0 ? ncols : nrows,
              PETSC_DETERMINE, PETSC_DETERMINE);
  MatSetUp(prealloc);

  return prealloc;
}


std::vector<Mat>
PETScMatrix::preAllocators () const
{
  const size_t blocks = solParams.getNoBlocks();
  std::vector<Mat> prealloc;
  prealloc.resize(blocks*blocks);
  auto itPre = prealloc.begin();
  for (size_t i = 0; i < blocks; ++i)
    for (size_t j = 0; j < blocks; ++j, ++itPre) {
      const int nrows = adm.dd.getMaxEq(i+1) - adm.dd.getMinEq(i+1) + 1;
      const int ncols = adm.dd.getMaxEq(j+1) - adm.dd.getMinEq(j+1) + 1;
      *itPre = preAllocator(nrows, ncols);
    }

  return prealloc;
}


void PETScMatrix::setupSparsity (const IntVec& elms,
                                 const SAM& sam)
{
  PetscInt neq = adm.dd.getMaxEq() - adm.dd.getMinEq() + 1;
  Mat prealloc = preAllocator(neq);

  std::swap(pA, prealloc);
  std::for_each(elms.begin(), elms.end(),
                [this, &sam](const int elm)
                {
                  IntVec meen;
                  sam.getElmEqns(meen, elm+1);
                  this->assemble(Matrix(meen.size(), meen.size()), sam, elm+1);
                });
  this->endAssembly();
  std::swap(pA, prealloc);

  MatPreallocatorPreallocate(prealloc, PETSC_TRUE, pA);

  MatDestroy(&prealloc);
  MatSetOption(pA, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  MatSetOption(pA, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
}


void PETScMatrix::setupBlockSparsity (const IntVec& elms,
                                      const SAM& sam)
{
  auto prealloc = preAllocators();
  this->setupGlb2Blk(sam);

  std::swap(matvec, prealloc);
  std::for_each(elms.begin(), elms.end(),
                [this, &sam](const int elm)
                {
                  IntVec meen;
                  sam.getElmEqns(meen, elm+1);
                  this->assemble(Matrix(meen.size(), meen.size()), sam, elm+1);
                });
  std::swap(matvec, prealloc);

  for (Mat& pmat : prealloc) {
    MatAssemblyBegin(pmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(pmat, MAT_FINAL_ASSEMBLY);
  }

  const size_t blocks = solParams.getNoBlocks();

  auto it = matvec.begin();
  auto itPre = prealloc.begin();
  for (size_t i = 0; i < blocks; ++i)
    for (size_t j = 0; j < blocks; ++j, ++it, ++itPre) {
      const int nrows = adm.dd.getMaxEq(i+1) - adm.dd.getMinEq(i+1) + 1;
      const int ncols = adm.dd.getMaxEq(j+1) - adm.dd.getMinEq(j+1) + 1;
      MatSetSizes(*it, nrows, ncols,
                  PETSC_DETERMINE, PETSC_DETERMINE);
      MatPreallocatorPreallocate(*itPre, PETSC_TRUE, *it);
      MatSetUp(*it);

      MatDestroy(&*itPre);
      MatSetOption(*it, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
      MatSetOption(*it, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    }
}


void PETScMatrix::setupGlb2Blk (const SAM& sam)
{
  // map from SAM indices to block matrix indices
  size_t blocks = solParams.getNoBlocks();
  glb2Blk.resize(sam.neq, {});
  const DomainDecomposition& dd = adm.dd;

  for (int ieq = 1; ieq <= sam.neq; ++ieq)
    for (size_t b = 0; b < blocks; ++b) {
      if (const auto it = dd.getG2LEQ(b+1).find(ieq);
          it != dd.getG2LEQ(b+1).end())
      {
        glb2Blk[ieq-1][0] = b;
        if (adm.isParallel())
          glb2Blk[ieq-1][1] = adm.dd.isPartitioned()
                                  ? adm.dd.getGlobalEq(ieq, b+1)
                                  : adm.dd.getGlobalEq(it->second, b+1);
        else
          glb2Blk[ieq-1][1] = it->second;
        break;
      }
    }
}


bool PETScMatrix::assemble (const Matrix& eM, const SAM& sam, int e)
{
  IntVec meen;
  if (!sam.getElmEqns(meen,e,eM.rows()))
    return false;

#pragma omp critical
  assemPETSc(eM,*this,nullptr,adm.dd,glb2Blk,meen,sam.meqn,sam.mpmceq,sam.mmceq,sam.ttcc);

  return this->flagNonZeroEqs(meen);
}


bool PETScMatrix::assemble (const Matrix& eM, const SAM& sam,
                            SystemVector& B, int e)
{
  StdVector* Bptr = dynamic_cast<StdVector*>(&B);
  if (!Bptr) return false;

  IntVec meen;
  if (!sam.getElmEqns(meen,e,eM.rows()))
    return false;

#pragma omp critical
  assemPETSc(eM,*this,Bptr,adm.dd,glb2Blk,meen,sam.meqn,sam.mpmceq,sam.mmceq,sam.ttcc);

  return this->flagNonZeroEqs(meen);
}


bool PETScMatrix::assemble (const Matrix& eM, const SAM& sam,
                            SystemVector& B, const IntVec& meq)
{
  PETScVector* Bptr = dynamic_cast<PETScVector*>(&B);
  if (!Bptr) return false;

  if (eM.rows() < meq.size() || eM.cols() < meq.size())
    return false;

#pragma omp critical
  assemPETSc(eM,*this,Bptr,adm.dd,glb2Blk,meq,sam.meqn,sam.mpmceq,sam.mmceq,sam.ttcc);

  return this->flagNonZeroEqs(meq);
}


bool PETScMatrix::assemble (const Matrix& eM, const IntVec& meq)
{
#pragma omp critical
  for (size_t i = 0; i < meq.size(); ++i)
    for (size_t j = 0; j < meq.size(); ++j)
      if (m_dd)
        MatSetValue(pA, m_dd->getGlobalEq(meq[i]+1)-1,
                     m_dd->getGlobalEq(meq[j]+1)-1, eM(i+1, j+1), ADD_VALUES);
      else
        MatSetValue(pA, meq[i], meq[j], eM(i+1, j+1), ADD_VALUES);

  return true;
}


bool PETScMatrix::endAssembly ()
{
  MatAssemblyBegin(pA,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pA,MAT_FINAL_ASSEMBLY);
  assembled = true;

  return true;
}


void PETScMatrix::init ()
{
  // Set all matrix elements to zero
  if (matvec.empty())
    MatZeroEntries(pA);
  else for (Mat& m : matvec)
    MatZeroEntries(m);

  assembled = false;
  factored = false;
}


void PETScMatrix::mult (Real alpha)
{
  MatScale(pA, alpha);
}


bool PETScMatrix::add (Real sigma, int ieq)
{
  if (ieq > static_cast<int>(nrow))
    return false;
  else if (ieq > 0) {
    const int eq = adm.dd.getGlobalEq(ieq) - 1;
    MatSetValue(pA, eq, eq, sigma, ADD_VALUES);
  } else
    MatShift(pA, sigma);

  return true;
}


bool PETScMatrix::multiply (const SystemVector& B, SystemVector& C) const
{
  const PETScVector* Bptr = dynamic_cast<const PETScVector*>(&B);
        PETScVector* Cptr = dynamic_cast<PETScVector*>(&C);

  if ((!Bptr) || (!Cptr))
    return false;

  MatMult(pA,Bptr->getVector(),Cptr->getVector());
  return true;
}


bool PETScMatrix::solve (SystemVector& B, Real*)
{
  PETScVector* Bptr = dynamic_cast<PETScVector*>(&B);
  if (!Bptr)
    return false;

  Vec x;
  VecDuplicate(Bptr->getVector(),&x);
  VecCopy(Bptr->getVector(),x);

  bool result = this->solve(x, Bptr->getVector(),
                            solParams.getStringValue("type") != "preonly");
  VecDestroy(&x);

  return result;
}


bool PETScMatrix::solve (const SystemVector& b, SystemVector& x)
{
  const PETScVector* Bptr = dynamic_cast<const PETScVector*>(&b);
  if (!Bptr)
    return false;

  PETScVector* Xptr = dynamic_cast<PETScVector*>(&x);
  if (!Xptr)
    return false;

  return this->solve(Bptr->getVector(),Xptr->getVector(),false);
}


bool PETScMatrix::solve (const Vec& b, Vec& x, bool knoll)
{
  // Reset linear solver
  if (nLinSolves && solParams.hasValue("reset_pc")) {
    const std::string string_val = solParams.getStringValue("reset_pc");
    int val = solParams.getIntValue("reset_pc");
    if (string_val == "all" ||
        (string_val == "first" && nLinSolves == 1) ||
        (val > 0 && nLinSolves % val == 0)) {
      KSPDestroy(&ksp);
      KSPCreate(*adm.getCommunicator(),&ksp);
      setParams = true;
      factored = false;
      adm.cout << "Resetting preconditioner" << std::endl;
    }
  }

#if PETSC_VERSION_MINOR < 5
    KSPSetOperators(ksp,pA,pA, factored ? SAME_PRECONDITIONER : SAME_NONZERO_PATTERN);
#else
    KSPSetOperators(ksp,pA,pA);
    KSPSetReusePreconditioner(ksp, factored ? PETSC_TRUE : PETSC_FALSE);
#endif

  if (setParams) {
    if (!setParameters(true))
      return false;
    setParams = false;
  }

  if (knoll)
    KSPSetInitialGuessKnoll(ksp,PETSC_TRUE);
  else
    KSPSetInitialGuessNonzero(ksp,solParams.getStringValue("type") == "preonly" ?
                                   PETSC_FALSE : PETSC_TRUE);
  KSPSolve(ksp,b,x);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp,&reason);
  if (reason < 0) {
    adm.cout << "\n Linear solve failed with reason " << KSPConvergedReasons[reason] << std::endl;
    return false;
  }

  if (solParams.getIntValue("verbosity") > 1) {
    PetscInt its;
    KSPGetIterationNumber(ksp,&its);
    adm.cout << "\n Iterations for " << solParams.getStringValue("type")
             << " = " << its << std::endl;
  }
  nLinSolves++;
  factored = true;

  return true;
}


bool PETScMatrix::solveMultipleRhs (PETScVectors& B, Matrix& sField)
{
  Vec xg;
  VecScatter ctx;
  if (m_dd && m_dd->isPartitioned())
    VecScatterCreateToAll(B.get(0), &ctx, &xg);

  Vec x;
  VecDuplicate(B.get(0), &x);

  for (size_t i = 0; i < B.size(); ++i) {
    VecAssemblyBegin(B.get(i));
    VecAssemblyEnd(B.get(i));

    if (!this->solve(B.get(i), x, false))
      return false;

    if (m_dd && m_dd->isPartitioned()) {
      VecScatterBegin(ctx, x, xg, INSERT_VALUES, SCATTER_FORWARD);
      VecScatterEnd(ctx, x, xg, INSERT_VALUES, SCATTER_FORWARD);
      PetscScalar* ga;
      VecGetArray(xg, &ga);
      const auto& g2leq = m_dd->getG2LEQ(0);
      for (const auto& it : g2leq)
        sField(i+1, it.second) = ga[it.first-1];
      VecRestoreArray(xg, &ga);
    } else {
      PetscScalar* xa;
      VecGetArray(x, &xa);
      for (size_t eq = 0; eq < B.dim(); ++eq)
        sField(i+1, eq+1) = xa[eq];
      VecRestoreArray(x, &xa);
    }
  }

  VecDestroy(&x);

  if (m_dd && m_dd->isPartitioned()) {
    VecDestroy(&xg);
    VecScatterDestroy(&ctx);
  }

  return true;
}


bool PETScMatrix::solveEig (PETScMatrix& B, RealArray& val,
                            Matrix& vec, int nv, Real shift, int iop)
{
#ifdef HAS_SLEPC
  EPS eps;
  EPSCreate(*adm.getCommunicator(),&eps);

  const auto slepc_mode = std::array{EPS_HEP, EPS_NHEP, EPS_GHEP, EPS_GNHEP};
  EPSSetOperators(eps, pA, iop > 2 ? B.pA : nullptr);
  EPSSetProblemType(eps, slepc_mode[iop-1]);

  EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE);
  EPSSetDimensions(eps, nv, PETSC_DETERMINE, PETSC_DETERMINE);
  EPSSetFromOptions(eps);

  ST st;
  EPSGetST(eps, &st);
  STSetShift(st, shift);

  KSP oldKsp = ksp;
  STGetKSP(st, &ksp);
  this->setParameters(false);
  ksp = oldKsp;

  if (solParams.getIntValue("verbosity") > 0)
    EPSView(eps, PETSC_VIEWER_STDOUT_WORLD);

  EPSSolve(eps);

  PetscInt nconv;
  EPSGetConverged(eps, &nconv);

  PetscInt m, n;
  MatGetSize(pA, &m, &n);
  if (m != n)
    return false;

  Vec xr, xi;
  MatCreateVecs(pA, nullptr, &xr);
  VecDuplicate(xr, &xi);

  val.resize(nv);
  vec.resize(n, nv);

  Vec gr;
  VecScatter ctx;
  if (adm.dd.isPartitioned())
    VecScatterCreateToAll(xr, &ctx, &gr);

  for (int i = 0; i < std::min(nv, nconv); ++i) {
    PetscScalar kr, ki;
    EPSGetEigenpair(eps, i, &kr, &ki, xr, xi);
    val[i] = kr;
    if (adm.dd.isPartitioned()) {
      VecScatterBegin(ctx, xr, gr, INSERT_VALUES, SCATTER_FORWARD);
      VecScatterEnd(ctx, xr, gr, INSERT_VALUES, SCATTER_FORWARD);
      PetscScalar* grarr;
      VecGetArray(gr, &grarr);
      if (adm.dd.isPartitioned())
        for (const auto& it : adm.dd.getG2LEQ(0))
          vec(it.second,i+1) = grarr[it.first-1];
      VecRestoreArray(gr, &grarr);
    } else {
      PetscScalar* xrarr;
      VecGetArray(xr, &xrarr);
      vec.fillColumn(i+1,xrarr);
      VecRestoreArray(xr, &xrarr);
    }
  }

  VecDestroy(&xi);
  VecDestroy(&xr);

  if (adm.dd.isPartitioned()) {
    VecDestroy(&gr);
    VecScatterDestroy(&ctx);
  }

  EPSDestroy(&eps);

  return true;
#else
  return false;
#endif
}


Real PETScMatrix::Linfnorm () const
{
  PetscReal norm;
  MatNorm(pA,NORM_INFINITY,&norm);
  return norm;
}


bool PETScMatrix::setParameters (bool setup)
{
  // Set linear solver method
  KSPSetType(ksp,
             !forcedKSPType.empty() ? forcedKSPType.c_str()
                                    : solParams.getStringValue("type").c_str());
  KSPSetTolerances(ksp,solParams.getDoubleValue("rtol"),
                   solParams.getDoubleValue("atol"),
                   solParams.getDoubleValue("dtol"),
                   solParams.getIntValue("maxits"));
  PC pc;
  KSPGetPC(ksp,&pc);

  if (matvec.empty())
    solParams.setupPC(pc, 0, "", IntSet(), setup);
  else if (matvec.size() > 4) {
    std::cerr << "** PETSCMatrix ** Only two blocks supported for now." << std::endl;
    return false;
  }
  else {
    PCSetType(pc,PCFIELDSPLIT);
    PetscInt nsplit;
    KSP* subksp;
    std::array<PC,2> subpc;

    PCFieldSplitSetIS(pc,"u",isvec[0]);
    PCFieldSplitSetIS(pc,"p",isvec[1]);
    PCFieldSplitSetType(pc,PC_COMPOSITE_SCHUR);
    if (solParams.getStringValue("schur") == "lower")
      PCFieldSplitSetSchurFactType(pc,PC_FIELDSPLIT_SCHUR_FACT_LOWER);
    else if (solParams.getStringValue("schur") == "full")
      PCFieldSplitSetSchurFactType(pc,PC_FIELDSPLIT_SCHUR_FACT_FULL);
    else if (solParams.getStringValue("schur") == "diag")
      PCFieldSplitSetSchurFactType(pc,PC_FIELDSPLIT_SCHUR_FACT_DIAG);
    else
      PCFieldSplitSetSchurFactType(pc,PC_FIELDSPLIT_SCHUR_FACT_UPPER);

    PCFieldSplitSetSchurPre(pc,PC_FIELDSPLIT_SCHUR_PRE_SELFP,nullptr);

    PCSetFromOptions(pc);
    if (setup)
      PCSetUp(pc);
    PCFieldSplitGetSubKSP(pc,&nsplit,&subksp);

    // Preconditioner for blocks
    char pchar = '1';
    for (PetscInt m = 0; m < nsplit; ++m, ++pchar) {
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
      if (solParams.getBlock(m).getStringValue("pc") == "schur")
        new PETScSchurPC(subpc[m], matvec, solParams.getBlock(m), adm);
      else
        solParams.setupPC(subpc[m], m, prefix, adm.dd.getBlockEqs(m), setup);
    }
  }

  KSPSetFromOptions(ksp);
  if (setup)
    KSPSetUp(ksp);

  if (setup && solParams.getIntValue("verbosity") >= 1)
    KSPView(ksp, PETSC_VIEWER_STDOUT_(*adm.getCommunicator()));

  return true;
}


const DomainDecomposition& PETScMatrix::getDD () const
{
  return m_dd ? *m_dd : adm.dd;
}


PETScVector operator* (const SystemMatrix& A, const PETScVector& b)
{
  PETScVector results(b.getAdm());
  A.multiply(b, results);
  return results;
}


PETScVector operator/ (SystemMatrix& A, const PETScVector& b)
{
  PETScVector results(b.getAdm());
  A.solve(b, results);
  return results;
}
