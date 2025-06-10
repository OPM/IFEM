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
#include <cassert>


PETScVector::PETScVector(const ProcessAdm& padm) : adm(padm)
{
  VecCreate(*padm.getCommunicator(),&x);
  VecSetFromOptions(x);
  LinAlgInit::increfs();
}


PETScVector::PETScVector(const ProcessAdm& padm, size_t n)
  : StdVector(n), adm(padm)
{
  if (adm.isParallel())
    n = adm.dd.getMaxEq() - adm.dd.getMinEq() + 1;

  VecCreate(*adm.getCommunicator(),&x);
  VecSetSizes(x,n,PETSC_DECIDE);
  VecSetFromOptions(x);
  LinAlgInit::increfs();
}


PETScVector::PETScVector(const ProcessAdm& padm, const Real* values, size_t n)
  : StdVector(values,n), adm(padm)
{
  if (adm.isParallel())
    n = adm.dd.getMaxEq() - adm.dd.getMinEq() + 1;

  VecCreate(*adm.getCommunicator(),&x);
  VecSetSizes(x,n,PETSC_DECIDE);
  VecSetFromOptions(x);
  LinAlgInit::increfs();
}


PETScVector::PETScVector(const PETScVector& vec) :
  StdVector(vec), adm(vec.adm)
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


void PETScVector::init(Real value)
{
  StdVector::init(value);
  VecSet(x,value);
}


void PETScVector::redim(size_t n)
{
  VecDestroy(&x);
  VecCreate(*adm.getCommunicator(),&x);
  VecSetSizes(x,adm.dd.getMaxEq()-adm.dd.getMinEq() + 1,PETSC_DECIDE);
  VecSetFromOptions(x);
  StdVector::redim(n);
}


bool PETScVector::endAssembly()
{
  // Poor man's assembleDirect
  if (!adm.isParallel() && adm.dd.getMaxDOF() == 0)
    for (size_t i = 0; i < this->size(); ++i)
      VecSetValue(x, i, (*this)[i], ADD_VALUES);
  else
    for (size_t i = 0; i < this->size(); ++i)
      VecSetValue(x, adm.dd.getGlobalEq(i+1)-1, (*this)[i], ADD_VALUES);

  VecAssemblyBegin(x);
  VecAssemblyEnd(x);

  return true;
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


PETScMatrix::PETScMatrix (const ProcessAdm& padm, const LinSolParams& spar)
  : SparseMatrix(SUPERLU, 1), nsp(nullptr), adm(padm), solParams(spar, adm)
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
}


PETScMatrix::PETScMatrix(const ProcessAdm& padm, const PETScSolParams& spar,
                         const SparseMatrix& A)
    : SparseMatrix(A), nsp(nullptr), adm(padm), solParams(spar)
{
  // Create linear solver object
  KSPCreate(*adm.getCommunicator(),&ksp);

  setParams = true;
  ISsize = 0;
  nLinSolves = 0;
  assembled = false;
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


SystemMatrix* PETScMatrix::copy() const
{
  PETScMatrix* result = new PETScMatrix(this->adm, this->solParams,
                                        static_cast<const SparseMatrix&>(*this));
  if (this->assembled) {
    MatDuplicate(this->pA, MAT_COPY_VALUES, &result->pA);
    result->assembled = true;
  }
  return result;
}


void PETScMatrix::initAssembly (const SAM& sam, char)
{
  this->resize(sam.neq,sam.neq);
  if (!adm.dd.isPartitioned())
    this->preAssemble(sam,false);

  // Get number of local equations in linear system
  PetscInt neq = adm.dd.getMaxEq() - adm.dd.getMinEq() + 1;
  // Set correct number of rows and columns for matrix.
  MatSetSizes(pA,neq,neq,PETSC_DETERMINE,PETSC_DETERMINE);

  // Allocate sparsity pattern
  std::vector<IntSet> dofc;
  if (!adm.dd.isPartitioned())
    sam.getDofCouplings(dofc);

  if (matvec.empty()) {
    MatSetFromOptions(pA);

    // Allocate sparsity pattern
    if (adm.isParallel() && !adm.dd.isPartitioned())
      this->setupSparsityDD(sam);
    else if (adm.dd.isPartitioned())
      this->setupSparsityPartitioned(sam);
    else
      this->setupSparsitySerial(sam);

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
      blockEq.reserve(adm.dd.getMaxEq(i+1)-adm.dd.getMinEq(i+1)+1);
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

    if (adm.isParallel() && !adm.dd.isPartitioned())
      this->setupBlockSparsityDD(sam);
    else if (adm.dd.isPartitioned())
      this->setupBlockSparsityPartitioned(sam);
    else
      this->setupBlockSparsitySerial(sam);

    MatCreateNest(*adm.getCommunicator(),solParams.getNoBlocks(),isvec.data(),
                  solParams.getNoBlocks(),isvec.data(),matvec.data(),&pA);

 #ifndef SP_DEBUG
    // Do not abort program for allocation error in release mode
    for (Mat& m : matvec)
      MatSetOption(m,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
 #endif
  }

  assembled = false;
}


void PETScMatrix::setupSparsityDD (const SAM& sam)
{
  // Allocate sparsity pattern
  std::vector<IntSet> dofc;
  sam.getDofCouplings(dofc);

  int ifirst = adm.dd.getMinEq();
  int ilast  = adm.dd.getMaxEq();
  PetscInt neq = ilast - ifirst + 1;
  PetscIntVec d_nnz(neq, 0);
  IntVec o_nnz_g(adm.dd.getNoGlbEqs(), 0);
  for (int i = 0; i < sam.neq; ++i) {
    int eq = adm.dd.getGlobalEq(i+1);
    if (eq >= adm.dd.getMinEq() && eq <= adm.dd.getMaxEq())
      for (int leq : dofc[i]) {
        int g = adm.dd.getGlobalEq(leq);
        if (g >= adm.dd.getMinEq() && g <= adm.dd.getMaxEq())
          ++d_nnz[eq-adm.dd.getMinEq()];
        else if (g > 0)
          ++o_nnz_g[eq-1];
      }
    else
      o_nnz_g[eq-1] += dofc[i].size();
  }

  adm.allReduceAsSum(o_nnz_g);

  PetscIntVec o_nnz(o_nnz_g.begin()+ifirst-1, o_nnz_g.begin()+ilast);

  // TODO: multiplier cause big overallocation due to no multiplicity handling
  for (PetscInt& nnz : o_nnz)
    nnz = std::min(nnz, adm.dd.getNoGlbEqs());

  MatMPIAIJSetPreallocation(pA,PETSC_DEFAULT,d_nnz.data(),
                            PETSC_DEFAULT,o_nnz.data());
}


void PETScMatrix::setupSparsityPartitioned (const SAM& sam)
{
  // Setup sparsity pattern for global matrix
  PetscInt neq = adm.dd.getMaxEq() - adm.dd.getMinEq() + 1;
  SparseMatrix* lA = new SparseMatrix(neq, sam.neq);
  int iMin = adm.dd.getMinEq(0);
  int iMax = adm.dd.getMaxEq(0);
  for (int elm = 1; elm <= sam.nel; ++elm) {
    IntSet meen;
    sam.getUniqueEqns(meen,elm);
    for (int i : meen)
      if (adm.dd.getGlobalEq(i) >= iMin && adm.dd.getGlobalEq(i) <= iMax)
        for (int j : meen)
          (*lA)(adm.dd.getGlobalEq(i)-iMin+1,adm.dd.getGlobalEq(j)) = 0.0;
  }
  IntVec iA, jA;
  lA->calcCSR(iA,jA);
  delete lA;
  MatMPIAIJSetPreallocationCSR(pA, iA.data(), jA.data(), nullptr);

  // Setup sparsity pattern for local matrix
  for (int elm : adm.dd.getElms()) {
    IntSet meen;
    sam.getUniqueEqns(meen,elm+1);
    for (int i : meen)
      for (int j : meen)
        (*this)(i,j) = 0.0;
  }
  this->optimiseCols();
}


void PETScMatrix::setupSparsitySerial (const SAM& sam)
{
  std::vector<IntSet> dofc;
  sam.getDofCouplings(dofc);
  PetscIntVec Nnz;
  for (const IntSet& dofs : dofc)
    Nnz.push_back(dofs.size());

  MatSeqAIJSetPreallocation(pA,PETSC_DEFAULT,Nnz.data());

  PetscIntVec col;
  for (const IntSet& dofs : dofc)
    for (int dof : dofs)
      col.push_back(dof-1);

  MatSeqAIJSetColumnIndices(pA,&col[0]);
  MatSetOption(pA, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  MatSetOption(pA, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
}


std::vector<std::array<int,2>> PETScMatrix::setupGlb2Blk (const SAM& sam)
{
  // map from sparse matrix indices to block matrix indices
  glb2Blk.resize(A.size());
  size_t blocks = solParams.getNoBlocks();
  const DomainDecomposition& dd = adm.dd;
  std::vector<std::array<int,2>> eq2b(sam.neq, {{-1, 0}}); // cache

  for (size_t j = 0; j < cols(); ++j)
    for (int i = IA[j]; i < IA[j+1]; ++i) {
      int iblk = -1;
      int jblk = -1;
      if (eq2b[JA[i]][0] != -1) {
        iblk = eq2b[JA[i]][0];
        glb2Blk[i][1] = eq2b[JA[i]][1];
      } if (eq2b[j][0] != -1) {
        jblk = eq2b[j][0];
        glb2Blk[i][2] = eq2b[j][1];
      }

      for (size_t b = 0; b < blocks && (iblk == -1 || jblk == -1); ++b) {
        std::map<int,int>::const_iterator it;
        if (iblk == -1 && (it = dd.getG2LEQ(b+1).find(JA[i]+1)) != dd.getG2LEQ(b+1).end()) {
          iblk = b;
          eq2b[JA[i]][0] = b;
          eq2b[JA[i]][1] = glb2Blk[i][1] = it->second-1;
        }

        if (jblk == -1 && (it = dd.getG2LEQ(b+1).find(j+1)) != dd.getG2LEQ(b+1).end()) {
          jblk = b;
          eq2b[j][0] = b;
          eq2b[j][1] = glb2Blk[i][2] = it->second-1;
        }
      }
      if (iblk == -1 || jblk == -1) {
        std::cerr << "Failed to map (" << JA[i]+1 << ", " << j+1 << ") " << std::endl;
        std::cerr << "iblk: " << iblk << ", jblk: " << jblk << std::endl;
        assert(0);
      }
      glb2Blk[i][0] = iblk*solParams.getNoBlocks() + jblk;
    }

  return eq2b;
}


void PETScMatrix::setupGlb2BlkPart (const SAM& sam)
{
  // map from sparse matrix indices to block matrix indices
  glb2Blk.resize(A.size());
  size_t blocks = solParams.getNoBlocks();
  const DomainDecomposition& dd = adm.dd;
  std::vector<std::array<int,2>> eq2b(sam.neq, {{-1, 0}}); // cache

  for (size_t j = 0; j < cols(); ++j)
    for (int i = IA[j]; i < IA[j+1]; ++i) {
      int iblk = -1;
      int jblk = -1;
      if (eq2b[JA[i]][0] != -1) {
        iblk = eq2b[JA[i]][0];
        glb2Blk[i][1] = eq2b[JA[i]][1];
      } if (eq2b[j][0] != -1) {
        jblk = eq2b[j][0];
        glb2Blk[i][2] = eq2b[j][1];
      }

      for (size_t b = 0; b < blocks && (iblk == -1 || jblk == -1); ++b) {
        if (iblk == -1 && dd.getBlockEqs(b).find(JA[i]+1) != dd.getBlockEqs(b).end()) {
          iblk = b;
          eq2b[JA[i]][0] = b;
          eq2b[JA[i]][1] = glb2Blk[i][1] = dd.getMLGEQ(b+1)[JA[i]]-1;
        }

        if (jblk == -1 && dd.getBlockEqs(b).find(j+1) != dd.getBlockEqs(b).end()) {
          jblk = b;
          eq2b[j][0] = b;
          eq2b[j][1] = glb2Blk[i][2] = dd.getMLGEQ(b+1)[j]-1;
        }
      }
      if (iblk == -1 || jblk == -1) {
        std::cerr << "Failed to map (" << JA[i]+1 << ", " << j+1 << ") " << std::endl;
        std::cerr << "iblk: " << iblk << ", jblk: " << jblk << std::endl;
        assert(0);
      }
      glb2Blk[i][0] = iblk*solParams.getNoBlocks() + jblk;
    }
}


void PETScMatrix::setupBlockSparsityDD (const SAM& sam)
{
  size_t blocks = solParams.getNoBlocks();
  const DomainDecomposition& dd = adm.dd;
  std::vector<IntSet> dofc;
  sam.getDofCouplings(dofc);
  std::vector<std::array<int,2>> eq2b = this->setupGlb2Blk(sam);

  std::vector<PetscIntVec> d_nnz(blocks*blocks);
  std::vector<IntVec> o_nnz_g(blocks*blocks);
  size_t k = 0;
  for (size_t i = 0; i < blocks; ++i)
    for (size_t j = 0; j < blocks; ++j, ++k) {
      d_nnz[k].resize(dd.getMaxEq(i+1)-dd.getMinEq(i+1)+1);
      o_nnz_g[k].resize(dd.getNoGlbEqs(i+1));
    }

  for (int i = 0; i < sam.neq; ++i) {
    int blk = eq2b[i][0]+1;
    int row = eq2b[i][1]+1;
    int grow = dd.getGlobalEq(row, blk);

    if (grow >= dd.getMinEq(blk) && grow <= dd.getMaxEq(blk))
      for (int dof : dofc[i]) {
        int cblk = eq2b[dof-1][0]+1;
        int col = eq2b[dof-1][1]+1;
        int gcol = dd.getGlobalEq(col, cblk);
        if (gcol >= dd.getMinEq(cblk) && gcol <= dd.getMaxEq(cblk))
          ++d_nnz[(blk-1)*blocks + cblk-1][grow-dd.getMinEq(blk)];
        else
          ++o_nnz_g[(blk-1)*blocks + cblk-1][grow-1];
      }
    else
      for (int dof : dofc[i]) {
        int cblk = eq2b[dof-1][0]+1;
        ++o_nnz_g[(blk-1)*blocks+cblk-1][grow-1];
      }
  }

  for (size_t i = k = 0; i < blocks; ++i)
    for (size_t j = 0; j < blocks; ++j, ++k) {
      int nrows = dd.getMaxEq(i+1)-dd.getMinEq(i+1)+1;
      int ncols = dd.getMaxEq(j+1)-dd.getMinEq(j+1)+1;
      MatSetSizes(matvec[k], nrows, ncols,
                  PETSC_DETERMINE, PETSC_DETERMINE);
      MatSetFromOptions(matvec[k]);
      adm.allReduceAsSum(o_nnz_g[k]);
      PetscIntVec o_nnz(o_nnz_g[k].begin()+dd.getMinEq(i+1)-1, o_nnz_g[k].begin()+dd.getMaxEq(i+1));

      // TODO: multiplier cause big overallocation due to no multiplicity handling
      for (PetscInt& nnz : o_nnz)
        nnz = std::min(nnz, dd.getNoGlbEqs(j+1));

      MatMPIAIJSetPreallocation(matvec[k],PETSC_DEFAULT,d_nnz[k].data(),
                                PETSC_DEFAULT,o_nnz.data());
      MatSetUp(matvec[k]);
    }
}


void PETScMatrix::setupBlockSparsityPartitioned (const SAM& sam)
{
  size_t blocks = solParams.getNoBlocks();
  const DomainDecomposition& dd = adm.dd;

  // Setup sparsity pattern for local matrix
  for (int elm : adm.dd.getElms()) {
    IntSet meen;
    sam.getUniqueEqns(meen,elm+1);
    for (int i : meen)
      for (int j : meen)
        (*this)(i,j) = 0.0;
  }
  this->optimiseCols();

  this->setupGlb2BlkPart(sam);

  std::vector<Mat> prealloc;
  prealloc.resize(blocks*blocks);
  auto itPre = prealloc.begin();
  for (size_t i = 0; i < blocks; ++i)
    for (size_t j = 0; j < blocks; ++j, ++itPre) {
      MatCreate(*adm.getCommunicator(), &(*itPre));
      MatSetType(*itPre, MATPREALLOCATOR);
      int nrows = dd.getMaxEq(i+1)-dd.getMinEq(i+1)+1;
      int ncols = dd.getMaxEq(j+1)-dd.getMinEq(j+1)+1;
      MatSetSizes(*itPre, nrows, ncols,
                  PETSC_DETERMINE, PETSC_DETERMINE);
      MatSetUp(*itPre);
    }

  Mat pBlock;
  MatCreate(*adm.getCommunicator(), &pBlock);
  MatCreateNest(*adm.getCommunicator(),solParams.getNoBlocks(),isvec.data(),
                solParams.getNoBlocks(),isvec.data(),prealloc.data(),&pBlock);

  std::swap(matvec, prealloc);
  std::swap(pBlock, pA);
  this->endAssembly();
  std::swap(pBlock, pA);
  std::swap(matvec, prealloc);
  MatDestroy(&pBlock);

  auto it = matvec.begin();
  itPre = prealloc.begin();
  for (size_t i = 0; i < blocks; ++i)
    for (size_t j = 0; j < blocks; ++j, ++it, ++itPre) {
      int nrows = dd.getMaxEq(i+1)-dd.getMinEq(i+1)+1;
      int ncols = dd.getMaxEq(j+1)-dd.getMinEq(j+1)+1;
      MatSetSizes(*it, nrows, ncols,
                  PETSC_DETERMINE, PETSC_DETERMINE);

      MatPreallocatorPreallocate(*itPre, PETSC_TRUE, *it);

      MatSetUp(*it);
    }

  for (Mat& pmat : prealloc)
    MatDestroy(&pmat);
}


void PETScMatrix::setupBlockSparsitySerial (const SAM& sam)
{
  size_t blocks = solParams.getNoBlocks();
  const DomainDecomposition& dd = adm.dd;
  std::vector<IntSet> dofc;
  sam.getDofCouplings(dofc);
  this->setupGlb2Blk(sam);

  auto it = matvec.begin();
  for (size_t i = 0; i < blocks; ++i)
    for (size_t j = 0; j < blocks; ++j, ++it) {
      std::vector<PetscInt> nnz;
      nnz.reserve(dd.getBlockEqs(i).size());
      for (int leq : dd.getBlockEqs(i))
        nnz.push_back(std::min(dofc[leq-1].size(), dd.getBlockEqs(j).size()));

      int nrows = dd.getMaxEq(i+1)-dd.getMinEq(i+1)+1;
      int ncols = dd.getMaxEq(j+1)-dd.getMinEq(j+1)+1;
      MatSetSizes(*it, nrows, ncols,
                  PETSC_DETERMINE, PETSC_DETERMINE);

      MatSeqAIJSetPreallocation(*it, PETSC_DEFAULT, nnz.data());
      MatSetUp(*it);
    }
}


bool PETScMatrix::endAssembly ()
{
  if (!this->SparseMatrix::endAssembly())
    return false;

  if (IA.empty() && !assembled) {
    this->optimiseCols();
    return this->assembleDirect();
  }

  for (size_t j = 0; j < cols(); ++j)
    for (int i = IA[j]; i < IA[j+1]; ++i)
      if (matvec.empty())
        MatSetValue(pA,
                    adm.dd.getGlobalEq(JA[i]+1)-1,
                    adm.dd.getGlobalEq(j+1)-1,
                    A[i], ADD_VALUES);
      else if (adm.dd.isPartitioned())
        MatSetValue(matvec[glb2Blk[i][0]],
                    glb2Blk[i][1],
                    glb2Blk[i][2],
                    A[i], ADD_VALUES);
      else
      {
        int rblock = glb2Blk[i][0] / adm.dd.getNoBlocks() + 1;
        int cblock = glb2Blk[i][0] % adm.dd.getNoBlocks() + 1;
        MatSetValue(matvec[glb2Blk[i][0]],
                    adm.dd.getGlobalEq(glb2Blk[i][1]+1, rblock)-1,
                    adm.dd.getGlobalEq(glb2Blk[i][2]+1, cblock)-1,
                    A[i], ADD_VALUES);
      }

  MatAssemblyBegin(pA,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pA,MAT_FINAL_ASSEMBLY);

  assembled = true;

  return true;
}


void PETScMatrix::init ()
{
  this->SparseMatrix::init();

  // Set all matrix elements to zero
  if (matvec.empty())
    MatZeroEntries(pA);
  else for (Mat& m : matvec)
    MatZeroEntries(m);

  assembled = false;
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

  if (A.empty() || !assembled)
    return this->solveDirect(*Bptr);

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

  if (setParams) {
#if PETSC_VERSION_MINOR < 5
    KSPSetOperators(ksp,pA,pA, factored ? SAME_PRECONDITIONER : SAME_NONZERO_PATTERN);
#else
    KSPSetOperators(ksp,pA,pA);
    KSPSetReusePreconditioner(ksp, factored ? PETSC_TRUE : PETSC_FALSE);
#endif
    if (!setParameters())
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


bool PETScMatrix::assembleDirect()
{
  MatSetSizes(pA, this->dim(1), this->dim(2),
              PETSC_DETERMINE, PETSC_DETERMINE);

  MatSeqAIJSetPreallocation(pA, 10, nullptr);
  MatSetOption(pA, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE);
  MatSetUp(pA);

  for (size_t j = 0; j < cols(); ++j)
    for (int i = IA[j]; i < IA[j+1]; ++i)
      MatSetValue(pA, JA[i], j, A[i], INSERT_VALUES);

  MatAssemblyBegin(pA,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pA,MAT_FINAL_ASSEMBLY);
  this->assembled = true;

  return true;
}


bool PETScMatrix::solveDirect(PETScVector& B)
{
  // the sparsity pattern has been grown in-place, we need to init PETsc state.
  // this is currently only used for patch-global L2 systems.
  if (A.empty() && !this->optimiseCols())
    return false;

  // Set correct number of rows and columns for matrix.
  size_t nrow = IA.size()-1;
  if (nrow == 0 || IA.empty())
    return false;

  MatSetSizes(pA, nrow, nrow, PETSC_DECIDE, PETSC_DECIDE);
  MatSetFromOptions(pA);
  PetscInt max = 0;
  for (size_t i = 0; i < nrow; ++i) // symmetric so row/column sizes should be the same
    if (IA[i+1]-IA[i] > max)
      max = IA[i+1]-IA[i];
  MatSeqAIJSetPreallocation(pA, max, nullptr);
  MatSetOption(pA, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE);
  MatSetUp(pA);

  for (size_t j = 0; j < nrow; ++j)
    for (int i = IA[j]; i < IA[j+1]; ++i)
      MatSetValue(pA, JA[i], j, A[i], INSERT_VALUES);

  MatAssemblyBegin(pA,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pA,MAT_FINAL_ASSEMBLY);

  Vec B1, x;
  VecCreate(PETSC_COMM_SELF, &B1);
  VecCreate(PETSC_COMM_SELF, &x);
  VecSetSizes(B1, nrow, PETSC_DECIDE);
  VecSetSizes(x, nrow, PETSC_DECIDE);
  VecSetFromOptions(B1);
  VecSetFromOptions(x);

  size_t nrhs = B.dim() / nrow;
  PetscScalar* bv;
  VecGetArray(B.getVector(), &bv);
  for (size_t i = 0; i < nrhs; ++i) {
    for (size_t j = 0; j < nrow; ++j)
      VecSetValue(B1, j, bv[j*nrow+i], INSERT_VALUES);

    VecAssemblyBegin(B1);
    VecAssemblyEnd(B1);

    if (!this->solve(B1, x, false))
      return false;
    PetscScalar* aa;
    VecGetArray(x, &aa);
    std::copy(aa, aa+nrow, B.getPtr()+i*nrow);
    std::copy(aa, aa+nrow, bv+i*nrow);
    VecRestoreArray(x, &aa);
  }
  VecRestoreArray(B.getVector(), &bv);

  VecDestroy(&x);
  VecDestroy(&B1);

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

  MatAssemblyBegin(pA,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pA,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(B.pA,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B.pA,MAT_FINAL_ASSEMBLY);

  EPSSetOperators(eps,pA,B.pA);
  EPSSetProblemType(eps,EPS_GHEP);
  EPSSetType(eps,EPSKRYLOVSCHUR);
  EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE);
  EPSGetST(eps,&st);
  STSetShift(st,shift);
  EPSSetDimensions(eps,nv,4*nv,PETSC_DETERMINE);
  EPSSetFromOptions(eps);
  EPSSolve(eps);
  EPSGetConverged(eps,&nconv);

  MatGetSize(pA,&m,&n);
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
  MatNorm(pA,NORM_INFINITY,&norm);
  return norm;
}


bool PETScMatrix::setParameters(PETScMatrix* P, PETScVector* Pb)
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
    solParams.setupPC(pc, 0, "", IntSet());
  else if (matvec.size() > 4) {
    std::cerr << "** PETSCMatrix ** Only two blocks supported for now." << std::endl;
    return false;
  }
  else {
    PCSetType(pc,PCFIELDSPLIT);
    PetscInt nsplit;
    KSP  *subksp;
    PC   subpc[2];

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
    PCSetUp(pc);
    PCFieldSplitGetSubKSP(pc,&nsplit,&subksp);

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
      if (solParams.getBlock(m).getStringValue("pc") == "schur")
        new PETScSchurPC(subpc[m], matvec, solParams.getBlock(m), adm);
      else
        solParams.setupPC(subpc[m], m, prefix, adm.dd.getBlockEqs(m));
    }
  }

  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  if (solParams.getIntValue("verbosity") >= 1)
    KSPView(ksp, PETSC_VIEWER_STDOUT_(*adm.getCommunicator()));

  return true;
}


PETScVector operator*(const SystemMatrix& A, const PETScVector& b)
{
  PETScVector results(b.getAdm());
  A.multiply(b, results);
  return results;
}


PETScVector operator/(SystemMatrix& A, const PETScVector& b)
{
  PETScVector results(b.getAdm());
  A.solve(b, results);
  return results;
}
