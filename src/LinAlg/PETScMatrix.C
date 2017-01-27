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
#include "PETScSolParams.h"
#include "LinSolParams.h"
#include "LinAlgInit.h"
#include "SAMpatchPETSc.h"
#include "ProcessAdm.h"
#include "SIMenums.h"
#include "ASMstruct.h"
#include "DomainDecomposition.h"
#include "Utilities.h"
#include <cassert>


PETScVector::PETScVector(const ProcessAdm& padm) : adm(padm)
{
  VecCreate(*padm.getCommunicator(),&x);
  VecSetFromOptions(x);
  LinAlgInit::increfs();
}


PETScVector::PETScVector(const ProcessAdm& padm, size_t n) :
  StdVector(n), adm(padm)
{
  VecCreate(*adm.getCommunicator(),&x);
  VecSetSizes(x,adm.dd.getMaxEq()-adm.dd.getMinEq()+1,PETSC_DECIDE);
  VecSetFromOptions(x);
  LinAlgInit::increfs();
}


PETScVector::PETScVector(const ProcessAdm& padm, const Real* values, size_t n) :
  StdVector(values, n), adm(padm)
{
  VecCreate(*adm.getCommunicator(),&x);
  VecSetSizes(x,adm.dd.getMaxEq()-adm.dd.getMinEq() + 1,PETSC_DECIDE);
  VecSetFromOptions(x);
  this->restore(values);
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


bool PETScVector::beginAssembly()
{
  for (size_t i = 0; i < size(); ++i)
    VecSetValue(x , adm.dd.getGlobalEq(i+1)-1, (*this)[i], ADD_VALUES);

  VecAssemblyBegin(x);
  return true;
}


bool PETScVector::endAssembly()
{
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


PETScMatrix::PETScMatrix (const ProcessAdm& padm, const LinSolParams& spar,
                          LinAlg::LinearSystemType ltype) :
 SparseMatrix(SUPERLU, 1), nsp(nullptr), adm(padm), solParams(spar, adm),
 linsysType(ltype)
{
  // Create matrix object, by default the matrix type is AIJ
  MatCreate(*adm.getCommunicator(),&A);

  // Create linear solver object
  KSPCreate(*adm.getCommunicator(),&ksp);

  LinAlgInit::increfs();

  if (spar.getNoBlocks() > 1) {
    matvec.resize(spar.getNoBlocks()*spar.getNoBlocks());
    for (auto& it : matvec) {
      MatCreate(*adm.getCommunicator(), &it);
      MatSetFromOptions(it);
    }
  }

  setParams = true;
  ISsize = 0;
  nLinSolves = 0;
}


PETScMatrix::~PETScMatrix ()
{
  // Deallocation of linear solver object.
  KSPDestroy(&ksp);

  // Deallocation of matrix object.
  MatDestroy(&A);
  LinAlgInit::decrefs();
  for (auto& it : matvec)
    MatDestroy(&it);

  for (auto& it : isvec)
    ISDestroy(&it);

  matvec.clear();
}


void PETScMatrix::initAssembly (const SAM& sam, bool b)
{
  SparseMatrix::initAssembly(sam, b);
  SparseMatrix::preAssemble(sam, b);

  const SAMpatchPETSc* samp = dynamic_cast<const SAMpatchPETSc*>(&sam);
  if (!samp)
    return;

  // Get number of local equations in linear system
  const PetscInt neq  = adm.dd.getMaxEq()- adm.dd.getMinEq() + 1;

  // Set correct number of rows and columns for matrix.
  MatSetSizes(A,neq,neq,PETSC_DECIDE,PETSC_DECIDE);

  // Allocate sparsity pattern
  std::vector<std::set<int>> dofc;
  sam.getDofCouplings(dofc);

  if (matvec.empty()) {
    // Set correct number of rows and columns for matrix.
    MatSetSizes(A,neq,neq,PETSC_DECIDE,PETSC_DECIDE);

    MatSetFromOptions(A);

    // Allocate sparsity pattern
    if (adm.isParallel()) {
      int ifirst = adm.dd.getMinEq();
      int ilast  = adm.dd.getMaxEq();
      PetscIntVec d_nnz(neq, 0);
      IntVec o_nnz_g(adm.dd.getNoGlbEqs(), 0);
      for (int i = 0; i < samp->getNoEquations(); ++i) {
        int eq = adm.dd.getGlobalEq(i+1);
        if (eq >= adm.dd.getMinEq() && eq <= adm.dd.getMaxEq()) {
          for (const auto& it : dofc[i]) {
            int g = adm.dd.getGlobalEq(it);
            if (g > 0) {
              if (g < adm.dd.getMinEq() || g > adm.dd.getMaxEq())
                ++o_nnz_g[eq-1];
              else
                ++d_nnz[eq-adm.dd.getMinEq()];
            }
          }
        } else
          o_nnz_g[eq-1] += dofc[i].size();
      }

      adm.allReduceAsSum(o_nnz_g);

      PetscIntVec o_nnz(o_nnz_g.begin()+ifirst-1, o_nnz_g.begin()+ilast);

      // TODO: multiplier cause big overallocation due to no multiplicity handling
      for (auto& it : o_nnz)
        it = std::min(it, adm.dd.getNoGlbEqs());

      MatMPIAIJSetPreallocation(A,PETSC_DEFAULT,d_nnz.data(),
                                  PETSC_DEFAULT,o_nnz.data());
    } else {
      PetscIntVec Nnz;
      for (const auto& it : dofc)
        Nnz.push_back(it.size());

      MatSeqAIJSetPreallocation(A,PETSC_DEFAULT,Nnz.data());

      PetscIntVec col;
      for (const auto& it2 : dofc)
        for (const auto& it : it2)
          col.push_back(it-1);

      MatSeqAIJSetColumnIndices(A,&col[0]);
      MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
      MatSetOption(A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    }

    MatSetUp(A);

    if (linsysType == LinAlg::SPD)
      MatSetOption(A, MAT_SPD, PETSC_TRUE);
    if (linsysType == LinAlg::SYMMETRIC)
      MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE);

#ifndef SP_DEBUG
    // Do not abort program for allocation error in release mode
    MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
#endif
  } else {
    const DomainDecomposition& dd = adm.dd;
    size_t blocks = solParams.getNoBlocks();

    // map from sparse matrix indices to block matrix indices
    glb2Blk.resize(SparseMatrix::A.size());
    std::vector<std::array<int,2>> eq2b(sam.getNoEquations(), {-1, 0}); // cache
    for (size_t j = 0; j < cols(); ++j) {
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
          if (iblk == -1 && (it = dd.getG2LEQ(b).find(JA[i]+1)) != dd.getG2LEQ(b).end()) {
            iblk = b;
            eq2b[JA[i]][0] = b;
            eq2b[JA[i]][1] = glb2Blk[i][1] = it->second-1;
          }

          if (jblk == -1 && (it = dd.getG2LEQ(b).find(j+1)) != dd.getG2LEQ(b).end()) {
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
    }

    if (adm.isParallel()) {
      std::vector<PetscIntVec> d_nnz(blocks*blocks);
      std::vector<IntVec> o_nnz_g(blocks*blocks);
      size_t k = 0;
      for (size_t i = 0; i < blocks; ++i)
        for (size_t j = 0; j < blocks; ++j, ++k) {
          d_nnz[k].resize(dd.getMaxEq(i+1)-dd.getMinEq(i+1)+1);
          o_nnz_g[k].resize(dd.getNoGlbEqs(i+1));
        }

      for (int i = 0; i < sam.getNoEquations(); ++i) {
        int blk = eq2b[i][0]+1;
        int row = eq2b[i][1]+1;
        int grow = dd.getGlobalEq(row, blk);

        if (grow >= dd.getMinEq(blk) && grow <= dd.getMaxEq(blk)) {
          for (const auto& it : dofc[i]) {
            int cblk = eq2b[it-1][0]+1;
            int col = eq2b[it-1][1]+1;
            int gcol = dd.getGlobalEq(col, cblk);
            if (gcol >= dd.getMinEq(cblk) && gcol <= dd.getMaxEq(cblk))
              ++d_nnz[(blk-1)*blocks + cblk-1][grow-dd.getMinEq(blk)];
            else
              ++o_nnz_g[(blk-1)*blocks + cblk-1][grow-1];
          }
        } else {
          for (const auto& it : dofc[i]) {
            int cblk = eq2b[it-1][0]+1;
            ++o_nnz_g[(blk-1)*blocks+cblk-1][grow-1];
          }
        }
      }

      k = 0;
      for (size_t i = 0; i < blocks; ++i) {
        for (size_t j = 0; j < blocks; ++j, ++k) {
          int nrows = dd.getMaxEq(i+1)-dd.getMinEq(i+1)+1;
          int ncols = dd.getMaxEq(j+1)-dd.getMinEq(j+1)+1;
          MatSetSizes(matvec[k], nrows, ncols,
                      PETSC_DETERMINE, PETSC_DETERMINE);
          MatSetFromOptions(matvec[k]);
          adm.allReduceAsSum(o_nnz_g[k]);
          PetscIntVec o_nnz(o_nnz_g[k].begin()+dd.getMinEq(i+1)-1, o_nnz_g[k].begin()+dd.getMaxEq(i+1));

          // TODO: multiplier cause big overallocation due to no multiplicity handling
          for (auto& it : o_nnz)
            it = std::min(it, dd.getNoGlbEqs(j+1));

          MatMPIAIJSetPreallocation(matvec[k],PETSC_DEFAULT,d_nnz[k].data(),
                                              PETSC_DEFAULT,o_nnz.data());
          MatSetUp(matvec[k]);
        }
      }
    } else {
      auto it = matvec.begin();
      for (size_t i = 0; i < blocks; ++i) {
        for (size_t j = 0; j < blocks; ++j, ++it) {
          std::vector<PetscInt> nnz;
          nnz.reserve(dd.getBlockEqs(i).size());
          for (const auto& it2 : dd.getBlockEqs(i))
            nnz.push_back(std::min(dofc[it2-1].size(), dd.getBlockEqs(j).size()));

          int nrows = dd.getMaxEq(i+1)-dd.getMinEq(i+1)+1;
          int ncols = dd.getMaxEq(j+1)-dd.getMinEq(j+1)+1;
          MatSetSizes(*it, nrows, ncols,
                      PETSC_DETERMINE, PETSC_DETERMINE);

          MatSeqAIJSetPreallocation(*it, PETSC_DEFAULT, nnz.data());
          MatSetUp(*it);
        }
      }
    }

    isvec.resize(dd.getNoBlocks());
    // index sets
    for (size_t i = 0; i < isvec.size(); ++i) {
      std::vector<int> blockEq;
      blockEq.reserve(dd.getMaxEq(i+1)-dd.getMinEq(i+1)+1);
      for (auto& it : dd.getBlockEqs(i)) {
        int eq = dd.getGlobalEq(it);
        if (eq >= dd.getMinEq())
          blockEq.push_back(eq-1);
      }

      ISCreateGeneral(*adm.getCommunicator(),blockEq.size(),
                      blockEq.data(),PETSC_COPY_VALUES,&isvec[i]);
    }

    MatCreateNest(*adm.getCommunicator(),solParams.getNoBlocks(),isvec.data(),
                  solParams.getNoBlocks(),isvec.data(),matvec.data(),&A);

 #ifndef SP_DEBUG
    // Do not abort program for allocation error in release mode
    for (auto& it : matvec)
      MatSetOption(it,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
 #endif
  }
}


bool PETScMatrix::beginAssembly()
{
  if (matvec.empty()) {
    for (size_t j = 0; j < cols(); ++j)
      for (int i = IA[j]; i < IA[j+1]; ++i)
        MatSetValue(A, adm.dd.getGlobalEq(JA[i]+1)-1,
                    adm.dd.getGlobalEq(j+1)-1, SparseMatrix::A[i], ADD_VALUES);
  } else {
    for (size_t j = 0; j < cols(); ++j) {
      for (int i = IA[j]; i < IA[j+1]; ++i) {
        int rblock = glb2Blk[i][0] / adm.dd.getNoBlocks() + 1;
        int cblock = glb2Blk[i][0] % adm.dd.getNoBlocks() + 1;
        MatSetValue(matvec[glb2Blk[i][0]],
                    adm.dd.getGlobalEq(glb2Blk[i][1]+1, rblock)-1,
                    adm.dd.getGlobalEq(glb2Blk[i][2]+1, cblock)-1,
                    SparseMatrix::A[i], ADD_VALUES);
      }
    }
  }
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
  SparseMatrix::init();

  // Set all matrix elements to zero
  if (matvec.empty())
    MatZeroEntries(A);
  else {
    for (auto& it : matvec)
      MatZeroEntries(it);
  }
}


bool PETScMatrix::multiply (const SystemVector& B, SystemVector& C) const
{
  const PETScVector* Bptr = dynamic_cast<const PETScVector*>(&B);
        PETScVector* Cptr = dynamic_cast<PETScVector*>(&C);

  if ((!Bptr) || (!Cptr))
    return false;

  MatMult(A,Bptr->getVector(),Cptr->getVector());
  return true;
}


bool PETScMatrix::solve (SystemVector& B, bool newLHS, Real*)
{
  PETScVector* Bptr = dynamic_cast<PETScVector*>(&B);
  if (!Bptr)
    return false;

  StdVector* Bsptr = dynamic_cast<StdVector*>(&B);
  if (!Bsptr)
    return false;

  Vec x;
  VecDuplicate(Bptr->getVector(),&x);
  VecCopy(Bptr->getVector(),x);

  bool result = this->solve(x,Bptr->getVector(),newLHS,true);
  VecDestroy(&x);

  return result;
}


bool PETScMatrix::solve (const SystemVector& b, SystemVector& x, bool newLHS)
{
  const PETScVector* Bptr = dynamic_cast<const PETScVector*>(&b);
  if (!Bptr)
    return false;

  PETScVector* Xptr = dynamic_cast<PETScVector*>(&x);
  if (!Xptr)
    return false;

  return this->solve(Bptr->getVector(),Xptr->getVector(),newLHS,false);
}


bool PETScMatrix::solve (const Vec& b, Vec& x, bool newLHS, bool knoll)
{
  // Reset linear solver
  if (nLinSolves && solParams.getIntValue("reset_solves"))
    if (nLinSolves%solParams.getIntValue("reset_solves") == 0) {
      KSPDestroy(&ksp);
      KSPCreate(*adm.getCommunicator(),&ksp);
      setParams = true;
    }

  if (setParams) {
#if PETSC_VERSION_MINOR < 5
    KSPSetOperators(ksp,A,A, newLHS ? SAME_NONZERO_PATTERN : SAME_PRECONDITIONER);
#else
    KSPSetOperators(ksp,A,A);
    KSPSetReusePreconditioner(ksp, newLHS ? PETSC_FALSE : PETSC_TRUE);
#endif
    if (!setParameters())
      return false;
    setParams = false;
  }
  if (knoll)
    KSPSetInitialGuessKnoll(ksp,PETSC_TRUE);
  else
    KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
  KSPSolve(ksp,b,x);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp,&reason);
  if (reason < 0) {
    PetscPrintf(PETSC_COMM_WORLD, "\n Linear solve failed with reason %s",KSPConvergedReasons[reason]);
    return false;
  }

  if (solParams.getIntValue("verbosity") > 1) {
    PetscInt its;
    KSPGetIterationNumber(ksp,&its);
    PetscPrintf(PETSC_COMM_WORLD,"\n Iterations for %s = %D\n",solParams.getStringValue("type").c_str(),its);
  }
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


bool PETScMatrix::setParameters(PETScMatrix* P, PETScVector* Pb)
{
  // Set linear solver method
  KSPSetType(ksp,solParams.getStringValue("type").c_str());
  KSPSetTolerances(ksp,solParams.getDoubleValue("rtol"),
                   solParams.getDoubleValue("atol"),
                   solParams.getDoubleValue("dtol"),
                   solParams.getIntValue("maxits"));
  PC pc;
  KSPGetPC(ksp,&pc);

  if (matvec.empty()) {
    solParams.setupPC(pc, 0, "", std::set<int>());
  } else {
    if (matvec.size() > 4) {
      std::cerr << "** PETSCMatrix ** Only two blocks supported for now." << std::endl;
      return false;
    }
    solParams.setupSchurComplement(matvec);
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

//    MatCreateSchurComplement(matvec[0],matvec[0],matvec[1],matvec[2],matvec[3],&Sp);
    PCSetFromOptions(pc);
    PCSetUp(pc);
    PCFieldSplitGetSubKSP(pc,&nsplit,&subksp);
#if PETSC_VERSION_MINOR < 5
    KSPSetOperators(subksp[1],solParams.getSchurComplement(),solParams.getSchurComplement(),SAME_PRECONDITIONER);
#else
    KSPSetOperators(subksp[1],solParams.getSchurComplement(),solParams.getSchurComplement());
    KSPSetReusePreconditioner(subksp[1], PETSC_TRUE);
#endif

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
