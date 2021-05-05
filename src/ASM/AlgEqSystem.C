// $Id$
//==============================================================================
//!
//! \file AlgEqSystem.C
//!
//! \date Nov 11 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Storage of an algebraic equation system for a FEM problem.
//!
//==============================================================================

#include "AlgEqSystem.h"
#include "SparseMatrix.h"
#include "ElmMats.h"
#include "SAM.h"
#include "IFEM.h"
#ifdef USE_OPENMP
#include <omp.h>
#endif


AlgEqSystem::AlgEqSystem (const SAM& s, const ProcessAdm* a) : sam(s), adm(a)
{
  d = &c;
}


bool AlgEqSystem::init (LinAlg::MatrixType mtype, const LinSolParams* spar,
                        size_t nmat, size_t nvec, size_t nscl,
                        bool withReactions, int num_threads_SLU)
{
  // Using the sign of the num_threads_SLU argument to flag this (convenience)
  bool dontLockSparsityPattern = num_threads_SLU < 0;

  size_t i;
  for (i = nmat; i < A.size(); i++)
    delete A[i]._A;

  for (i = nvec; i < b.size(); i++)
    delete b[i];

  A.resize(nmat);
  b.resize(nvec,nullptr);
  c.resize(nscl,0.0);
  R.clear();

  for (i = 0; i < A.size(); i++)
  {
    if (!A[i]._A)
    {
      if (spar)
        A[i]._A = SystemMatrix::create(adm,mtype,*spar);
      else
        A[i]._A = SystemMatrix::create(adm,mtype,abs(num_threads_SLU));
      if (!A[i]._A) return false;
    }

    A[i]._A->initAssembly(sam,dontLockSparsityPattern);
    A[i]._b = nullptr;
  }

  for (i = 0; i < b.size(); i++)
    if (!b[i])
    {
      b[i] = SystemVector::create(adm,mtype);
      if (!b[i]) return false;
    }

  bool ok = true;
  if (A.size() == 1 && !b.empty())
  {
    A.front()._b = b.front();
    ok = sam.initForAssembly(*b.front(), withReactions ? &R : nullptr);
  }

  for (i = 0; i < b.size(); i++)
    b[i]->redim(sam.getNoEquations());

  return ok;
}


void AlgEqSystem::clear ()
{
  size_t i;

  for (i = 0; i < A.size(); i++)
    delete A[i]._A;

  for (i = 0; i < b.size(); i++)
    delete b[i];

  if (d && d != &c)
    delete[] d;
  d = nullptr;

  A.clear();
  b.clear();
  c.clear();
  R.clear();
}


bool AlgEqSystem::setAssociatedVector (size_t imat, size_t ivec)
{
  if (imat < A.size() && ivec < b.size())
    A[imat]._b = b[ivec];
  else
    return false;

  return true;
}


void AlgEqSystem::initialize (bool initLHS)
{
  size_t i;

  if (initLHS)
    for (i = 0; i < A.size(); i++)
      A[i]._A->init();

  for (i = 0; i < b.size(); i++)
    b[i]->init();

  for (i = 0; i < c.size(); i++)
    c[i] = 0.0;

#ifdef USE_OPENMP
  size_t nthread = omp_get_max_threads();
  if (nthread > 1 && !c.empty())
  {
    d = new std::vector<double>[nthread];
    for (i = 0; i < nthread; i++)
      d[i].resize(c.size(),0.0);
  }
#endif

  R.fill(0.0);
}


bool AlgEqSystem::assemble (const LocalIntegral* elmObj, int elmId)
{
  const ElmMats* elMat = dynamic_cast<const ElmMats*>(elmObj);
  if (!elMat)
    return false; // Logic error, shouldn't happen...
  else if (elMat->empty())
    return true; // Silently ignore if no element matrices

  size_t i;
  bool status = true;
  if (A.size() == 1 && !b.empty())
  {
    // The algebraic system consists of one system matrix and one RHS-vector.
    // Extract the element-level Newton matrix and associated RHS-vector for
    // general time-dependent and/or nonlinear problems.
    Vector* reac = R.empty() ? nullptr : &R;
    status = sam.assembleSystem(*b.front(), elMat->getRHSVector(), elmId, reac);
#if SP_DEBUG > 2
    for (i = 1; i < b.size() && i < elMat->b.size(); i++)
      std::cout <<"\nElement right-hand-side vector "<< i+1 << elMat->b[i];
#endif

    if (status && elMat->withLHS) // we have LHS element matrices
    {
      if (elMat->rhsOnly) // we only want the RHS system vector
	status = sam.assembleSystem(*b.front(),
				    elMat->getNewtonMatrix(), elmId, reac);
      else // we want both the LHS system matrix and the RHS system vector
	status = sam.assembleSystem(*A.front()._A, *b.front(),
				    elMat->getNewtonMatrix(), elmId, reac);
    }

    // Assembly of additional system right-hand-side vectors
    for (i = 1; i < b.size() && i < elMat->b.size() && status; i++)
      status = sam.assembleSystem(*b[i], elMat->b[i], elmId);
  }
  else
  {
#if SP_DEBUG > 2
    if (elMat->withLHS && !elMat->rhsOnly)
      for (i = 0; i < elMat->A.size() && i < A.size(); i++)
	std::cout <<"Coefficient matrix A"<< i <<" for element "
		  << elmId << elMat->A[i] << std::endl;
    for (i = 0; i < elMat->b.size() && i < b.size(); i++)
      std::cout <<"Right-hand-side vector b"<< i <<" for element "
		<< elmId << elMat->b[i] << std::endl;
#endif

    // Assembly of system right-hand-side vectors
    for (i = 0; i < b.size() && i < elMat->b.size() && status; i++)
      status = sam.assembleSystem(*b[i], elMat->b[i], elmId);

    // Assembly of system coefficient matrices, possibly with right-hand-side
    // contributions too, due to multi-point constraints
    if (elMat->withLHS)
      for (i = 0; i < A.size() && i < elMat->A.size() && status; i++)
	if (A[i]._b)
	{
	  if (elMat->rhsOnly) // we only want the RHS system vectors
	    status = sam.assembleSystem(*A[i]._b, elMat->A[i], elmId);
	  else // we want both LHS system matrices and RHS system vectors
	    status = sam.assembleSystem(*A[i]._A, *A[i]._b, elMat->A[i], elmId);
	}
	else if (!elMat->rhsOnly) // we want LHS system matrices only
	  status = sam.assembleSystem(*A[i]._A, elMat->A[i], elmId);
  }

#if SP_DEBUG > 2
  for (i = 0; i < elMat->c.size() && i < c.size(); i++)
    std::cout <<"Scalar "<< i <<" for element "
              << elmId <<": "<< elMat->c[i] << std::endl;
#endif

  // Assembly of scalar quantities
  size_t it = 0;
#ifdef USE_OPENMP
  it = omp_get_thread_num();
#endif
  for (i = 0; i < c.size() && i < elMat->c.size(); i++)
    d[it][i] += elMat->c[i];

  if (!status)
    std::cerr <<" *** AlgEqSystem::assemble: Failure for element "<< elmId
	      <<": size(A)="<< A.size() <<","<< elMat->A.size()
	      <<" size(b)="<< b.size() <<","<< elMat->b.size() << std::endl;
  return status;
}


bool AlgEqSystem::finalize (bool newLHS)
{
  // Communication of matrix and vector assembly (for PETSc matrices only)
  if (newLHS)
    for (size_t i = 0; i < A.size(); i++)
      if (!A[i]._A->beginAssembly())
	return false;
      else if (!A[i]._A->endAssembly())
	return false;
#if SP_DEBUG > 2
      else if (A[i]._A->dim() < 100)
	std::cout <<"\nSystem coefficient matrix:"<< *A[i]._A;
#endif

  for (size_t i = 0; i < b.size(); i++)
    if (!b[i]->beginAssembly())
      return false;
    else if (!b[i]->endAssembly())
      return false;
#if SP_DEBUG > 2
    else
      std::cout <<"\nSystem right-hand-side vector:"<< *b[i];
#endif

  if (c.empty()) return true;

#ifdef USE_OPENMP
  size_t nthread = omp_get_max_threads();
  if (nthread > 1)
  {
    for (size_t i = 0; i < nthread; i++)
      for (size_t j = 0; j < c.size(); j++)
        c[j] += d[i][j];

    delete[] d;
    d = nullptr;
  }
#endif
#if SP_DEBUG > 2
  std::cout <<"\nScalar quantities:";
  for (size_t i = 0; i < c.size(); i++)
    std::cout <<" "<< c[i];
  std::cout << std::endl;
#endif

  return true;
}


bool AlgEqSystem::staticCondensation (Matrix& Ared, Vector& bred,
                                      const IntVec& extNodes, size_t imat) const
{
  if (imat > A.size()) return false;

  const SparseMatrix* Amat = dynamic_cast<const SparseMatrix*>(A[imat]._A);
  const double*       Rvec = A[imat]._b ? A[imat]._b->getRef() : nullptr;
  if (!Amat || !Rvec) return false;

  // Find the equation numbers to retain
  IntVec mnen, meqn2;
  for (int inod : extNodes)
    if (sam.getNodeEqns(mnen,inod))
      meqn2.insert(meqn2.end(),mnen.begin(),mnen.end());
    else
      return false;

  // Check that all retain DOFs are free
  size_t nspdof = 0;
  for (int ieq : meqn2)
    if (ieq <= 0) ++nspdof;
  if (nspdof > 0)
  {
    std::cerr <<" *** AlgEqSystem::staticCondensation: There are "<< nspdof
              <<" specified DOFs among the DOFs to retain."<< std::endl;
    return false;
  }

  // Extract the sub-matrices associated with internal (1) and external (2) DOFs
  IFEM::cout <<"\nExtracting sub-matrices A11, A12, A22 ..."<< std::endl;
  std::sort(meqn2.begin(),meqn2.end());
  std::array<SparseMatrix,4> Asub;
  // Asub = { A11, A21, A12, A22 }
  if (!Amat->split(Asub,meqn2))
    return false;

  size_t neq0 = Amat->rows();
  size_t neq1 = Asub[0].rows();
  size_t neq2 = Asub[1].rows();

  // Extract the associated sub-vectors
  StdVector R1(neq1); // Internal DOFs
  Vector    R2(neq2); // External (retained) DOFs
  size_t ieq, ip2, ieq1 = 0, ieq2 = 0;
  for (ieq = ip2 = 0; ieq < neq0; ieq++)
    if ((int)ieq+1 < meqn2[ip2])
      R1(++ieq1) = Rvec[ieq];
    else
    {
      R2(++ieq2) = Rvec[ieq];
      ++ip2;
    }

  IFEM::cout <<"\nPerforming static condensation ["<< neq0 <<"x"<< neq0
             <<"] --> ["<< neq2 <<"x"<< neq2 <<"]"<< std::endl;

  // Calculate [A11]^-1*{R1} => R1
  if (!Asub[0].solve(R1,true))
    return false;

  // Calculate [A21]*([A11]^-1*{R1}) => Rtmp
  StdVector Rtmp;
  if (!Asub[1].multiply(R1,Rtmp))
    return false;

  // Calculate the reduced right-hand-side vector
  bred = R2 - Rtmp;

  // Initialize the reduced coefficient matrix
  Ared.resize(neq2,neq2);

  // Loop over the external DOFs
  for (ieq2 = 1; ieq2 <= neq2; ieq2++)
  {
    // Extract the ieq2'th column of [A12] => R1
    Asub[2].getColumn(ieq2,R1);
    // Extract the ieq2'th column of [A22] => R2
    Asub[3].getColumn(ieq2,R2);

    // Calculate [A11]^-1*{R1} => R1, reusing the factorization of [A11]
    if (!Asub[0].solve(R1,false))
      return false;

    // Calculate [A21]*([A11]^-1*{R1}) => Rtmp
    if (!Asub[1].multiply(R1,Rtmp))
      return false;

    // Insert {R2}-{Rtmp} into the reduced matrix Ared
    Ared.fillColumn(ieq2,R2-Rtmp);
  }

  return true;
}
