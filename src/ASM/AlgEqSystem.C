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
#include "ElmMats.h"
#include "SAM.h"


bool AlgEqSystem::init (SystemMatrix::Type mtype, const LinSolParams* spar,
			size_t nmat, size_t nvec, bool withReactions,
			LinAlg::LinearSystemType ltype, int num_threads_SLU)
{
  // Using the sign of the num_threads_SLU argument to flag this (convenience)
  bool dontLockSparsityPattern = num_threads_SLU < 0;

  size_t i;
  for (i = nmat; i < A.size(); i++)
    if (A[i]._A) delete A[i]._A;

  for (i = nvec; i < b.size(); i++)
    if (b[i]) delete b[i];

  A.resize(nmat);
  b.resize(nvec,nullptr);
  R.clear();

  for (i = 0; i < A.size(); i++)
  {
    if (!A[i]._A)
    {
      if (spar)
        A[i]._A = SystemMatrix::create(adm,mtype,*spar,ltype);
      else
        A[i]._A = SystemMatrix::create(adm,mtype,ltype,abs(num_threads_SLU));
      if (!A[i]._A) return false;
    }

    A[i]._A->initAssembly(sam,dontLockSparsityPattern);
    A[i]._b = nullptr;
  }

  SystemVector::Type vtype = SystemVector::STD;
  if (mtype == SystemMatrix::PETSC)
    vtype = SystemVector::PETSC;
  if (mtype == SystemMatrix::ISTL)
    vtype = SystemVector::ISTL;

  for (i = 0; i < b.size(); i++)
    if (!b[i])
    {
      b[i] = SystemVector::create(adm,vtype);
      if (!b[i]) return false;
    }

  bool ok = true;
  if (A.size() == 1 && !b.empty())
    ok = sam.initForAssembly(*b.front(), withReactions ? &R : nullptr);

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

  A.clear();
  b.clear();
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

  return true;
}
