// $Id: AlgEqSystem.C,v 1.13 2011-02-08 09:32:18 kmo Exp $
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
#include "LinSolParams.h"


void AlgEqSystem::init (SystemMatrix::Type mtype, const LinSolParams* spar,
			size_t nmat, size_t nvec, int num_threads_SLU)
{
  size_t i;
  for (i = nmat; i < A.size(); i++)
    if (A[i]._A) delete A[i]._A;

  for (i = nvec; i < b.size(); i++)
    if (b[i]) delete b[i];

  A.resize(nmat);
  b.resize(nvec,0);
  R.clear();

  for (i = 0; i < A.size(); i++)
  {
    if (!A[i]._A)
      if (spar)
	A[i]._A = SystemMatrix::create(mtype,*spar);
      else
	A[i]._A = SystemMatrix::create(mtype,num_threads_SLU);
    A[i]._b = 0;
  }

  for (i = 0; i < b.size(); i++)
    if (!b[i])
      if (mtype == SystemMatrix::PETSC)
	b[i] = SystemVector::create(SystemVector::PETSC);
      else
	b[i] = SystemVector::create(SystemVector::STD);
}


void AlgEqSystem::init (bool initLHS)
{
  size_t i;

  if (initLHS)
    for (i = 0; i < A.size(); i++)
      if (A[i]._A) A[i]._A->init();

  for (i = 0; i < b.size(); i++)
    if (b[i]) b[i]->init();

  R.fill(0.0);
}


void AlgEqSystem::clear ()
{
  size_t i;

  for (i = 0; i < A.size(); i++)
    if (A[i]._A) delete A[i]._A;

  for (i = 0; i < b.size(); i++)
    if (b[i]) delete b[i];

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


void AlgEqSystem::initAssembly ()
{
  size_t i;

  for (i = 0; i < A.size(); i++)
    if (A[i]._A) A[i]._A->initAssembly(sam);

  if (A.size() == 1 && b.size() == 1)
    sam.initForAssembly(*b.front(),&R);
  else for (i = 0; i < b.size(); i++)
    if (b[i])
      b[i]->redim(sam.getNoEquations());
}


bool AlgEqSystem::assemble (const LocalIntegral* elmObj, int elmId)
{
  const ElmMats* elMat = dynamic_cast<const ElmMats*>(elmObj);
  if (!elMat) return false;

  bool status = true;
  if (A.size() == 1 && b.size() == 1)
  {
    // The algebraic system consists of one system matrix and one RHS-vector.
    // Extract the element-level Newton matrix and associated RHS-vector for
    // general time-dependent and/or nonlinear problems.
    status = sam.assembleSystem(*b.front(), elMat->getRHSVector(), elmId, &R);
    if (status && elMat->withLHS) // we have LHS element matrices
      if (elMat->rhsOnly) // we only want the RHS system vector
	status = sam.assembleSystem(*b.front(),
				    elMat->getNewtonMatrix(), elmId, &R);
      else // we want both the LHS system matrix and the RHS system vector
	status = sam.assembleSystem(*A.front()._A, *b.front(),
				    elMat->getNewtonMatrix(), elmId, &R);
  }
  else
  {
    size_t i;
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
	  if (elMat->rhsOnly) // we only want the RHS system vectors
	    status = sam.assembleSystem(*A[i]._b, elMat->A[i], elmId);
	  else // we want both LHS system matrices and RHS system vectors
	    status = sam.assembleSystem(*A[i]._A, *A[i]._b, elMat->A[i], elmId);
	else if (!elMat->rhsOnly) // we want LHS system matrices only
	  status = sam.assembleSystem(*A[i]._A, elMat->A[i], elmId);
  }

  if (!status)
    std::cerr <<" *** AlgEqSystem::assemble: Failure for element "<< elmId
	      <<": size(A)="<< A.size() <<","<< elMat->A.size()
	      <<" size(b)="<< b.size() <<","<< elMat->b.size() << std::endl;
  return status;
}
