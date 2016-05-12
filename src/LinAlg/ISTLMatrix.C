// $Id$
//==============================================================================
//!
//! \file ISTLMatrix.C
//!
//! \date Mar 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Representation of the system matrix in ISTL format.
//!
//==============================================================================

#include "ISTLMatrix.h"
#include "LinSolParams.h"
#include "LinAlgInit.h"
#include "ProcessAdm.h"
#include "SIMenums.h"
#include "ASMstruct.h"
#include "Profiler.h"
#include "SAMpatch.h"
#include "DomainDecomposition.h"


ISTLVector::ISTLVector(const ProcessAdm& padm) : adm(padm)
{
  LinAlgInit::increfs();
}


ISTLVector::ISTLVector(const ProcessAdm& padm, size_t n) : adm(padm)
{
  x.resize(n);
  LinAlgInit::increfs();
}


ISTLVector::ISTLVector(const ProcessAdm& padm, const Real* values, size_t n) : adm(padm)
{
  x.resize(n);
  this->restore(values);
  LinAlgInit::increfs();
}


ISTLVector::ISTLVector(const ISTLVector& vec) : adm(vec.adm)
{
  x = vec.x;
  LinAlgInit::increfs();
}


ISTLVector::~ISTLVector()
{
  LinAlgInit::decrefs();
}


void ISTLVector::init(Real value)
{
  StdVector::init(value);
  x = value;
}

size_t ISTLVector::dim() const
{
  return x.size();
}


void ISTLVector::redim(size_t n)
{
  x.resize(n);
  StdVector::redim(n);
}


bool ISTLVector::beginAssembly()
{
  for (size_t i = 0; i < size(); ++i)
    x[i] = (*this)[i];

  return true;
}


bool ISTLVector::endAssembly()
{
  return true;
}


Real ISTLVector::L1norm() const
{
  return x.one_norm();
}


Real ISTLVector::L2norm() const
{
  return x.two_norm();
}


Real ISTLVector::Linfnorm() const
{
  return x.infinity_norm();
}


ISTLMatrix::ISTLMatrix (const ProcessAdm& padm, const LinSolParams& spar,
                        LinAlg::LinearSystemType ltype) :
 SparseMatrix(SUPERLU, 1), adm(padm),
 solParams(spar, adm), linsysType(ltype)
{
  LinAlgInit::increfs();

  setParams = true;
  nLinSolves = 0;
}


ISTLMatrix::ISTLMatrix (const ISTLMatrix& B) :
  adm(B.adm), solParams(B.solParams.get(), B.adm)
{
  A = B.A;

  LinAlgInit::increfs();

  setParams = true;
}


ISTLMatrix::~ISTLMatrix ()
{
  LinAlgInit::decrefs();
}


void ISTLMatrix::initAssembly (const SAM& sam, bool b)
{
  SparseMatrix::initAssembly(sam, b);
  SparseMatrix::preAssemble(sam, b);

  std::vector<std::set<int>> dofc;
  sam.getDofCouplings(dofc);

  // Set correct number of rows and columns for matrix.
  size_t sum = 0;
  for (const auto& it : dofc)
    sum += it.size();

  A.setSize(rows(), cols(), sum);
  A.setBuildMode(ISTL::Mat::random);

  for (size_t i = 0; i < dofc.size(); ++i)
    A.setrowsize(i,dofc[i].size());
  A.endrowsizes();

  for (size_t i = 0; i < dofc.size(); ++i)
    for (const auto& it : dofc[i])
      A.addindex(i, it-1);

  A.endindices();

  A = 0;
}

bool ISTLMatrix::beginAssembly()
{
  for (size_t j = 0; j < cols(); ++j)
    for (int i = IA[j]; i < IA[j+1]; ++i)
      A[JA[i]][j] = SparseMatrix::A[i];

  return true;
}


bool ISTLMatrix::endAssembly()
{
  return true;
}


void ISTLMatrix::init ()
{
  SparseMatrix::init();

  // Set all matrix elements to zero
  A = 0;
}



bool ISTLMatrix::solve (SystemVector& B, bool newLHS, Real*)
{
  if (!pre)
    std::tie(solver, pre, op) = solParams.setupPC(A);

  ISTLVector* Bptr = dynamic_cast<ISTLVector*>(&B);
  if (!Bptr || !solver || !pre)
    return false;

  try {
    Dune::InverseOperatorResult r;
    ISTL::Vec b(Bptr->getVector());
    Bptr->getVector() = 0;
    solver->apply(Bptr->getVector(), b, r);
  } catch (Dune::ISTLError& e) {
    std::cerr << "ISTL exception " << e << std::endl;
    return false;
  }

  for (size_t i = 0; i < rows(); ++i)
    (*Bptr)(i+1) = Bptr->getVector()[i];

  return true;
}


bool ISTLMatrix::solve (const SystemVector& b, SystemVector& x, bool newLHS)
{
  if (!pre)
    std::tie(solver, pre, op) = solParams.setupPC(A);

  const ISTLVector* Bptr = dynamic_cast<const ISTLVector*>(&b);
  if (!Bptr || ! solver || !pre)
    return false;

  ISTLVector* Xptr = dynamic_cast<ISTLVector*>(&x);
  if (!Xptr)
    return false;

  try {
    Dune::InverseOperatorResult r;
    solver->apply(Xptr->getVector(),
                  const_cast<ISTL::Vec&>(Bptr->getVector()), r);
  } catch (Dune::ISTLError& e) {
    std::cerr << "ISTL exception " << e << std::endl;
    return false;
  }

  for (size_t i = 0; i < rows(); ++i)
    (*Xptr)(i+1) = Xptr->getVector()[i];

  return true;
}


Real ISTLMatrix::Linfnorm () const
{
  return A.infinity_norm();
}
