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
#include "SAM.h"


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
 SparseMatrix(SUPERLU, 1), op(A), adm(padm),
 solParams(spar), linsysType(ltype)
{
  LinAlgInit::increfs();

  setParams = true;
  nLinSolves = 0;
}


ISTLMatrix::ISTLMatrix (const ISTLMatrix& B) :
  op(A), adm(B.adm), solParams(B.solParams)
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

  std::vector<std::set<int>> dofc;
  sam.getDofCouplings(dofc);

  // Set correct number of rows and columns for matrix.
  size_t sum = 0;
  for (const auto& it : dofc)
    sum += it.size();

  A.setSize(rows(), cols(), sum);
  A.setBuildMode(Mat::random);

  for (size_t i = 0; i < dofc.size(); ++i)
    A.setrowsize(i,dofc[i].size());
  A.endrowsizes();

  for (size_t i = 0; i < dofc.size(); ++i)
    for (const auto& it : dofc[i])
      A.addindex(i, it-1);
  A.endindices();
}

bool ISTLMatrix::beginAssembly()
{
  this->optimiseSLU();
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


/*! \brief Helper template for setting up a solver with the appropriate
           preconditioner type.
    \details We cannot instance using dynamic polymorphism, the solver need
             access to the real type for the preconditioner. We can however
             call the solver in the interface class scope afterwards.
 */
template<class Prec>
static Dune::InverseOperator<ISTLMatrix::Vec,ISTLMatrix::Vec>*
  setupWithPreType(const LinSolParams& solParams,
                   Dune::MatrixAdapter<ISTLMatrix::Mat,ISTLMatrix::Vec,ISTLMatrix::Vec>& op,
                   Dune::Preconditioner<ISTLMatrix::Vec,ISTLMatrix::Vec>& prec)
{
  Prec& pre = static_cast<Prec&>(prec);

  if (solParams.getMethod() == "bcgs")
    return new Dune::BiCGSTABSolver<ISTLMatrix::Vec>(op, pre, solParams.getRelTolerance(),
                                                     solParams.getMaxIterations(),1);
  else if (solParams.getMethod() == "cg")
    return new Dune::CGSolver<ISTLMatrix::Vec>(op, pre, solParams.getRelTolerance(),
                                               solParams.getMaxIterations(),1);
  else if (solParams.getMethod() == "minres")
    return new Dune::MINRESSolver<ISTLMatrix::Vec>(op, pre, solParams.getRelTolerance(),
                                                   solParams.getMaxIterations(),1);
  else if (solParams.getMethod() == "gmres")
    return new Dune::RestartedGMResSolver<ISTLMatrix::Vec>(op, pre, solParams.getRelTolerance(), 100,
                                                           solParams.getMaxIterations(),1);

  return nullptr;
}


void ISTLMatrix::setupSolver()
{
  if (!pre) {
    if (solParams.getBlock(0).prec == "ilu") {
      if (solParams.getBlock(0).ilu_fill_level == 0) {
        pre.reset(new Dune::SeqILU0<Mat,Vec,Vec>(A, 0.9));
        solver.reset(setupWithPreType<Dune::SeqILU0<Mat,Vec,Vec>>(solParams, op, *pre));
      } else {
        pre.reset(new Dune::SeqILUn<Mat,Vec,Vec>(A, solParams.getBlock(0).ilu_fill_level, 0.9));
        solver.reset(setupWithPreType<Dune::SeqILUn<Mat,Vec,Vec>>(solParams, op, *pre));
      }
    } else if (solParams.getBlock(0).prec == "sor") {
      pre.reset(new Dune::SeqSOR<Mat,Vec,Vec>(A, 1, 1.0));
      solver.reset(setupWithPreType<Dune::SeqSOR<Mat,Vec,Vec>>(solParams, op, *pre));
    } else if (solParams.getBlock(0).prec == "ssor") {
      pre.reset(new Dune::SeqSOR<Mat,Vec,Vec>(A, 1, 1.0));
      solver.reset(setupWithPreType<Dune::SeqSSOR<Mat,Vec,Vec>>(solParams, op, *pre));
    } else if (solParams.getBlock(0).prec == "jacobi") {
      pre.reset(new Dune::SeqJac<Mat,Vec,Vec>(A, 1, 1.0));
      solver.reset(setupWithPreType<Dune::SeqJac<Mat,Vec,Vec>>(solParams, op, *pre));
    } else if (solParams.getBlock(0).prec == "gs") {
      pre.reset(new Dune::SeqGS<Mat,Vec,Vec>(A, 1, 1.0));
      solver.reset(setupWithPreType<Dune::SeqGS<Mat,Vec,Vec>>(solParams, op, *pre));
    }
  }
}


bool ISTLMatrix::solve (SystemVector& B, bool newLHS, Real*)
{
  ISTLVector* Bptr = dynamic_cast<ISTLVector*>(&B);
  if (!Bptr)
    return false;

  setupSolver();

  try {
    Dune::InverseOperatorResult r;
    Vec b(Bptr->getVector());
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
  const ISTLVector* Bptr = dynamic_cast<const ISTLVector*>(&b);
  if (!Bptr)
    return false;

  ISTLVector* Xptr = dynamic_cast<ISTLVector*>(&x);
  if (!Xptr)
    return false;

  setupSolver();

  try {
    Dune::InverseOperatorResult r;
    solver->apply(Xptr->getVector(),
                  const_cast<Vec&>(Bptr->getVector()), r);
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
