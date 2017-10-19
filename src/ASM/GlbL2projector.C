// $Id$
//==============================================================================
//!
//! \file GlbL2projector.C
//!
//! \date Oct 16 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief General integrand for L2-projection of secondary solutions.
//!
//==============================================================================

#include "GlbL2projector.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "FiniteElement.h"
#include "TimeDomain.h"
#include "ASMbase.h"
#include "ElmMats.h"
#include "IntegrandBase.h"
#include "Function.h"
#include "Profiler.h"
#include "SparseMatrix.h"
#ifdef HAS_PETSC
#include "PETScMatrix.h"
#include "LinSolParams.h"
#include "ProcessAdm.h"
#endif


LinAlg::MatrixType GlbL2::MatrixType   = LinAlg::SPARSE;
LinSolParams*      GlbL2::SolverParams = nullptr;


/*!
  \brief Local integral container class for L2-projections.
*/

class L2Mats : public ElmMats
{
public:
  //! \brief The constructor initializes pointers and references.
  //! \param[in] p The global L2 integrand object containing projection matrices
  //! \param[in] nen Number of element nodes
  //! \param[in] nf Number of field components
  //! \param[in] q Pointer to element data associated with the problem integrand
  L2Mats(GlbL2& p, size_t nen, size_t nf, LocalIntegral* q = nullptr)
    : gl2Int(p), elmData(q)
  {
    this->resize(1,nf);
    this->redim(nen);
  }

  //! \brief Empty destructor.
  virtual ~L2Mats() {}

  //! \brief Destruction method to clean up after numerical integration.
  virtual void destruct() { delete elmData; delete this; }

  GlbL2&         gl2Int;  //!< The global L2-projection integrand
  LocalIntegral* elmData; //!< Element data associated with problem integrand

  IntVec  mnpc;        //!< Matrix of element nodal correspondance
  uIntVec elem_sizes;  //!< Size of each basis on the element
  uIntVec basis_sizes; //!< Size of each basis on the patch
};


/*!
  \brief Global integral container class for L2-projections.
*/

class L2GlobalInt : public GlobalIntegral
{
public:
  //! \brief The constructor initializes the system matrix references.
  L2GlobalInt(SparseMatrix& A_, StdVector& B_) : A(A_), B(B_) {}

  //! \brief Empty destructor.
  virtual ~L2GlobalInt() {}

  //! \brief Adds a LocalIntegral object into a corresponding global object.
  virtual bool assemble(const LocalIntegral* elmObj, int)
  {
    const L2Mats* elm = static_cast<const L2Mats*>(elmObj);
    for (size_t i = 0; i < elm->mnpc.size(); i++)
    {
      int inod = elm->mnpc[i]+1;
      for (size_t j = 0; j < elm->mnpc.size(); j++)
      {
        int jnod = elm->mnpc[j]+1;
        A(inod,jnod) += elm->A.front()(i+1,j+1);
      }
      for (const Vector& b : elm->b)
      {
        B(inod) += b[i];
        inod += A.dim();
      }
    }
    return true;
  }

private:
  SparseMatrix& A; //!< Reference to left-hand-side matrix
  StdVector&    B; //!< Reference to right-hand-side vector
};


GlbL2::GlbL2 (IntegrandBase* p, size_t n) : problem(p)
{
  nrhs = p->getNoFields(2);
  this->allocate(n);
}


GlbL2::GlbL2 (FunctionBase* f, size_t n) : problem(nullptr), functions({f})
{
  nrhs = f->dim();
  this->allocate(n);
}


GlbL2::GlbL2 (const FunctionVec& f, size_t n) : problem(nullptr), functions(f)
{
  nrhs = 0;
  for (FunctionBase* func : f)
    nrhs += func->dim();
  this->allocate(n);
}


GlbL2::~GlbL2()
{
  delete pA;
  delete pB;
#ifdef HAS_PETSC
  delete adm;
#endif
}


void GlbL2::allocate (size_t n)
{
#ifdef HAS_PETSC
  adm = nullptr;
  if (GlbL2::MatrixType == LinAlg::PETSC && GlbL2::SolverParams)
  {
    adm = new ProcessAdm();
    pA = new PETScMatrix(*adm,*GlbL2::SolverParams);
    pB = new PETScVector(*adm,n*nrhs);
  }
  else
#endif
  {
    pA = new SparseMatrix(SparseMatrix::SUPERLU);
    pB = new StdVector(n*nrhs);
  }

  pA->redim(n,n);
}


int GlbL2::getIntegrandType () const
{
  if (problem)
    // Mask off the element interface flag
    return problem->getIntegrandType() & ~INTERFACE_TERMS;
  else
    return STANDARD;
}


LocalIntegral* GlbL2::getLocalIntegral (size_t nen, size_t iEl,
                                        bool neumann) const
{
  if (problem)
    return new L2Mats(*const_cast<GlbL2*>(this),nen,nrhs,
                      problem->getLocalIntegral(nen,iEl,neumann));
  else
    return new L2Mats(*const_cast<GlbL2*>(this),nen,nrhs);
}


bool GlbL2::initElement (const IntVec& MNPC, const FiniteElement& fe,
                         const Vec3& Xc, size_t nPt,
                         LocalIntegral& elmInt)
{
  L2Mats& gl2 = static_cast<L2Mats&>(elmInt);

  gl2.mnpc = MNPC;
  if (problem && gl2.elmData)
    return problem->initElement(MNPC,fe,Xc,nPt,*gl2.elmData);
  else
    return true;
}


bool GlbL2::initElement (const IntVec& MNPC1,
                         const uIntVec& elem_sizes, const uIntVec& basis_sizes,
                         LocalIntegral& elmInt)
{
  L2Mats& gl2 = static_cast<L2Mats&>(elmInt);

  gl2.mnpc = MNPC1;
  gl2.elem_sizes = elem_sizes;
  gl2.basis_sizes = basis_sizes;
  if (problem && gl2.elmData)
    return problem->initElement(MNPC1,elem_sizes,basis_sizes,*gl2.elmData);
  else
    return true;
}


bool GlbL2::evalInt (LocalIntegral& elmInt,
                     const FiniteElement& fe,
                     const Vec3& X) const

{
  L2Mats& gl2 = static_cast<L2Mats&>(elmInt);

  Vector solPt;
  solPt.reserve(nrhs);
  if (problem)
  {
    if (!problem->evalSol(solPt,fe,X,gl2.mnpc))
      if (!problem->diverged(fe.iGP+1))
        return false;
  }
  else if (functions.size() == 1)
    solPt = functions.front()->getValue(X);
  else for (FunctionBase* func : functions)
  {
    RealArray funcPt = func->getValue(X);
    solPt.insert(solPt.end(),funcPt.begin(),funcPt.end());
  }

  gl2.A.front().outer_product(fe.N,fe.N,true,fe.detJxW);
  for (size_t j = 0; j < solPt.size(); j++)
    gl2.b[j].add(fe.N,solPt[j]*fe.detJxW);

  return true;
}


bool GlbL2::evalIntMx (LocalIntegral& elmInt,
                       const MxFiniteElement& fe,
                       const Vec3& X) const

{
  if (!problem)
    return this->evalInt(elmInt,fe,X);

  L2Mats& gl2 = static_cast<L2Mats&>(elmInt);

  Vector solPt;
  if (!problem->evalSol(solPt,fe,X,gl2.mnpc,gl2.elem_sizes,gl2.basis_sizes))
    if (!problem->diverged(fe.iGP+1))
      return false;

  gl2.A.front().outer_product(fe.N,fe.N,true,fe.detJxW);
  for (size_t j = 0; j < solPt.size(); j++)
    gl2.b[j].add(fe.N,solPt[j]*fe.detJxW);

  return true;
}


void GlbL2::preAssemble (const std::vector<IntVec>& MMNPC, size_t nel)
{
  pA->preAssemble(MMNPC,nel);
}


bool GlbL2::solve (Matrix& sField)
{
  SparseMatrix& A = *pA;
  StdVector&    B = *pB;

  // Insert a 1.0 value on the diagonal for equations with no contributions.
  // Needed in immersed boundary calculations with "totally outside" elements.
  size_t i, j, nnod = A.dim();
  for (i = 1; i <= nnod; i++)
    if (A(i,i) == 0.0) A(i,i) = 1.0;

#if SP_DEBUG > 1
  std::cout <<"\nGlobal L2-projection matrix:\n"<< A;
  std::cout <<"\nGlobal L2-projection RHS:"<< B;
#endif

  // Solve the patch-global equation system
  if (!A.solve(B)) return false;

  // Store the nodal values of the projected field
  sField.resize(nrhs,nnod);
  for (i = 1; i <= nnod; i++)
    for (j = 1; j <= nrhs; j++)
      sField(j,i) = B(i+(j-1)*nnod);

#if SP_DEBUG > 1
  std::cout <<"\nSolution:"<< sField;
#endif
  return true;
}


bool GlbL2::solve (const std::vector<Matrix*>& sField)
{
  if (sField.size() != functions.size())
  {
    std::cerr <<" *** GlbL2::solve: Logic error, size(sField)="<< sField.size()
              <<" != size(functions)="<< functions.size() << std::endl;
    return false;
  }

  SparseMatrix& A = *pA;
  StdVector&    B = *pB;

  // Insert a 1.0 value on the diagonal for equations with no contributions.
  // Needed in immersed boundary calculations with "totally outside" elements.
  size_t i, j, nnod = A.dim();
  for (i = 1; i <= nnod; i++)
    if (A(i,i) == 0.0) A(i,i) = 1.0;

#if SP_DEBUG > 1
  std::cout <<"\nGlobal L2-projection matrix:\n"<< A;
  std::cout <<"\nGlobal L2-projection RHS:"<< B;
#endif

  // Solve the patch-global equation system
  if (!A.solve(B)) return false;

  // Store the nodal values of the projected fields
  size_t offset = 0;
  for (size_t k = 0; k < sField.size(); k++)
  {
    size_t ncomp = functions[k]->dim();
    sField[k]->resize(ncomp,nnod);
    for (i = 1; i <= nnod; i++)
      for (j = 1; j <= ncomp; j++)
        (*sField[k])(j,i) = B(i+(offset+j-1)*nnod);
    offset += ncomp;

#if SP_DEBUG > 1
    std::cout <<"\nSolution "<< k <<":"<< *sField[k];
#endif
  }

  return true;
}


bool ASMbase::L2projection (Matrix& sField,
                            IntegrandBase* integrand,
                            const TimeDomain& time)
{
  PROFILE2("ASMbase::L2projection");

  GlbL2 gl2(integrand,this->getNoNodes(1));
  L2GlobalInt dummy(*gl2.pA,*gl2.pB);

  gl2.preAssemble(MNPC,this->getNoElms(true));
  return this->integrate(gl2,dummy,time) && gl2.solve(sField);
}


bool ASMbase::L2projection (Matrix& sField, FunctionBase* function, double t)
{
  PROFILE2("ASMbase::L2projection");

  GlbL2 gl2(function,this->getNoNodes(1));
  L2GlobalInt dummy(*gl2.pA,*gl2.pB);
  TimeDomain time; time.t = t;

  gl2.preAssemble(MNPC,this->getNoElms(true));
  return this->integrate(gl2,dummy,time) && gl2.solve(sField);
}


bool ASMbase::L2projection (const std::vector<Matrix*>& sField,
                            const FunctionVec& function, double t)
{
  PROFILE2("ASMbase::L2projection");

  GlbL2 gl2(function,this->getNoNodes(1));
  L2GlobalInt dummy(*gl2.pA,*gl2.pB);
  TimeDomain time; time.t = t;

  gl2.preAssemble(MNPC,this->getNoElms(true));
  return this->integrate(gl2,dummy,time) && gl2.solve(sField);
}


bool ASMbase::globalL2projection (Matrix& sField,
                                  const IntegrandBase& integrand,
                                  bool continuous) const
{
  if (this->empty()) return true; // silently ignore empty patches

  PROFILE2("ASMbase::globalL2");

  // Assemble the projection matrices
  size_t i, nnod = this->getNoProjectionNodes();
  size_t j, ncomp = integrand.getNoFields(2);
  SparseMatrix* A;
  StdVector* B;
  switch (GlbL2::MatrixType) {
  case LinAlg::UMFPACK:
    A = new SparseMatrix(SparseMatrix::UMFPACK);
    B = new StdVector(nnod*ncomp);
    break;
#ifdef HAS_PETSC
  case LinAlg::PETSC:
    if (GlbL2::SolverParams)
    {
      A = new PETScMatrix(ProcessAdm(), *GlbL2::SolverParams);
      B = new PETScVector(ProcessAdm(), nnod*ncomp);
      break;
    }
#endif
  default:
    A = new SparseMatrix(SparseMatrix::SUPERLU);
    B = new StdVector(nnod*ncomp);
  }
  A->redim(nnod,nnod);

  if (!this->assembleL2matrices(*A,*B,integrand,continuous))
  {
    delete A;
    delete B;
    return false;
  }

#if SP_DEBUG > 1
  std::cout <<"---- Matrix A -----\n"<< *A
            <<"-------------------"<< std::endl;
  std::cout <<"---- Vector B -----"<< *B
            <<"-------------------"<< std::endl;
#endif

  // Solve the patch-global equation system
  if (!A->solve(*B)) return false;

  // Store the control-point values of the projected field
  sField.resize(ncomp,nnod);
  for (i = 1; i <= nnod; i++)
    for (j = 1; j <= ncomp; j++)
      sField(j,i) = (*B)(i+(j-1)*nnod);

#if SP_DEBUG > 1
  std::cout <<"- Solution Vector -"<< sField
            <<"-------------------"<< std::endl;
#endif
  delete A;
  delete B;
  return true;
}
