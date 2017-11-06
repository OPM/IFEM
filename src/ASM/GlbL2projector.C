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
#include "IntegrandBase.h"
#include "Function.h"
#include "Profiler.h"
#ifdef HAS_PETSC
#include "PETScMatrix.h"
#include "LinSolParams.h"
#include "ProcessAdm.h"
#endif


SystemMatrix::Type GlbL2::MatrixType = SystemMatrix::SPARSE;
LinSolParams* GlbL2::SolverParams = nullptr;


/*!
  \brief Local integral container class for L2-projections.
*/

class L2Mats : public LocalIntegral
{
public:
  //! \brief The constructor initializes pointers and references.
  //! \param[in] p The global L2 integrand object containing projection matrices
  //! \param[in] q Pointer to element data associated with the problem integrand
  L2Mats(GlbL2& p, LocalIntegral* q = nullptr) : gl2Int(p), elmData(q) {}
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


GlbL2::GlbL2 (IntegrandBase* p, size_t n) : A(SparseMatrix::SUPERLU)
{
  problem = p;
  function = nullptr;

  A.redim(n,n);
  B.redim(n*p->getNoFields(2));
}


GlbL2::GlbL2 (FunctionBase* f, size_t n) : A(SparseMatrix::SUPERLU)
{
  problem = nullptr;
  function = f;

  A.redim(n,n);
  B.redim(n*f->dim());
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
    return new L2Mats(*const_cast<GlbL2*>(this),
                      problem->getLocalIntegral(nen,iEl,neumann));
  else
    return new L2Mats(*const_cast<GlbL2*>(this));
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
  if (problem)
  {
    if (!problem->evalSol(solPt,fe,X,gl2.mnpc))
      if (!problem->diverged(fe.iGP+1))
        return false;
  }
  else if (function)
    solPt = function->getValue(X);

  return this->formL2Mats(gl2.mnpc,solPt,fe,X);
}


bool GlbL2::evalIntMx (LocalIntegral& elmInt,
                       const MxFiniteElement& fe,
                       const Vec3& X) const

{
  L2Mats& gl2 = static_cast<L2Mats&>(elmInt);

  Vector solPt;
  if (problem)
  {
    if (!problem->evalSol(solPt,fe,X,gl2.mnpc,gl2.elem_sizes,gl2.basis_sizes))
      if (!problem->diverged(fe.iGP+1))
        return false;
  }
  else if (function)
    solPt = function->getValue(X);

  return this->formL2Mats(gl2.mnpc,solPt,fe,X);
}


bool GlbL2::formL2Mats (const IntVec& mnpc, const Vector& solPt,
                        const FiniteElement& fe, const Vec3& X) const
{
  size_t a, b, nnod = A.dim();
  for (a = 0; a < fe.N.size(); a++)
  {
    int inod = mnpc[a]+1;
    for (b = 0; b < fe.N.size(); b++)
    {
      int jnod = mnpc[b]+1;
      A(inod,jnod) += fe.N[a]*fe.N[b]*fe.detJxW;
    }
    for (b = 0; b < solPt.size(); b++)
      B(inod+b*nnod) += fe.N[a]*solPt[b]*fe.detJxW;
  }

  return true;
}


void GlbL2::preAssemble (const std::vector<IntVec>& MMNPC, size_t nel)
{
  A.preAssemble(MMNPC,nel);
}


bool GlbL2::solve (Matrix& sField)
{
  // Insert a 1.0 value on the diagonal for equations with no contributions.
  // Needed in immersed boundary calculations with "totally outside" elements.
  size_t i, nnod = A.dim();
  for (i = 1; i <= nnod; i++)
    if (A(i,i) == 0.0) A(i,i) = 1.0;

#if SP_DEBUG > 1
  std::cout <<"\nGlobal L2-projection matrix:\n"<< A;
  std::cout <<"\nGlobal L2-projection RHS:"<< B;
#endif

  // Solve the patch-global equation system
  if (!A.solve(B)) return false;

  // Store the nodal values of the projected field
  size_t j, ncomp = 0;
  if (problem)
    ncomp = problem->getNoFields(2);
  else if (function)
    ncomp = function->dim();
  sField.resize(ncomp,nnod);
  for (i = 1; i <= nnod; i++)
    for (j = 1; j <= ncomp; j++)
      sField(j,i) = B(i+(j-1)*nnod);

#if SP_DEBUG > 1
  std::cout <<"\nSolution:"<< sField;
#endif
  return true;
}


bool ASMbase::L2projection (Matrix& sField,
                            IntegrandBase* integrand,
                            const TimeDomain& time)
{
  PROFILE2("ASMbase::L2projection");

  GlbL2 gl2(integrand,this->getNoNodes(1));
  GlobalIntegral dummy;

  gl2.preAssemble(MNPC,this->getNoElms(true));
  return this->integrate(gl2,dummy,time) && gl2.solve(sField);
}


bool ASMbase::L2projection (Matrix& sField, FunctionBase* function, double t)
{
  PROFILE2("ASMbase::L2projection");

  GlbL2 gl2(function,this->getNoNodes(1));
  GlobalIntegral dummy;
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
#ifdef HAS_PETSC
  if (GlbL2::MatrixType == SystemMatrix::PETSC && GlbL2::SolverParams)
  {
    A = new PETScMatrix(ProcessAdm(), *GlbL2::SolverParams, LinAlg::SYMMETRIC);
    B = new PETScVector(ProcessAdm(), nnod*ncomp);
  }
  else
#endif
  {
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
