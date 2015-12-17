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
#include "Profiler.h"


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

  GlbL2&         gl2Int;       //!< The global L2 projection integrand
  LocalIntegral* elmData;      //!< Element data associated with problem integrand
  IntVec         mnpc;         //!< Matrix of element nodal correspondance for bases
  std::vector<size_t> elem_sizes;   //!< Size of each basis on the element
  std::vector<size_t> basis_sizes;   //!< Size of each basis on the patch
};


GlbL2::GlbL2 (IntegrandBase& p, size_t n) : problem(p), A(SparseMatrix::SUPERLU)
{
  A.redim(n,n);
  B.redim(n*p.getNoFields(2));
}


int GlbL2::getIntegrandType () const
{
  // Mask off the element interface flag
  return problem.getIntegrandType() & ~INTERFACE_TERMS;
}


LocalIntegral* GlbL2::getLocalIntegral (size_t nen, size_t iEl,
                                        bool neumann) const
{
  return new L2Mats(*const_cast<GlbL2*>(this),
		    problem.getLocalIntegral(nen,iEl,neumann));
}


bool GlbL2::initElement (const IntVec& MNPC, const FiniteElement& fe,
                         const Vec3& Xc, size_t nPt,
                         LocalIntegral& elmInt)
{
  L2Mats& gl2 = static_cast<L2Mats&>(elmInt);

  gl2.mnpc = MNPC;
  return problem.initElement(MNPC,fe,Xc,nPt,*gl2.elmData);
}


bool GlbL2::initElement (const IntVec& MNPC1,
                         const std::vector<size_t>& elem_sizes,
                         const std::vector<size_t>& basis_sizes,
                         LocalIntegral& elmInt)
{
  L2Mats& gl2 = static_cast<L2Mats&>(elmInt);

  gl2.mnpc  = MNPC1;
  gl2.elem_sizes = elem_sizes;
  gl2.basis_sizes = basis_sizes;
  return problem.initElement(MNPC1,elem_sizes,basis_sizes,*gl2.elmData);
}


bool GlbL2::evalInt (LocalIntegral& elmInt,
                     const FiniteElement& fe,
                     const Vec3& X) const

{
  L2Mats& gl2 = static_cast<L2Mats&>(elmInt);

  Vector solPt;
  if (!problem.evalSol(solPt,fe,X,gl2.mnpc))
    if (!problem.diverged(fe.iGP+1))
      return false;

  size_t a, b, nnod = A.dim();
  for (a = 0; a < fe.N.size(); a++)
  {
    int inod = gl2.mnpc[a]+1;
    for (b = 0; b < fe.N.size(); b++)
    {
      int jnod = gl2.mnpc[b]+1;
      A(inod,jnod) += fe.N[a]*fe.N[b]*fe.detJxW;
    }
    for (b = 0; b < solPt.size(); b++)
      B(inod+b*nnod) += fe.N[a]*solPt[b]*fe.detJxW;
  }

  return true;
}


bool GlbL2::evalIntMx (LocalIntegral& elmInt,
                       const MxFiniteElement& fe,
                       const Vec3& X) const

{
  L2Mats& gl2 = static_cast<L2Mats&>(elmInt);

  Vector solPt;
  if (!problem.evalSol(solPt,fe,X,gl2.mnpc,gl2.elem_sizes,gl2.basis_sizes))
    if (!problem.diverged(fe.iGP+1))
      return false;

  size_t a, b, nnod = A.dim();
  for (a = 0; a < fe.basis(1).size(); a++)
  {
    int inod = gl2.mnpc[a]+1;
    for (b = 0; b < fe.basis(1).size(); b++)
    {
      int jnod = gl2.mnpc[b]+1;
      A(inod,jnod) += fe.basis(1)[a]*fe.basis(1)[b]*fe.detJxW;
    }
    for (b = 0; b < solPt.size(); b++)
      B(inod+b*nnod) += fe.basis(1)[a]*solPt[b]*fe.detJxW;
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
  size_t nnod = A.dim();
  for (size_t j = 1; j <= nnod; j++)
    if (A(j,j) == 0.0) A(j,j) = 1.0;

  // Solve the patch-global equation system
  if (!A.solve(B)) return false;

  // Store the nodal values of the projected field
  size_t ncomp = B.dim() / nnod;
  sField.resize(ncomp,nnod);
  for (size_t i = 1; i <= nnod; i++)
    for (size_t j = 1; j <= ncomp; j++)
      sField(j,i) = B(i+(j-1)*nnod);

  return true;
}


bool ASMbase::L2projection (Matrix& sField,
			    const IntegrandBase& integrand,
			    const TimeDomain& time)
{
  PROFILE2("ASMbase::L2projection");

  GlbL2 gl2(const_cast<IntegrandBase&>(integrand),this->getNoNodes(1));
  GlobalIntegral dummy;

  gl2.preAssemble(MNPC,this->getNoElms(true));
  return this->integrate(gl2,dummy,time) && gl2.solve(sField);
}
