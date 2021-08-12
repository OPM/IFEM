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
#include "ASMunstruct.h"
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
  \brief Expands a 2D tensor parametrization point to an unstructured one.
  \details Takes as input a tensor mesh, for instance
     in[0] = {0,1,2}
     in[1] = {2,3,5}
   and expands this to an unstructured representation, i.e.,
     out[0] = {0,1,2,0,1,2,0,1,2}
     out[1] = {2,2,2,3,3,3,5,5,5}
*/

static void expandTensorGrid2 (const RealArray* in, RealArray* out)
{
  out[0].resize(in[0].size()*in[1].size());
  out[1].resize(in[0].size()*in[1].size());

  size_t i, j, ip = 0;
  for (j = 0; j < in[1].size(); j++)
    for (i = 0; i < in[0].size(); i++, ip++) {
      out[0][ip] = in[0][i];
      out[1][ip] = in[1][j];
    }
}


/*!
  \brief Expands a 3D tensor parametrization point to an unstructured one.
  \details Takes as input a tensor mesh, for instance
     in[0] = {0,1}
     in[1] = {2,3}
     in[2] = {7,9}
   and expands this to an unstructured representation, i.e.,
     out[0] = {0,1,0,1,0,1,0,1}
     out[1] = {2,2,3,3,2,2,3,3}
     out[2] = {7,7,7,7,9,9,9,9}
*/

static void expandTensorGrid3 (const RealArray* in, RealArray* out)
{
  out[0].resize(in[0].size()*in[1].size()*in[2].size());
  out[1].resize(in[0].size()*in[1].size()*in[2].size());
  out[2].resize(in[0].size()*in[1].size()*in[2].size());

  size_t i, j, k, ip = 0;
  for (k = 0; k < in[2].size(); k++)
    for (j = 0; j < in[1].size(); j++)
      for (i = 0; i < in[0].size(); i++, ip++) {
        out[0][ip] = in[0][i];
        out[1][ip] = in[1][j];
        out[2][ip] = in[2][k];
      }
}


L2ProbIntegrand::L2ProbIntegrand (const ASMbase& patch, const IntegrandBase& itg) :
  L2Integrand(patch), m_itg(itg)
{
}


bool L2ProbIntegrand::evaluate (Matrix& sField, const RealArray* gpar) const
{
  return m_patch.evalSolution(sField, m_itg, gpar);
}


size_t L2ProbIntegrand::dim () const
{
  return m_itg.getNoFields(2);
}


L2FuncIntegrand::L2FuncIntegrand (const ASMbase& patch, const FunctionBase& func) :
  L2Integrand(patch), m_func(func)
{
}


bool L2FuncIntegrand::evaluate (Matrix& sField, const RealArray* gpar) const
{
  std::array<RealArray,3> expanded;
  if (dynamic_cast<const ASMunstruct*>(&m_patch) == nullptr) {
    if (m_patch.getNoSpaceDim() == 2) {
      expandTensorGrid2(gpar, expanded.data());
      gpar = expanded.data();
    } else if (m_patch.getNoSpaceDim() == 3) {
      expandTensorGrid3(gpar, expanded.data());
      gpar = expanded.data();
    }
  }

  sField.resize(m_func.dim(), gpar[0].size());

  Real2DMat u;
  m_patch.getParameterDomain(u);
  for (size_t i = 0; i < gpar[0].size(); ++i) {
    Vec3 X;
    std::vector<double> xi(3), param(3);
    xi[0] = u[0][0] + gpar[0][i] / (u[0][1] - u[0][0]);
    if (m_patch.getNoSpaceDim() > 1)
      xi[1] = u[1][0] + gpar[1][i] / (u[1][1] - u[1][0]);
    if (m_patch.getNoSpaceDim() > 2)
      xi[2] = u[2][0] + gpar[2][i] / (u[2][1] - u[2][0]);
    m_patch.evalPoint(xi.data(), param.data(), X);
    std::vector<Real> vals = m_func.getValue(X);
    for (size_t j = 1; j <= m_func.dim(); ++j)
      sField(j,i+1) = vals[j-1];
  }

  return true;
}


size_t L2FuncIntegrand::dim () const
{
  return m_func.dim();
}


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
    this->resize(p.pA ? 1 : 0,nf);
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
  L2GlobalInt(int dim_, SparseMatrix* A_, StdVector& B_)
    : dim(dim_), A(A_), B(B_) {}

  //! \brief Empty destructor.
  virtual ~L2GlobalInt() {}

  //! \brief Adds a LocalIntegral object into a corresponding global object.
  virtual bool assemble(const LocalIntegral* elmObj, int)
  {
    const L2Mats* elm = static_cast<const L2Mats*>(elmObj);
    for (size_t i = 0; i < elm->mnpc.size(); i++)
    {
      int inod = elm->mnpc[i]+1;
      if (A) {
        for (size_t j = 0; j < elm->mnpc.size(); j++)
        {
          int jnod = elm->mnpc[j]+1;
          (*A)(inod,jnod) += elm->A.front()(i+1,j+1);
        }
      }
      for (const Vector& b : elm->b)
      {
        B(inod) += b[i];
        inod += dim;
      }
    }
    return true;
  }

private:
  int dim; //!< Dimension of matrix
  SparseMatrix* A; //!< Reference to left-hand-side matrix
  StdVector&    B; //!< Reference to right-hand-side vector
};


GlbL2::GlbL2 (IntegrandBase* p, size_t n, SparseMatrix* A) :
  pA(A), problem(p)
{
  nrhs = p->getNoFields(2);
  this->allocate(n);
}


GlbL2::GlbL2 (FunctionBase* f, size_t n, SparseMatrix* A) :
  pA(A), problem(nullptr), functions({f})
{
  nrhs = f->dim();
  this->allocate(n);
}


GlbL2::GlbL2 (const FunctionVec& f, size_t n, SparseMatrix* A) :
  pA(A), problem(nullptr), functions(f)
{
  nrhs = 0;
  for (FunctionBase* func : f)
    nrhs += func->dim();
  this->allocate(n);
}


GlbL2::~GlbL2()
{
  delete pB;
}


void GlbL2::allocate (size_t n)
{
#ifdef HAS_PETSC
  if (GlbL2::MatrixType == LinAlg::PETSC && GlbL2::SolverParams)
  {
    static ProcessAdm adm;
    if (!pA) {
      pA = new PETScMatrix(adm,*GlbL2::SolverParams);
      pA->redim(n,n);
    }

    pB = new PETScVector(adm,n*nrhs);
  }
  else
#endif
  {
    if (!pA) {
      pA = new SparseMatrix(SparseMatrix::SUPERLU);
      pA->redim(n,n);
    }
    pB = new StdVector(n*nrhs);
  }
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

  if (!gl2.A.empty())
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

  if (!gl2.A.empty())
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
  GlbL2 gl2(integrand,this->getNoNodes(1),glbL2_A);
  L2GlobalInt dummy(this->getNoNodes(1),glbL2_A ? nullptr : gl2.pA,*gl2.pB);
  if (!glbL2_A)
    gl2.preAssemble(MNPC,this->getNoElms(true));
  glbL2_A = gl2.pA;
  return this->integrate(gl2,dummy,time) && gl2.solve(sField);
}


bool ASMbase::L2projection (Matrix& sField, FunctionBase* function, double t)
{
  PROFILE2("ASMbase::L2projection");
  GlbL2 gl2(function,this->getNoNodes(1),glbL2_A);
  L2GlobalInt dummy(this->getNoNodes(1),glbL2_A ? nullptr : gl2.pA,*gl2.pB);
  if (!glbL2_A)
    gl2.preAssemble(MNPC,this->getNoElms(true));
  glbL2_A = gl2.pA;
  TimeDomain time; time.t = t;
  return this->integrate(gl2,dummy,time) && gl2.solve(sField);
}


bool ASMbase::L2projection (const std::vector<Matrix*>& sField,
                            const FunctionVec& function, double t)
{
  PROFILE2("ASMbase::L2projection");
  GlbL2 gl2(function,this->getNoNodes(1),glbL2_A);
  L2GlobalInt dummy(this->getNoNodes(1),glbL2_A ? nullptr : gl2.pA,*gl2.pB);
  if (!glbL2_A)
    gl2.preAssemble(MNPC,this->getNoElms(true));
  glbL2_A = gl2.pA;
  TimeDomain time; time.t = t;
  return this->integrate(gl2,dummy,time) && gl2.solve(sField);
}


bool ASMbase::globalL2projection (Matrix& sField,
                                  const L2Integrand& integrand,
                                  bool continuous) const
{
  if (this->empty()) return true; // silently ignore empty patches

  PROFILE2("ASMbase::globalL2");

  // Assemble the projection matrices
  size_t i, nnod = this->getNoProjectionNodes();
  size_t j, ncomp = integrand.dim();
  StdVector* B;
  switch (GlbL2::MatrixType) {
  case LinAlg::UMFPACK:
    if (!glbL2_A)
      glbL2_A = new SparseMatrix(SparseMatrix::UMFPACK);
    B = new StdVector(nnod*ncomp);
    break;
#ifdef HAS_PETSC
  case LinAlg::PETSC:
    if (GlbL2::SolverParams)
    {
      static ProcessAdm adm;
      if (!glbL2_A)
        glbL2_A = new PETScMatrix(adm, *GlbL2::SolverParams);
      B = new PETScVector(adm, nnod*ncomp);
      break;
    }
#endif
  default:
    glbL2_A = new SparseMatrix(SparseMatrix::SUPERLU);
    B = new StdVector(nnod*ncomp);
  }
  glbL2_A->redim(nnod,nnod);

  if (!this->assembleL2matrices(*glbL2_A,*B,integrand,continuous))
  {
    delete B;
    return false;
  }

#if SP_DEBUG > 1
  std::cout <<"---- Matrix A -----\n"<< *glbL2_A
            <<"-------------------"<< std::endl;
  std::cout <<"---- Vector B -----"<< *B
            <<"-------------------"<< std::endl;
#endif

  // Solve the patch-global equation system
  if (!glbL2_A->solve(*B)) return false;

  // Store the control-point values of the projected field
  sField.resize(ncomp,nnod);
  for (i = 1; i <= nnod; i++)
    for (j = 1; j <= ncomp; j++)
      sField(j,i) = (*B)(i+(j-1)*nnod);

#if SP_DEBUG > 1
  std::cout <<"- Solution Vector -"<< sField
            <<"-------------------"<< std::endl;
#endif
  delete B;
  return true;
}
