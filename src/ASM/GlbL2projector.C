// $Id$
//==============================================================================
//!
//! \file GlbL2projector.C
//!
//! \date Oct 16 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief General integrands for global L2-projection.
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
#include <numeric>


namespace GlbL2
{
  //! Matrix type for global L2-projection (version 1 only)
  LinAlg::MatrixType MatrixType = LinAlg::SPARSE;
  //! Linear solver parameters for the global L2-projection (version 1 only)
  LinSolParams* SolverParams = nullptr;
}


namespace
{
  /*!
    \brief Expands a 2D tensor parametrization point to an unstructured one.
    \details Takes as input a tensor mesh, for instance
    \code
      in[0] = {0,1,2};
      in[1] = {2,3,5};
    \endcode
    and expands this to an unstructured representation, i.e.,
    \code
      out[0] = {0,1,2,0,1,2,0,1,2};
      out[1] = {2,2,2,3,3,3,5,5,5};
    \endcode
  */

  const RealArray* expandTensorGrid2 (const RealArray* in, RealArray* out)
  {
    out[0].resize(in[0].size()*in[1].size());
    out[1].resize(in[0].size()*in[1].size());

    size_t i, j, ip = 0;
    for (j = 0; j < in[1].size(); j++)
      for (i = 0; i < in[0].size(); i++, ip++) {
        out[0][ip] = in[0][i];
        out[1][ip] = in[1][j];
      }

    return out;
  }


  /*!
    \brief Expands a 3D tensor parametrization point to an unstructured one.
    \details Takes as input a tensor mesh, for instance
    \code
      in[0] = {0,1}
      in[1] = {2,3}
      in[2] = {7,9}
    \endcode
    and expands this to an unstructured representation, i.e.,
    \code
      out[0] = {0,1,0,1,0,1,0,1}
      out[1] = {2,2,3,3,2,2,3,3}
      out[2] = {7,7,7,7,9,9,9,9}
    \endcode
  */

  const RealArray* expandTensorGrid3 (const RealArray* in, RealArray* out)
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

    return out;
  }
}


L2ProbIntegrand::L2ProbIntegrand (const ASMbase& patch,
                                  const IntegrandBase& itg,
                                  const ProcessAdm& adm) :
  L2Integrand(patch, adm), m_itg(itg)
{
  myTime = itg.getTimeLevel();
}


bool L2ProbIntegrand::evaluate (Matrix& sField, const RealArray* gpar) const
{
  return m_patch.evalSolution(sField, m_itg, gpar);
}


size_t L2ProbIntegrand::dim () const
{
  return m_itg.getNoFields(2);
}


L2FuncIntegrand::L2FuncIntegrand (const ASMbase& patch,
                                  const FunctionBase& func,
                                  const ProcessAdm& adm, double time) :
  L2Integrand(patch, adm), m_func(func)
{
  myTime = time;
}


bool L2FuncIntegrand::evaluate (Matrix& sField, const RealArray* gpar) const
{
  std::array<RealArray,3> expanded;
  if (!dynamic_cast<const ASMunstruct*>(&m_patch))
  {
    if (m_patch.getNoSpaceDim() == 2)
      gpar = expandTensorGrid2(gpar, expanded.data());
    else if (m_patch.getNoSpaceDim() == 3)
      gpar = expandTensorGrid3(gpar, expanded.data());
  }

  sField.resize(m_func.dim(), gpar[0].size());

  Real2DMat u;
  m_patch.getParameterDomain(u);

  double xi[3];    // Normalized spline domain parameters always in [0,1]
  double param[3]; // Actual domain parameters to be used for evaluation
  Vec4   X(param,myTime);
  for (size_t i = 0; i < gpar[0].size(); ++i) {
    for (size_t j = 0; j < 3 && j < m_patch.getNoSpaceDim(); ++j)
      xi[j] = u[j][0] + gpar[j][i] / (u[j][1] - u[j][0]);
    m_patch.evalPoint(xi, param, X);
    sField.fillColumn(i+1, m_func.getValue(X));
  }

  return true;
}


size_t L2FuncIntegrand::dim () const
{
  return m_func.dim();
}


namespace
{
  using uIntVec     = std::vector<size_t>;
  using FunctionVec = std::vector<FunctionBase*>;

  class GL2;


  /*!
    \brief Local integral container class for L2-projections (version 2).
  */

  class L2Mats : public ElmMats
  {
  public:
    //! \brief The constructor initializes pointers and references.
    //! \param[in] p Global L2 integrand object containing projection matrices
    //! \param[in] iEl Global element number (1-based)
    //! \param[in] nen Number of element nodes
    //! \param[in] nf Number of field components
    //! \param[in] q The main problem integrand
    L2Mats(GL2& p, size_t iEl, size_t nen, size_t nf, Integrand* q)
      : gl2Int(p)
    {
      elmData = q ? q->getLocalIntegral(nen,iEl) : nullptr;
      this->resize(1,nf);
      this->redim(nen);
    }

    //! \brief Destruction method to clean up after numerical integration.
    void destruct() override { delete elmData; delete this; }

    GL2&           gl2Int;  //!< The global L2-projection integrand
    LocalIntegral* elmData; //!< Element data associated with problem integrand

    IntVec  mnpc;        //!< Matrix of element nodal correspondance
    uIntVec elem_sizes;  //!< Size of each basis on the element
    uIntVec basis_sizes; //!< Size of each basis on the patch
  };


  /*!
    \brief General integrand for L2-projection of secondary solutions.
    \details This class is used only for version 2 of the L2-projection.
  */

  class GL2 : public Integrand
  {
  public:
    //! \brief The constructor initializes the projection matrices.
    //! \param[in] prb The main problem integrand
    //! \param[in] n Dimension of the L2-projection matrices (number of nodes)
    GL2(IntegrandBase* prb, size_t n) : A(SparseMatrix::SUPERLU)
    {
      problem = prb;
      nrhs = prb->getNoFields(2);
      A.redim(n,n);
      B.redim(n*nrhs);
    }

    //! \brief Alternative constructor for projection of an explicit function.
    //! \param[in] func The function to do L2-projection on
    //! \param[in] n Dimension of the L2-projection matrices (number of nodes)
    GL2(FunctionBase* func, size_t n) : A(SparseMatrix::SUPERLU)
    {
      problem = nullptr;
      functions = { func };
      nrhs = func->dim();
      A.redim(n,n);
      B.redim(n*nrhs);
    }

    //! \brief Alternative constructor for projection of explicit functions.
    //! \param[in] funcs The functions to do L2-projection on
    //! \param[in] n Dimension of the L2-projection matrices (number of nodes)
    GL2(const FunctionVec& funcs, size_t n) : A(SparseMatrix::SUPERLU)
    {
      problem = nullptr;
      functions = funcs;
      nrhs = std::accumulate(funcs.begin(), funcs.end(), 0u,
                             [](size_t a, const FunctionBase* f)
                             { return a + f->dim(); });
      A.redim(n,n);
      B.redim(n*nrhs);
    }

    //! \brief Returns current solution mode.
    SIM::SolutionMode getMode(bool) const override { return SIM::RECOVERY; }

    //! \brief Defines which FE quantities are needed by the integrand.
    int getIntegrandType() const override
    {
      if (problem) // Mask off the element interface flag
        return problem->getIntegrandType() & ~INTERFACE_TERMS;
      else
        return STANDARD;
    }

    using Integrand::getLocalIntegral;
    //! \brief Returns a local integral object for the given element.
    //! \param[in] nen Number of nodes on element
    //! \param[in] iEl Global element number (1-based)
    LocalIntegral* getLocalIntegral(size_t nen, size_t iEl, bool) const override
    {
      return new L2Mats(*const_cast<GL2*>(this),iEl,nen,nrhs,problem);
    }

    //! \brief Initializes current element for numerical integration.
    //! \param[in] MNPC Matrix of nodal point correspondance for current element
    //! \param[in] fe Nodal and integration point data for current element
    //! \param[in] X0 Cartesian coordinates of the element center
    //! \param[in] nPt Number of integration points in this element
    //! \param elmInt Local integral for element
    bool initElement(const IntVec& MNPC, const FiniteElement& fe,
                     const Vec3& X0, size_t nPt,
                     LocalIntegral& elmInt) override
    {
      L2Mats& gl2 = static_cast<L2Mats&>(elmInt);

      gl2.mnpc = MNPC;
      if (problem && gl2.elmData)
        return problem->initElement(MNPC,fe,X0,nPt,*gl2.elmData);
      else
        return true;
    }

    //! \brief Initializes current element for numerical integration (mixed).
    //! \param[in] MNPC1 Matrix of nodal point correspondance for the element
    //! \param[in] fe Nodal and integration point data for current element
    //! \param[in] elem_sizes Size of each basis on the element
    //! \param[in] basis_sizes Size of each basis on the patch
    //! \param elmInt Local integral for element
    bool initElement(const IntVec& MNPC1, const MxFiniteElement& fe,
                     const uIntVec& elem_sizes, const uIntVec& basis_sizes,
                     LocalIntegral& elmInt) override
    {
      L2Mats& gl2 = static_cast<L2Mats&>(elmInt);

      gl2.mnpc = MNPC1;
      gl2.elem_sizes = elem_sizes;
      gl2.basis_sizes = basis_sizes;
      if (problem && gl2.elmData)
        return problem->initElement(MNPC1,fe,elem_sizes,basis_sizes,
                                    *gl2.elmData);
      else
        return true;
    }

    //! \brief Dummy implementation.
    bool initElement(const IntVec&, LocalIntegral&) override { return false; }
    //! \brief Dummy implementation.
    bool initElement(const IntVec&, const uIntVec&, const uIntVec&,
                     LocalIntegral&) override { return false; }

    //! \brief Dummy implementation.
    bool initElementBou(const IntVec&,
                        LocalIntegral&) override { return false; }
    //! \brief Dummy implementation.
    bool initElementBou(const IntVec&, const uIntVec&, const uIntVec&,
                        LocalIntegral&) override { return false; }

    using Integrand::evalInt;
    //! \brief Evaluates the integrand at an interior point.
    //! \param elmInt The local integral object to receive the contributions
    //! \param[in] fe Finite element data of current integration point
    //! \param[in] X Cartesian coordinates of current integration point
    bool evalInt(LocalIntegral& elmInt,
                 const FiniteElement& fe, const Vec3& X) const override
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
        solPt.push_back(funcPt.begin(),funcPt.end());
      }

      gl2.A.front().outer_product(fe.N,fe.N,true,fe.detJxW);
      for (size_t j = 0; j < solPt.size(); j++)
        gl2.b[j].add(fe.N,solPt[j]*fe.detJxW);

      return true;
    }

    using Integrand::evalIntMx;
    //! \brief Evaluates the integrand at an interior point.
    //! \param elmInt The local integral object to receive the contributions
    //! \param[in] fe Mixed finite element data of current integration point
    //! \param[in] X Cartesian coordinates of current integration point
    bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                   const Vec3& X) const override
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

    //! \brief Pre-computes the sparsity pattern of the projection matrix \b A.
    //! \param[in] MMNPC Matrix of matrices of nodal point correspondances
    //! \param[in] nel Number of elements
    //! \param[in] checkNZ If \e true, flag non-zero contributions in assembly
    void preAssemble(const std::vector<IntVec>& MMNPC, size_t nel, bool checkNZ)
    {
      A.preAssemble(MMNPC,nel);
      if (checkNZ)
        A.initNonZeroEqs();
    }

    //! \brief Solves the projection equation system and evaluates nodal values.
    //! \param[out] sField Nodal/control-point values of the projected results.
    bool solve(Matrix& sField)
    {
      // Insert a 1.0 value on the diagonal for equations with no contributions.
      // Needed when immersed boundaries with "totally outside" elements.
      const size_t nnod = A.dim();
      for (size_t i = 1; i <= nnod; i++)
        if (A(i,i) == 0.0) A(i,i) = 1.0;

#if SP_DEBUG > 1
      std::cout <<"\nGlobal L2-projection matrix:\n"<< A;
      std::cout <<"\nGlobal L2-projection RHS:"<< B;
#endif

      // Solve the patch-global equation system
      if (!A.solve(B,sField)) return false;

#if SP_DEBUG > 1
      std::cout <<"\nSolution:"<< sField;
#endif
      return true;
    }

    //! \brief Solves the projection equation system and evaluates nodal values.
    //! \param[out] sField Nodal/control-point values of the projected results.
    bool solve(const std::vector<Matrix*>& sField)
    {
      if (sField.size() != functions.size())
      {
        std::cerr <<" *** GL2::solve: Logic error, size(field)="<< sField.size()
                  <<" != size(functions)="<< functions.size() << std::endl;
        return false;
      }

      // Insert a 1.0 value on the diagonal for equations with no contributions.
      // Needed when immersed boundaries with "totally outside" elements.
      const size_t nnod = A.dim();
      for (size_t i = 1; i <= nnod; i++)
        if (A(i,i) == 0.0) A(i,i) = 1.0;

#if SP_DEBUG > 1
      std::cout <<"\nGlobal L2-projection matrix:\n"<< A;
      std::cout <<"\nGlobal L2-projection RHS:"<< B;
#endif

      // Solve the patch-global equation system
      Matrix& solutions = *sField.front();
      if (!A.solve(B,solutions)) return false;

      if (sField.size() > 1)
      {
        // Split the solution vectors into one for each function
        size_t offset = functions.front()->dim();
        for (size_t k = 1; k < sField.size(); k++)
        {
          size_t ncomp = functions[k]->dim();
          sField[k]->resize(ncomp,nnod);
          for (size_t i = 1; i <= nnod; i++)
            for (size_t j = 1; j <= ncomp; j++)
              (*sField[k])(j,i) = solutions(offset+j,i);
          offset += ncomp;
        }
        sField.front()->expandRows(functions.front()->dim()-solutions.rows());
      }
#if SP_DEBUG > 1
      for (size_t k = 0; k < sField.size(); k++)
        std::cout <<"\nSolution "<< k <<":"<< *sField[k];
#endif
      return true;
    }

  private:
    SparseMatrix A; //!< Left-hand-side matrix of the L2-projection
    StdVector    B; //!< Right-hand-side vectors of the L2-projection

    IntegrandBase* problem; //!< The main problem integrand
    FunctionVec  functions; //!< Explicit functions to L2-project
    size_t            nrhs; //!< Number of right-hand-size vectors

    friend class L2GlobalInt;
  };


  /*!
    \brief Global integral container class for L2-projections (version 2).
  */

  class L2GlobalInt : public GlobalIntegral
  {
  public:
    //! \brief The constructor initializes the system matrix references.
    L2GlobalInt(GL2& gl2) : A(gl2.A), B(gl2.B) {}

    //! \brief Adds a LocalIntegral object into a corresponding global object.
    bool assemble(const LocalIntegral* elmObj, int) override
    {
      const L2Mats* elm = static_cast<const L2Mats*>(elmObj);
      for (size_t i = 0; i < elm->mnpc.size(); i++)
      {
        int inod = elm->mnpc[i]+1;
        A.flagNonZeroEq(inod);
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
}


/*!
  This method implements version 2 of the global L2-projection
  for secondary solution variables of the provided \a integrand.
*/

bool ASMbase::L2projection (Matrix& sField,
                            IntegrandBase* integrand,
                            const TimeDomain& time)
{
  PROFILE2("ASMbase::L2projection");

  GL2 gl2(integrand,this->getNoNodes(1));
  L2GlobalInt dummy(gl2);

  gl2.preAssemble(MNPC,this->getNoElms(true),this->getElementActivator());
  return this->integrate(gl2,dummy,time) && gl2.solve(sField);
}


/*!
  This method implements version 2 of the global L2-projection
  of an explicit function.
*/

bool ASMbase::L2projection (Matrix& sField, FunctionBase* function, double t)
{
  PROFILE2("ASMbase::L2projection");

  GL2 gl2(function,this->getNoNodes(1));
  L2GlobalInt dummy(gl2);

  gl2.preAssemble(MNPC,this->getNoElms(true),this->getElementActivator());
  return this->integrate(gl2,dummy,TimeDomain(t)) && gl2.solve(sField);
}


/*!
  This method implements version 2 of the global L2-projection
  of a set of explicit functions.
*/

bool ASMbase::L2projection (const std::vector<Matrix*>& sField,
                            const FunctionVec& function, double t)
{
  PROFILE2("ASMbase::L2projection");

  GL2 gl2(function,this->getNoNodes(1));
  L2GlobalInt dummy(gl2);

  gl2.preAssemble(MNPC,this->getNoElms(true),this->getElementActivator());
  return this->integrate(gl2,dummy,TimeDomain(t)) && gl2.solve(sField);
}


/*!
  This method implements version 1 of the global L2-projection,
  and supports both discrete and continuous projections.
*/

bool ASMbase::globalL2projection (Matrix& sField,
                                  const L2Integrand& integrand,
                                  bool continuous, bool enforceEnds) const
{
  if (this->empty()) return true; // silently ignore empty patches

  PROFILE2("ASMbase::globalL2");

  const size_t npnod = this->getNoProjectionNodes();
  const size_t ncomp = integrand.dim();
#ifdef HAS_PETSC
  PETScMatrix* Ap = nullptr;
  ProcessAdm singleAdm;
#endif

  // Assemble the projection matrices
  SystemMatrix* A;
  SystemVector* B;
  switch (GlbL2::MatrixType) {
  case LinAlg::UMFPACK:
    A = new SparseMatrix(npnod,npnod,SparseMatrix::UMFPACK);
    B = new StdVector(npnod*ncomp);
    break;

#ifdef HAS_PETSC
  case LinAlg::PETSC:
    if (GlbL2::SolverParams)
    {
      const IntMat nodes = this->getElmNodes(ASM::PROJECTION_BASIS);
      const bool useAdm = neighbors.empty() && nodes.size() == nel;
      const bool isPart = integrand.getAdm().dd.isPartitioned();
      const ProcessAdm& adm = isPart ? integrand.getAdm() : singleAdm;

      IntMat neighs;
      if (adm.dd.isPartitioned() && !useAdm && adm.getProcId() == 0) {
        neighs.resize(nodes.size());
        this->getElmConnectivities(neighs, ASM::PROJECTION_BASIS);
      }

      A = Ap = new PETScMatrix(adm, *GlbL2::SolverParams);
      Ap->init(npnod, &nodes, &neighs, useAdm ? &adm.dd.getElms() : nullptr);
      if (adm.dd.isPartitioned())
        const_cast<ASMbase*>(this)->generateProjThreadGroupsFromElms(Ap->getDD().getElms());
      B = new PETScVectors(*Ap, ncomp);
      break;
    }
#endif

  default:
    A = new SparseMatrix(npnod,npnod,SparseMatrix::SUPERLU);
    B = new StdVector(npnod*ncomp);
  }

  if (this->getElementActivator())
    A->initNonZeroEqs();

  if (!this->assembleL2matrices(*A,*B,integrand,continuous))
  {
    delete A;
    delete B;
    return false;
  }

  A->endAssembly();
  B->endAssembly();

#if SP_DEBUG > 1
  std::cout <<"---- Matrix A -----\n"<< *A
            <<"-------------------"<< std::endl;
  const Vector& bVec = B->vec();
  if (!bVec.empty())
  {
    std::cout <<"---- Vector B -----";
    for (size_t i = 1; i <= npnod; i++)
      for (size_t j = 0; j < ncomp; j++)
	std::cout << (j > 0 ? " " : "\n") << bVec(i+npnod*j);
    std::cout <<"\n-------------------"<< std::endl;
  }
#endif

  // Solve the patch-global equation system
  if (!A->solve(*B,sField))
  {
    delete A;
    delete B;
    sField.clear();
    return false;
  }

#if SP_DEBUG > 1
  std::cout <<"- Solution Vector -"<< sField
            <<"-------------------"<< std::endl;
#endif
  delete A;
  delete B;
  if (!enforceEnds) return true;

  // Get parameter and node numbers for the domain corners
  Real2DMat u;
  IntVec corners;
  if (!this->getParameterDomain(u,&corners))
    return true; // Silently ignore if corner points are not provided

  // Evaluate the solution at the corners
  Matrix sCorner;
  if (!integrand.evaluate(sCorner,u.data()))
    return false;

  // Enforce the corner values in the projected field
  for (size_t i = 0; i < corners.size(); i++)
  {
#if SP_DEBUG > 1
    std::cout <<"Replacing end/corner-point values of projected field at node "
              << corners[i] <<"\nfrom"<< Vector(sField.getColumn(corners[i]))
              <<"to"<< Vector(sCorner.getColumn(1+i));
#endif
    sField.fillColumn(corners[i],sCorner.getColumn(1+i));
  }

  return true;
}
