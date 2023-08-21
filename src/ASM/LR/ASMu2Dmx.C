// $Id$
//==============================================================================
//!
//! \file ASMu2Dmx.C
//!
//! \date May 7 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for assembly of unstructured 2D spline mixed FE models.
//!
//==============================================================================

#include "ASMu2Dmx.h"

#include "GoTools/geometry/SplineSurface.h"

#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Element.h"
#include "LRSpline/Basisfunction.h"

#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "IFEM.h"
#include "LocalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "SplineUtils.h"
#include "Point.h"
#include "Profiler.h"
#include "Vec3.h"

#include <array>
#include <fstream>
#include <numeric>


ASMu2Dmx::ASMu2Dmx (unsigned char n_s, const CharVec& n_f)
  : ASMu2D(n_s, *std::max_element(n_f.begin(),n_f.end())), ASMmxBase(n_f)
{
  threadBasis = nullptr;
}


ASMu2Dmx::ASMu2Dmx (const ASMu2Dmx& patch, const CharVec& n_f)
  : ASMu2D(patch), ASMmxBase(n_f[0]==0?patch.nfx:n_f),
    m_basis(patch.m_basis)
{
  threadBasis = patch.threadBasis;
  nfx = patch.nfx;
  nb =  patch.nb;
}


const LR::LRSplineSurface* ASMu2Dmx::getBasis (int basis) const
{
  if (basis < 1 || basis > (int)m_basis.size())
    return nullptr;

  return m_basis[basis-1].get();
}


LR::LRSplineSurface* ASMu2Dmx::getBasis (int basis)
{
  if (basis < 1 || basis > (int)m_basis.size())
    return nullptr;

  return m_basis[basis-1].get();
}


bool ASMu2Dmx::read (std::istream& is, int basis)
{
  if (basis == 0)
    return this->ASMu2D::read(is);

  if (basis < 0 || basis > static_cast<int>(nfx.size()))
    return false;

  if (m_basis.empty()) {
    m_basis.resize(nfx.size());
    nb.resize(nfx.size(), 0);
  }

  m_basis[basis-1] = std::make_shared<LR::LRSplineSurface>();
  is >> *m_basis[basis-1];
  nb[basis-1] = m_basis[basis-1]->nBasisFunctions();
  m_basis[basis-1]->generateIDs();

  return true;
}


bool ASMu2Dmx::write (std::ostream& os, int basis) const
{
  if (basis == -1)
    os << *projBasis;
  else
    os << *m_basis[basis-1];

  return os.good();
}


void ASMu2Dmx::clear (bool retainGeometry)
{
  if (!retainGeometry) {
    // Erase the spline data
    for (SplinePtr& patch : m_basis)
      patch.reset();

    m_basis.clear();
  }

  // Erase the FE data
  this->ASMu2D::clear(retainGeometry);
}


size_t ASMu2Dmx::getNoNodes (int basis) const
{
  if (basis > (int)nb.size() || basis < 1)
    return this->ASMbase::getNoNodes(basis);

  return nb[basis-1];
}


unsigned char ASMu2Dmx::getNoFields (int basis) const
{
  if (basis > (int)m_basis.size() || basis < 0)
    basis = 0;

  if (basis == 0)
    return std::accumulate(nfx.begin(), nfx.end(), 0);

  return nfx[basis-1];
}


unsigned char ASMu2Dmx::getNodalDOFs (size_t inod) const
{
  if (this->isLMn(inod)) return nLag;
  size_t nbc=0;
  for (size_t i=0;i<nb.size();++i)
    if (inod <= (nbc+=nb[i]))
      return nfx[i];

  return nfx[0];
}


char ASMu2Dmx::getNodeType (size_t inod) const
{
  if (this->isLMn(inod))
    return 'L';
  size_t nbc=nb.front();
  if (inod <= nbc)
    return 'D';
  else for (size_t i = 1; i < nb.size(); i++)
    if (inod <= (nbc += nb[i]))
      return 'O'+i;

  return 'X';
}


void ASMu2Dmx::initMADOF (const int* sysMadof)
{
  this->initMx(MLGN,sysMadof);
}


void ASMu2Dmx::extractNodeVec (const RealArray& globRes, RealArray& nodeVec,
                               unsigned char, int basis) const
{
  this->extractNodeVecMx(globRes,nodeVec,basis);
}


bool ASMu2Dmx::injectNodeVec (const RealArray& nodeRes, RealArray& globRes,
                              unsigned char, int basis) const
{
  this->injectNodeVecMx(globRes,nodeRes,basis);
  return true;
}


bool ASMu2Dmx::getSolution (Matrix& sField, const Vector& locSol,
                            const IntVec& nodes) const
{
  return this->getSolutionMx(sField,locSol,nodes);
}


bool ASMu2Dmx::generateFEMTopology ()
{
  if (!myMLGN.empty())
    return true;

  auto createLR = [this](Go::SplineSurface& srf)
                  {
                    return srf.rational() ? this->createLRNurbs(srf)
                                          : new LR::LRSplineSurface(&srf);
                  };

  if (tensorPrjBas)
  {
    projBasis.reset(createLR(*tensorPrjBas));
    delete tensorPrjBas;
    tensorPrjBas = nullptr;
  }

  if (m_basis.empty()) {
    SurfaceVec svec = ASMmxBase::establishBases(tensorspline, ASMmxBase::Type);
    m_basis.resize(svec.size());
    for (size_t b = 0; b < svec.size(); b++)
      m_basis[b].reset(createLR(*svec[b]));

    // we need to project on something that is not one of our bases
    if (ASMmxBase::Type == ASMmxBase::REDUCED_CONT_RAISE_BASIS1 ||
        ASMmxBase::Type == ASMmxBase::REDUCED_CONT_RAISE_BASIS2 ||
        ASMmxBase::Type == ASMmxBase::DIV_COMPATIBLE ||
        ASMmxBase::Type == ASMmxBase::SUBGRID) {
      Go::SplineSurface* otherBasis = nullptr;
      if (!projBasis)
        otherBasis = ASMmxBase::raiseBasis(tensorspline);

      if (ASMmxBase::Type == ASMmxBase::SUBGRID) {
        refBasis.reset(createLR(*otherBasis));
        if (!projBasis)
          projBasis = m_basis.front();
        altProjBasis = refBasis;
      }
      else {
        if (!projBasis)
          projBasis.reset(createLR(*otherBasis));
        refBasis = projBasis;
      }
      delete otherBasis;
    }
    else {
      if (!projBasis)
        projBasis = m_basis[2-ASMmxBase::elmBasis];
      refBasis = projBasis;
    }

    is_rational = tensorspline->rational();
    delete tensorspline;
    tensorspline = nullptr;
  }
  projBasis->generateIDs();
  refBasis->generateIDs();
  lrspline = m_basis[elmBasis-1];

  nb.clear();
  nb.reserve(m_basis.size());
  for (const SplinePtr& it : m_basis) {
    nb.push_back(it->nBasisFunctions());
#ifdef SP_DEBUG
    std::cout <<"Basis "<< nb.size()
              <<":\nnumCoefs: "<< nb.back()
              <<"\norder: "<< it->order(0) <<" "<< it->order(1) << std::endl;
#endif
  }

  if (shareFE == 'F') return true;

  nel = m_basis[elmBasis-1]->nElements();
  nnod = std::accumulate(nb.begin(), nb.end(), 0);

  myMLGE.resize(nel,0);
  myMLGN.resize(nnod);
  myMNPC.resize(nel);
  for (SplinePtr& it : m_basis)
    it->generateIDs();

  size_t iel = 0;
  for (const LR::Element* el1 : m_basis[elmBasis-1]->getAllElements())
  {
    size_t nfunc = 0;
    for (const SplinePtr& it : m_basis)
      nfunc += it->getElement(it->getElementContaining(el1->midpoint()))->nBasisFunctions();
    myMNPC[iel].resize(nfunc);

    size_t lnod = 0, ofs = 0;
    for (const SplinePtr& it : m_basis) {
      const LR::Element* el2 = it->getElement(it->getElementContaining(el1->midpoint()));
      for (LR::Basisfunction* b : el2->support())
        myMNPC[iel][lnod++] = b->getId() + ofs;
      ofs += it->nBasisFunctions();
    }

    myMLGE[iel++] = ++gEl; // global element number over all patches
  }

  for (size_t inod = 0; inod < nnod; ++inod)
    myMLGN[inod] = ++gNod;

#ifdef SP_DEBUG
  std::cout <<"NEL = "<< nel <<" NNOD = "<< nnod << std::endl;
#endif

  geo = m_basis[elmBasis-1].get();
  this->generateBezierBasis();
  this->generateBezierExtraction();

  return true;
}


bool ASMu2Dmx::integrate (Integrand& integrand,
                          GlobalIntegral& glInt,
                          const TimeDomain& time)
{
  if (m_basis.empty())
    return true; // silently ignore empty patches

  PROFILE2("ASMu2Dmx::integrate(I)");

  if (integrand.getReducedIntegration(nGauss))
  {
    std::cerr <<" *** Reduced integration not available for mixed LR splines"
              << std::endl;
    return false;
  }

  bool use2ndDer = integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES;

  if (myCache.empty()) {
    myCache.emplace_back(std::make_unique<BasisFunctionCache>(*this, cachePolicy, 1));
    for (size_t b = 2; b <= this->getNoBasis(); ++b) {
      const BasisFunctionCache& c = static_cast<const BasisFunctionCache&>(*myCache.front());
      myCache.emplace_back(std::make_unique<BasisFunctionCache>(c,b));
    }
  }

  for (std::unique_ptr<ASMu2D::BasisFunctionCache>& cache : myCache) {
    cache->setIntegrand(&integrand);
    cache->init(use2ndDer ? 2 : 1);
  }

  ASMu2D::BasisFunctionCache& cache = *myCache.front();

  const std::array<int,2>& ng = cache.nGauss();
  const std::array<const double*,2>& xg = cache.coord();
  const std::array<const double*,2>& wg = cache.weight();

  ThreadGroups oneGroup;
  if (glInt.threadSafe()) oneGroup.oneGroup(nel);
  const IntMat& groups = glInt.threadSafe() ? oneGroup[0] : threadGroups[0];


  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t t = 0; t < groups.size() && ok; ++t)
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < groups[t].size(); ++e)
    {
      if (!ok)
        continue;

      std::vector<int>    els;
      std::vector<size_t> elem_sizes;
      this->getElementsAt(threadBasis->getElement(groups[t][e])->midpoint(),els,elem_sizes);

      MxFiniteElement fe(elem_sizes);
      Matrix   Xnod, Jac;
      Matrix3D Hess;
      double   dXidu[2];
      double   param[3] = { 0.0, 0.0, 0.0 };
      Vec4     X(param,time.t);

      int geoEl = els[elmBasis-1];
      fe.iel = MLGE[geoEl-1];

      // Get element area in the parameter space
      double dA = 0.25*this->getParametricArea(geoEl);
      if (dA < 0.0)
      {
        ok = false; // topology error (probably logic error)
        continue;
      }

      // Set up control point (nodal) coordinates for current element
      if (!this->getElementCoordinates(Xnod,geoEl))
      {
        ok = false;
        continue;
      }

      if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
        fe.h = this->getElementCorners(geoEl,fe.XC);

      if (integrand.getIntegrandType() & Integrand::G_MATRIX)
      {
        // Element size in parametric space
        dXidu[0] = geo->getElement(geoEl-1)->umax()-geo->getElement(geoEl-1)->umin();
        dXidu[1] = geo->getElement(geoEl-1)->vmax()-geo->getElement(geoEl-1)->vmin();
      }

      // Initialize element quantities
      LocalIntegral* A = integrand.getLocalIntegral(elem_sizes,fe.iel);
      if (!integrand.initElement(MNPC[geoEl-1],fe,elem_sizes,nb,*A))
      {
        A->destruct();
        ok = false;
        continue;
      }

      // --- Integration loop over all Gauss points in each direction ----------

      int jp = (geoEl-1)*ng[0]*ng[1];
      fe.iGP = firstIp + jp; // Global integration point counter

      size_t ig = 0;
      for (int j = 0; j < ng[1]; j++)
        for (int i = 0; i < ng[0]; i++, fe.iGP++, ig++)
        {
          // Local element coordinates of current integration point
          fe.xi  = xg[0][i];
          fe.eta = xg[1][j];

          // Parameter values of current integration point
          fe.u = param[0] = cache.getParam(0,geoEl-1,i);
          fe.v = param[1] = cache.getParam(1,geoEl-1,j);

          std::vector<const BasisFunctionVals*> bfs(this->getNoBasis());
          for (size_t b = 0; b < m_basis.size(); ++b) {
            bfs[b] = &myCache[b]->getVals(geoEl-1,ig);
            fe.basis(b+1) = bfs[b]->N;
          }

          // Compute Jacobian inverse of the coordinate mapping and
          // basis function derivatives w.r.t. Cartesian coordinates
          if (!fe.Jacobian(Jac,Xnod,elmBasis,&bfs))
            continue; // skip singular points

          // Compute Hessian of coordinate mapping and 2nd order derivatives
          if (use2ndDer && !fe.Hessian(Hess,Jac,Xnod,elmBasis,&bfs))
            ok = false;

          // Compute G-matrix
          if (integrand.getIntegrandType() & Integrand::G_MATRIX)
            utl::getGmat(Jac,dXidu,fe.G);

          // Cartesian coordinates of current integration point
          X.assign(Xnod * fe.basis(elmBasis));

          // Evaluate the integrand and accumulate element contributions
          fe.detJxW *= dA*wg[0][i]*wg[1][j];
          if (!integrand.evalIntMx(*A,fe,time,X))
            ok = false;
        }

      // Finalize the element quantities
      if (ok && !integrand.finalizeElement(*A,fe,time,firstIp+jp))
        ok = false;

      // Assembly of global system integral
      if (ok && !glInt.assemble(A->ref(),fe.iel))
        ok = false;

      A->destruct();
    }

  for (std::unique_ptr<ASMu2D::BasisFunctionCache>& cache : myCache)
    cache->finalizeAssembly();
  return ok;
}


bool ASMu2Dmx::integrate (Integrand& integrand, int lIndex,
                          GlobalIntegral& glInt,
                          const TimeDomain& time)
{
  if (m_basis.empty())
    return true; // silently ignore empty patches

  PROFILE2("ASMu2Dmx::integrate(B)");

  // Get Gaussian quadrature points and weights
  int nGP = integrand.getBouIntegrationPoints(nGauss);
  const double* xg = GaussQuadrature::getCoord(nGP);
  const double* wg = GaussQuadrature::getWeight(nGP);
  if (!xg || !wg) return false;

  // Find the parametric direction of the edge normal {-2,-1, 1, 2}
  const int edgeDir = (lIndex%10+1)/((lIndex%2) ? -2 : 2);

  const int t1 = abs(edgeDir);   // Tangent direction normal to the patch edge
  const int t2 = 3-abs(edgeDir); // Tangent direction along the patch edge

  std::array<Vector,2> gpar;
  for (int d = 0; d < 2; d++)
    if (-1-d == edgeDir)
    {
      gpar[d].resize(nGP);
      gpar[d].fill(d == 0 ? lrspline->startparam(0) : lrspline->startparam(1));
    }
    else if (1+d == edgeDir)
    {
      gpar[d].resize(nGP);
      gpar[d].fill(d == 0 ? lrspline->endparam(0) : lrspline->endparam(1));
    }

  // Extract the Neumann order flag (1 or higher) for the integrand
  integrand.setNeumannOrder(1 + lIndex/10);

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex%10);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;

  std::vector<Matrix> dNxdu(m_basis.size());
  Matrix Xnod, Jac;
  double param[3] = { 0.0, 0.0, 0.0 };
  Vec4   X(param,time.t);
  Vec3   normal;


  // === Assembly loop over all elements on the patch edge =====================

  int iel = 0;
  for (const LR::Element* el1 : m_basis[elmBasis-1]->getAllElements())
  {
    ++iel;

    // Skip elements that are not on current boundary edge
    bool skipMe = false;
    switch (edgeDir)
    {
      case -1: if (el1->umin() != m_basis[elmBasis-1]->startparam(0)) skipMe = true; break;
      case  1: if (el1->umax() != m_basis[elmBasis-1]->endparam(0)  ) skipMe = true; break;
      case -2: if (el1->vmin() != m_basis[elmBasis-1]->startparam(1)) skipMe = true; break;
      case  2: if (el1->vmax() != m_basis[elmBasis-1]->endparam(1)  ) skipMe = true; break;
    }
    if (skipMe) continue;

    if (!myElms.empty() && !glInt.threadSafe() &&
        std::find(myElms.begin(), myElms.end(), iel-1) == myElms.end())
      continue;

    std::vector<int>    els;
    std::vector<size_t> elem_sizes;
    this->getElementsAt(el1->midpoint(),els,elem_sizes);

    MxFiniteElement fe(elem_sizes,firstp);
    int geoEl = els[elmBasis-1];
    fe.iel = MLGE[geoEl-1];
    fe.xi = fe.eta = edgeDir < 0 ? -1.0 : 1.0;
    firstp += nGP; // Global integration point counter

    // Get element edge length in the parameter space
    double dS = 0.5*this->getParametricLength(geoEl,t2);
    if (dS < 0.0) return false; // topology error (probably logic error)

    // Set up control point coordinates for current element
    if (!this->getElementCoordinates(Xnod,geoEl))
      return false;

    if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
      fe.h = this->getElementCorners(iel,fe.XC);

    // Initialize element quantities
    LocalIntegral* A = integrand.getLocalIntegral(elem_sizes,fe.iel,true);
    bool ok = integrand.initElementBou(MNPC[geoEl-1],elem_sizes,nb,*A);

    // Get integration gauss points over this element
    this->getGaussPointParameters(gpar[t2-1],t2-1,nGP,geoEl,xg);


    // --- Integration loop over all Gauss points along the edge ---------------

    for (int i = 0; i < nGP && ok; i++, fe.iGP++)
    {
      // Local element coordinates and parameter values
      // of current integration point
      if (t1 == 2)
        fe.xi = xg[i];
      else
        fe.eta = xg[i];
      fe.u = param[0] = gpar[0][i];
      fe.v = param[1] = gpar[1][i];

      // Evaluate basis function derivatives at current integration points
      std::vector<Go::BasisDerivsSf> splinex(m_basis.size());
      for (size_t b=0; b < m_basis.size(); ++b) {
        this->computeBasis(fe.u, fe.v, splinex[b], els[b]-1, m_basis[b].get());
        SplineUtils::extractBasis(splinex[b],fe.basis(b+1),dNxdu[b]);
      }

      // Compute Jacobian inverse of the coordinate mapping and
      // basis function derivatives w.r.t. Cartesian coordinates
      fe.detJxW = utl::Jacobian(Jac,normal,fe.grad(elmBasis),Xnod,
                                dNxdu[elmBasis-1],t1,t2);
      if (fe.detJxW == 0.0) continue; // skip singular points

      for (size_t b = 1; b <= m_basis.size(); ++b)
        if ((int)b != elmBasis)
          fe.grad(b).multiply(dNxdu[b-1],Jac);

      if (edgeDir < 0)
        normal *= -1.0;

      // Cartesian coordinates of current integration point
      X.assign(Xnod * fe.basis(elmBasis));

      // Evaluate the integrand and accumulate element contributions
      fe.detJxW *= dS*wg[i];
      ok = integrand.evalBouMx(*A,fe,time,X,normal);
    }

    // Finalize the element quantities
    if (ok && !integrand.finalizeElementBou(*A,fe,time))
      ok = false;

    // Assembly of global system integral
    if (ok && !glInt.assemble(A->ref(),fe.iel))
      ok = false;

    A->destruct();

    if (!ok) return false;
  }

  return true;
}


bool ASMu2Dmx::integrate (Integrand& integrand,
                          GlobalIntegral& glInt,
                          const TimeDomain& time,
                          const ASM::InterfaceChecker& iChkgen)
{
  if (m_basis.empty())
    return true; // silently ignore empty patches

  if (!(integrand.getIntegrandType() & Integrand::INTERFACE_TERMS))
    return true; // No interface terms

  PROFILE2("ASMu2Dmx::integrate(J)");

  const ASMu2D::InterfaceChecker& iChk =
                          static_cast<const ASMu2D::InterfaceChecker&>(iChkgen);

  // Get Gaussian quadrature points and weights
  int nGP = integrand.getBouIntegrationPoints(nGauss);
  const double* xg = GaussQuadrature::getCoord(nGP);
  const double* wg = GaussQuadrature::getWeight(nGP);
  if (!xg || !wg) return false;

  Matrix   Xnod, Jac;
  Vec4     X(nullptr,time.t);
  Vec3     normal;

  int iel = 0;
  for (const LR::Element* el1 : m_basis.front()->getAllElements())
  {
    short int status = iChk.hasContribution(++iel);
    if (!status) continue; // no interface contributions for this element

    if (!myElms.empty() && !glInt.threadSafe() &&
        std::find(myElms.begin(), myElms.end(), iel-1) == myElms.end())
      continue;

    std::vector<int>    els;
    std::vector<size_t> elem_sizes;
    this->getElementsAt(el1->midpoint(),els,elem_sizes);
    // Replace first entry by hosting element
    els.front()        = iel;
    elem_sizes.front() = el1->nBasisFunctions();

    // Set up control point coordinates for current element
    if (!this->getElementCoordinates(Xnod,els[elmBasis-1]))
      return false;

    LocalIntegral* A = integrand.getLocalIntegral(elem_sizes, iel);
    bool ok = integrand.initElement(MNPC[els[elmBasis-1]-1],elem_sizes,nb,*A);
    size_t origSize = A->vec.size();

    int bit = 8;
    for (int iedge = 4; iedge > 0 && status > 0 && ok; iedge--, bit /= 2) {
      if (status & bit) {
        const int edgeDir = (iedge+1)/((iedge%2) ? -2 : 2);
        const int t1 = abs(edgeDir);   // Tangent direction normal to the edge
        const int t2 = 3-abs(edgeDir); // Tangent direction along the edge

        // Set up parameters
        double u1 = iedge != 2 ? el1->umin() : el1->umax();
        double v1 = iedge <  4 ? el1->vmin() : el1->vmax();
        double u2(u1);
        double v2(v1);
        const double epsilon = 1e-8;
        double epsu = 0.0, epsv = 0.0;
        if (iedge == 1)
          epsu = epsilon;
        if (iedge == 2)
          epsu = -epsilon;
        if (iedge == 3)
          epsv = epsilon;
        if (iedge == 4)
          epsv = -epsilon;

        std::vector<double> intersections = iChk.getIntersections(iel, iedge);
        for (size_t i = 0; i < intersections.size() && ok; ++i) {
          if (iedge == 1 || iedge == 2)
            v2 = intersections[i];
          else
            u2 = intersections[i];

          int el_neigh = this->getBasis(1)->getElementContaining({u1-epsu,v1-epsv})+1;
          const LR::Element* el2 = m_basis.front()->getElement(el_neigh-1);
          std::vector<int>    els2;
          std::vector<size_t> elem_sizes2;
          this->getElementsAt(el2->midpoint(),els2,elem_sizes2);
          els2.front()        = el_neigh;
          elem_sizes2.front() = el2->nBasisFunctions();

          LocalIntegral* A_neigh = integrand.getLocalIntegral(elem_sizes2, el_neigh);
          ok = integrand.initElement(MNPC[els2[elmBasis-1]-1],
                                     elem_sizes2,nb,*A_neigh);

          // Element sizes for both elements
          std::vector<size_t> elem_sizes3(elem_sizes);
          std::copy(elem_sizes2.begin(), elem_sizes2.end(),
                    std::back_inserter(elem_sizes3));

          MxFiniteElement fe(elem_sizes3);
          fe.h = this->getElementCorners(els2[elmBasis-1], fe.XC);

          if (!A_neigh->vec.empty()) {
            A->vec.resize(origSize+A_neigh->vec.size());
            std::copy(A_neigh->vec.begin(), A_neigh->vec.end(),
                      A->vec.begin()+origSize);
          }
          A_neigh->destruct();

          double dS = (iedge == 1 || iedge == 2) ? v2 - v1 : u2 - u1;

          std::array<Vector,2> gpar;
          if (iedge == 1 || iedge == 2) {
            gpar[0].resize(nGauss);
            gpar[0].fill(u1);
            gpar[1].resize(nGauss);
            for (int i = 0; i < nGauss; ++i)
              gpar[1][i] = 0.5*((v2-v1)*xg[i] + v2 + v1);
          } else {
            gpar[0].resize(nGauss);
            for (int i = 0; i < nGauss; ++i)
              gpar[0][i] = 0.5*((u2-u1)*xg[i] + u2 + u1);
            gpar[1].resize(nGauss);
            gpar[1].fill(v1);
          }
          Matrix Xnod2, Jac2;
          ok &= this->getElementCoordinates(Xnod2,els2[elmBasis-1]);

          for (int g = 0; g < nGP && ok; g++, ++fe.iGP)
          {
            // Local element coordinates and parameter values
            // of current integration point
            fe.xi = xg[g];
            fe.eta = xg[g];
            fe.u = gpar[0][g];
            fe.v = gpar[1][g];

            // Evaluate basis function derivatives at current integration points
            std::vector<Matrix> dNxdu(m_basis.size()*2);
            for (size_t b=0; b < m_basis.size(); ++b) {
              Go::BasisDerivsSf spline;
              this->computeBasis(fe.u+epsu, fe.v+epsv, spline, els[b]-1, m_basis[b].get());
              SplineUtils::extractBasis(spline,fe.basis(b+1),dNxdu[b]);
              this->computeBasis(fe.u-epsu, fe.v-epsv, spline, els2[b]-1, m_basis[b].get());
              SplineUtils::extractBasis(spline,fe.basis(b+1+m_basis.size()),
                                        dNxdu[b+m_basis.size()]);
            }

            // Compute Jacobian inverse of the coordinate mapping and
            // basis function derivatives w.r.t. Cartesian coordinates
            fe.detJxW = utl::Jacobian(Jac2,normal,
                                      fe.grad(elmBasis+m_basis.size()),Xnod2,
                                      dNxdu[elmBasis-1+m_basis.size()],t1,t2);
            fe.detJxW = utl::Jacobian(Jac,normal,
                                      fe.grad(elmBasis),Xnod,
                                      dNxdu[elmBasis-1],t1,t2);
            if (fe.detJxW == 0.0) continue; // skip singular points

            for (size_t b = 1; b <= m_basis.size(); ++b)
              if ((int)b != elmBasis) {
                fe.grad(b).multiply(dNxdu[b-1],Jac);
                fe.grad(b+m_basis.size()).multiply(dNxdu[b-1+m_basis.size()],Jac);
              }

            if (edgeDir < 0)
              normal *= -1.0;

            // Cartesian coordinates of current integration point
            X.assign(Xnod * fe.basis(elmBasis));

            // Evaluate the integrand and accumulate element contributions
            fe.detJxW *= 0.5*dS*wg[g];
            ok = integrand.evalIntMx(*A,fe,time,X,normal);
          }

          if (iedge == 1 || iedge == 2)
            v1 = v2;
          else
            u1 = u2;
        }
      }
    }

    // Finalize the element quantities
    if (ok && !integrand.finalizeElement(*A,FiniteElement(),time))
      ok = false;

    // Assembly of global system integral
    if (ok && !glInt.assemble(A,els[elmBasis-1]))
      ok = false;

    A->destruct();

    if (!ok) return false;
  }

  return true;
}


bool ASMu2Dmx::evalSolution (Matrix& sField, const Vector& locSol,
                             const RealArray* gpar, bool,
                             int deriv, int nf) const
{
  size_t nPoints = gpar[0].size();
  if (nPoints != gpar[1].size())
    return false;

  Vector   ptSol;
  std::vector<Matrix> dNxdu(m_basis.size()), dNxdX(m_basis.size());
  Matrix   Jac, Xnod, eSol, ptDer;

  std::vector<Go::BasisPtsSf> splinex(m_basis.size());

  std::vector<size_t> nc(nfx.size(), 0);
  if (nf)
    nc[0] = nf;
  else
    std::copy(nfx.begin(), nfx.end(), nc.begin());

  // Evaluate the primary solution field at each point
  sField.resize(std::accumulate(nc.begin(), nc.end(), 0), nPoints);
  for (size_t i = 0; i < nPoints; i++)
  {
    RealArray Ztmp;
    const double* locPtr = locSol.data();
    for (size_t j=0; j <  m_basis.size(); ++j) {
      if (nc[j] == 0)
        continue;

      // Fetch element containing evaluation point.
      // Sadly, points are not always ordered in the same way as the elements.
      int iel = m_basis[j]->getElementContaining(gpar[0][i],gpar[1][i]);
      const LR::Element* el = m_basis[j]->getElement(iel);

      // Evaluate basis function values/derivatives at current parametric point
      // and multiply with control point values to get the point-wise solution
      this->computeBasis(gpar[0][i],gpar[1][i],splinex[j],iel,m_basis[j].get());

      Matrix val1(nc[j], splinex[j].basisValues.size());
      size_t col = 0;
      for (LR::Basisfunction* b : el->support())
        val1.fillColumn(++col, locPtr+b->getId()*nc[j]);
      Vector Ytmp;
      val1.multiply(splinex[j].basisValues,Ytmp);
      Ztmp.insert(Ztmp.end(),Ytmp.begin(),Ytmp.end());
      locPtr += nb[j]*nc[j];
    }

    sField.fillColumn(i+1, Ztmp);
  }

  return true;
}


bool ASMu2Dmx::evalSolution (Matrix& sField, const IntegrandBase& integrand,
                             const RealArray* gpar, bool) const
{
#ifdef SP_DEBUG
  std::cout <<"ASMu2Dmx::evalSolution(Matrix&,const IntegrandBase&,const RealArray*,bool)\n";
#endif

  sField.resize(0,0);

  // TODO: investigate the possibility of doing "regular" refinement by
  //       uniform tesselation grid and ignoring LR mesh lines

  size_t nPoints = gpar[0].size();
  bool use2ndDer = integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES;
  if (nPoints != gpar[1].size())
    return false;

  // Evaluate the secondary solution field at each point
  for (size_t i = 0; i < nPoints; i++)
  {
    // Fetch element containing evaluation point
    // sadly, points are not always ordered in the same way as the elements
    std::vector<int>    els;
    std::vector<size_t> elem_sizes;
    this->getElementsAt({gpar[0][i],gpar[1][i]},els,elem_sizes);

    // Evaluate the basis functions at current parametric point
    MxFiniteElement       fe(elem_sizes,firstIp+i);
    std::vector<Matrix>   dNxdu(m_basis.size());
    std::vector<Matrix3D> d2Nxdu2(use2ndDer ? m_basis.size() : 0);
    Matrix Jac, Xnod;
    Matrix3D Hess;
    if (use2ndDer)
      for (size_t b = 0; b < m_basis.size(); ++b) {
        Go::BasisDerivsSf2 spline;
        this->computeBasis(gpar[0][i],gpar[1][i],spline,els[b]-1,m_basis[b].get());
        SplineUtils::extractBasis(spline,fe.basis(b+1),dNxdu[b],d2Nxdu2[b]);
      }
    else
      for (size_t b = 0; b < m_basis.size(); ++b) {
        Go::BasisDerivsSf spline;
        this->computeBasis(gpar[0][i],gpar[1][i],spline,els[b]-1,m_basis[b].get());
        SplineUtils::extractBasis(spline,fe.basis(b+1),dNxdu[b]);
      }

    // Set up control point (nodal) coordinates for current element
    if (!this->getElementCoordinates(Xnod,els[elmBasis-1])) return false;

    // Compute Jacobian inverse of the coordinate mapping and
    // basis function derivatives w.r.t. Cartesian coordinates
    if (!fe.Jacobian(Jac,Xnod,elmBasis,nullptr,&dNxdu))
      continue; // skip singular points

    // Compute Hessian of coordinate mapping and 2nd order derivatives
    if (use2ndDer && !fe.Hessian(Hess,Jac,Xnod,elmBasis,nullptr,&d2Nxdu2))
      return false;

    // Cartesian coordinates of current integration point
    fe.u = gpar[0][i];
    fe.v = gpar[1][i];
    utl::Point X4(Xnod*fe.basis(elmBasis),{fe.u,fe.v});

    // Now evaluate the solution field
    Vector solPt;
    if (!integrand.evalSol(solPt,fe,X4,
                           MNPC[els[elmBasis-1]-1],elem_sizes,nb))
      return false;
    else if (sField.empty())
      sField.resize(solPt.size(),nPoints,true);

    sField.fillColumn(1+i,solPt);
  }

  return true;
}


bool ASMu2Dmx::refine (const LR::RefineData& prm, Vectors& sol)
{
  if (shareFE)
    return true;

  if (prm.errors.empty() && prm.elements.empty())
    return true;

  for (Vector& solvec : sol)
    for (size_t j = 0; j < m_basis.size(); j++) {
      Vector bVec;
      this->extractNodeVec(solvec, bVec, 0, j+1);
      LR::extendControlPoints(m_basis[j].get(), bVec, nfx[j]);
    }

  if (doRefine(prm, refBasis.get())) {
    for (size_t j = 0; j < m_basis.size(); ++j)
      if (refBasis != m_basis[j]) {
        if ((j == 0 && ASMmxBase::Type == REDUCED_CONT_RAISE_BASIS1) ||
            (j == 1 && ASMmxBase::Type == REDUCED_CONT_RAISE_BASIS2))
          this->copyRefinement(m_basis[j].get(), 2);
        else
          this->copyRefinement(m_basis[j].get(), 1);
      }

    // Uniformly refine to find basis 1
    if (ASMmxBase::Type == ASMmxBase::SUBGRID) {
      m_basis[0].reset(refBasis->copy());
      projBasis = m_basis.front();
      size_t nFunc = refBasis->nBasisFunctions();
      IntVec elems(nFunc);
      std::iota(elems.begin(),elems.end(),0);
      m_basis[0]->refineBasisFunction(elems);
    }

    if (altProjBasis)
      altProjBasis->generateIDs();

    size_t len = 0;
    for (size_t j = 0; j< m_basis.size(); ++j) {
      m_basis[j]->generateIDs();
      nb[j] = m_basis[j]->nBasisFunctions();
      len += nfx[j]*nb[j];
    }

    size_t ofs = 0;
    for (int i = sol.size()-1; i >= 0; i--)
      for (size_t j = 0; j < m_basis.size(); ++j) {
        sol[i].resize(len);
        LR::contractControlPoints(m_basis[j].get(), sol[i], nfx[j], ofs);
        ofs += nfx[j]*nb[j];
      }

  #ifdef SP_DEBUG
    std::cout <<"Refined mesh: ";
    for (const SplinePtr& it : m_basis)
      std::cout << it->nElements() <<" ";
    std::cout <<"elements ";
    for (const SplinePtr& it : m_basis)
      std::cout << it->nBasisFunctions() <<" ";
    std::cout <<"nodes."<< std::endl;
    std::cout << "Projection basis: "
              << projBasis->nElements() << " elements "
              << projBasis->nBasisFunctions() << " nodes" << std::endl;
    std::cout << "Refinement basis: "
              << refBasis->nElements() << " elements "
              << refBasis->nBasisFunctions() << " nodes" << std::endl;
  #endif

    return true;
  }

  return false;
}


Vec3 ASMu2Dmx::getCoord (size_t inod) const
{
  size_t b = 0;
  size_t nbb = 0;
  while (b < nb.size() && nbb+nb[b] < inod)
    nbb += nb[b++];
  ++b;

  const LR::Basisfunction* basis = this->getBasis(b)->getBasisfunction(inod-nbb-1);
  if (!basis) {
    std::cerr << "Asked to get coordinate for node " << inod
              << ", but only have " << this->getBasis(b)->nBasisFunctions()
              << " nodes in basis " << b << std::endl;
    return Vec3();
  }
  return Vec3(&(*basis->cp()),nsd);
}


void ASMu2Dmx::generateThreadGroups (const Integrand& integrand, bool silence,
                                     bool ignoreGlobalLM)
{
  int p1 = 0;
  if (ASMmxBase::Type == ASMmxBase::DIV_COMPATIBLE)
    threadBasis = this->getBasis(3);
  else
    for (size_t i = 1; i <= m_basis.size(); ++i)
      if (this->getBasis(i)->order(0) > p1) {
        threadBasis = this->getBasis(i);
        p1 = threadBasis->order(0);
      }

  std::vector<LR::LRSpline*> secConstraint;
  if (ASMmxBase::Type == ASMmxBase::SUBGRID ||
      ASMmxBase::Type == REDUCED_CONT_RAISE_BASIS1)
    secConstraint = {this->getBasis(2)};
  if (ASMmxBase::Type == REDUCED_CONT_RAISE_BASIS2)
    secConstraint = {this->getBasis(1)};
  if (ASMmxBase::Type == ASMmxBase::DIV_COMPATIBLE)
    secConstraint = {this->getBasis(1), this->getBasis(2)};

  LR::generateThreadGroups(threadGroups,threadBasis,secConstraint);
  LR::generateThreadGroups(projThreadGroups,projBasis.get());
  if (altProjBasis)
    LR::generateThreadGroups(altProjThreadGroups,altProjBasis.get());

  std::vector<const LR::LRSpline*> bases;
  for (const SplinePtr& basis : m_basis)
    bases.push_back(basis.get());

  if (silence || threadGroups[0].size() < 2) return;

  this->checkThreadGroups(threadGroups[0], bases, threadBasis);

  IFEM::cout <<"\nMultiple threads are utilized during element assembly.";
#ifdef SP_DEBUG
  for (size_t i = 0; i < threadGroups[0].size(); i++)
    IFEM::cout <<"\n Color "<< i+1 <<": "
               << threadGroups[0][i].size() <<" elements";
  IFEM::cout << std::endl;
#else
  this->analyzeThreadGroups(threadGroups[0]);
#endif
}


bool ASMu2Dmx::connectPatch (int edge, ASM2D& neighbor, int nedge, bool revers,
                             int basis, bool coordCheck, int thick)
{
  ASMu2Dmx* neighMx = dynamic_cast<ASMu2Dmx*>(&neighbor);
  if (!neighMx) return false;

  for (size_t i = 1; i <= m_basis.size(); ++i)
    if (basis == 0 || i == (size_t)basis)
      if (!this->connectBasis(edge,*neighMx,nedge,revers,i,0,0,coordCheck,thick))
        return false;

  this->addNeighbor(neighMx);
  return true;
}


void ASMu2Dmx::getBoundaryNodes (int lIndex, IntVec& nodes, int basis,
                                 int thick, int orient, bool local) const
{
  if (basis > 0)
    this->ASMu2D::getBoundaryNodes(lIndex, nodes, basis, thick, orient, local);
  else
    for (size_t b = 1; b <= this->getNoBasis(); ++b)
      this->ASMu2D::getBoundaryNodes(lIndex, nodes, b, thick, orient, local);
}


void ASMu2Dmx::remapErrors (RealArray& errors,
                            const RealArray& origErr, bool elemErrors) const
{
  const LR::LRSplineSurface* geo = this->getBasis(ASMmxBase::elmBasis);
  for (const LR::Element* elm : geo->getAllElements()) {
    int rEl = refBasis->getElementContaining((elm->umin()+elm->umax())/2.0,
                                             (elm->vmin()+elm->vmax())/2.0);
    if (elemErrors)
      errors[rEl] += origErr[elm->getId()];
    else
      for (LR::Basisfunction* b : refBasis->getElement(rEl)->support())
        errors[b->getId()] += origErr[elm->getId()];
  }
}


size_t ASMu2Dmx::getNoRefineNodes() const
{
  return refBasis->nBasisFunctions();
}


size_t ASMu2Dmx::getNoRefineElms() const
{
  return refBasis->nElements();
}


void ASMu2Dmx::storeMesh (const std::string& fName, int fType) const
{
  auto&& writeBasis = [fName,fType,this](std::shared_ptr<LR::LRSplineSurface> patch,
                                         const std::string& tag)
  {
    std::string fileName = "_patch_" + tag + "_" + fName + ".eps";
    if (fType%2) {
      std::ofstream meshFile("param"+fileName);
      patch->writePostscriptMesh(meshFile);
    }
    if ((fType/2)%2) {
      std::ofstream meshFile("physical"+fileName);
      if (is_rational)
        const_cast<ASMu2Dmx*>(this)->writePostscriptElementsNurbs(patch, meshFile);
      else
        patch->writePostscriptElements(meshFile);
    }
    if ((fType/4)%2) {
      std::ofstream meshFile("param_dot"+fileName);
      patch->writePostscriptElements(meshFile);
    }
    if ((fType/8)%2) {
      std::ofstream meshFile("physical_dot"+fileName);
      if (is_rational)
        const_cast<ASMu2Dmx*>(this)->writePostscriptMeshWithControlPointsNurbs(patch, meshFile);
      else
        patch->writePostscriptMeshWithControlPoints(meshFile);
    }
  };

  std::string btag("basis1");
  for (const SplinePtr& patch : m_basis)
  {
    writeBasis(patch, btag);
    ++btag.back();
  }
  writeBasis(projBasis, "proj");
  writeBasis(refBasis, "ref");
}


void ASMu2Dmx::copyRefinement (LR::LRSplineSurface* basis,
                               int multiplicity) const
{
  for (const LR::Meshline* line : refBasis->getAllMeshlines()) {
    int mult = line->multiplicity_ > 1 ? line->multiplicity_ : multiplicity;
    if (line->span_u_line_)
      basis->insert_const_v_edge(line->const_par_,
                                 line->start_, line->stop_, mult);
    else
      basis->insert_const_u_edge(line->const_par_,
                                 line->start_, line->stop_, mult);
  }
}


void ASMu2Dmx::swapProjectionBasis ()
{
  if (altProjBasis) {
    ASMmxBase::elmBasis = ASMmxBase::elmBasis == 1 ? 2 : 1;
    std::swap(projBasis, altProjBasis);
    std::swap(projThreadGroups, altProjThreadGroups);
    lrspline = m_basis[ASMmxBase::elmBasis-1];
    geo = lrspline.get();
    this->generateBezierBasis();
    this->generateBezierExtraction();
  }
}


void ASMu2Dmx::getElementsAt (const RealArray& param,
                              std::vector<int>& elms,
                              std::vector<size_t>& sizes) const
{
  elms.clear();
  sizes.clear();
  elms.reserve(m_basis.size());
  sizes.reserve(m_basis.size());
  for (const SplinePtr& basis : m_basis)
  {
    int iel = basis->getElementContaining(param);
    elms.push_back(1+iel);
    sizes.push_back(basis->getElement(iel)->nBasisFunctions());
  }
}


BasisFunctionVals ASMu2Dmx::BasisFunctionCache::calculatePt (size_t el,
                                                             size_t gp,
                                                             bool reduced) const
{
  const std::array<size_t,2> gpIdx = this->gpIndex(gp,reduced);
  double u = this->getParam(0,el,gpIdx[0],reduced);
  double v = this->getParam(1,el,gpIdx[1],reduced);

  const ASMu2Dmx& pch = static_cast<const ASMu2Dmx&>(patch);

  const LR::Element* el1 = pch.getBasis(ASMmxBase::elmBasis)->getElement(el);
  size_t el_b = patch.getBasis(basis)->getElementContaining(el1->midpoint());

  BasisFunctionVals result;
  if (nderiv == 1 || reduced) {
    Go::BasisDerivsSf spline;
    pch.computeBasis(u,v,spline,el_b,patch.getBasis(basis));
    SplineUtils::extractBasis(spline,result.N,result.dNdu);
  } else if (nderiv == 2) {
    Go::BasisDerivsSf2 spline;
    pch.computeBasis(u,v,spline,el_b,patch.getBasis(basis));
    SplineUtils::extractBasis(spline,result.N,result.dNdu,result.d2Ndu2);
  }

  return result;
}
