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

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SurfaceInterpolator.h"

#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Element.h"
#include "LRSpline/Basisfunction.h"

#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Profiler.h"
#include "Vec3Oper.h"
#include "Vec3.h"
#include "LRSplineFields2D.h"

#include <array>
#include <fstream>
#include <numeric>


ASMu2Dmx::ASMu2Dmx (unsigned char n_s, const CharVec& n_f)
  : ASMu2D(n_s, *std::max_element(n_f.begin(),n_f.end())), ASMmxBase(n_f)
{
  threadBasis = nullptr;
}


ASMu2Dmx::ASMu2Dmx (const ASMu2Dmx& patch, const CharVec& n_f)
  : ASMu2D(patch), ASMmxBase(n_f[0]==0?patch.nfx:n_f)
{
  m_basis = patch.m_basis;
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


bool ASMu2Dmx::write (std::ostream& os, int basis) const
{
  os << *m_basis[basis-1];
  return os.good();
}


void ASMu2Dmx::clear (bool retainGeometry)
{
  if (!retainGeometry) {
    // Erase the spline data
    for (auto& patch : m_basis)
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


void ASMu2Dmx::extractNodeVec (const Vector& globRes, Vector& nodeVec,
                               unsigned char, int basis) const
{
  this->extractNodeVecMx(globRes,nodeVec,basis);
}


bool ASMu2Dmx::injectNodeVec (const Vector& nodeRes, Vector& globRes,
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

  if (m_basis.empty()) {
    auto vec = ASMmxBase::establishBases(tensorspline, ASMmxBase::Type);
    m_basis.resize(vec.size());
    for (size_t i=0;i<vec.size();++i)
      m_basis[i].reset(new LR::LRSplineSurface(vec[i].get()));

    // we need to project on something that is not one of our bases
    if (ASMmxBase::Type == ASMmxBase::REDUCED_CONT_RAISE_BASIS1 ||
        ASMmxBase::Type == ASMmxBase::DIV_COMPATIBLE ||
        ASMmxBase::Type == ASMmxBase::SUBGRID) {
      std::shared_ptr<Go::SplineSurface> otherBasis =
          ASMmxBase::establishBases(tensorspline, ASMmxBase::FULL_CONT_RAISE_BASIS1).front();
      if (ASMmxBase::Type == ASMmxBase::SUBGRID) {
        projBasis = m_basis[0];
        refBasis.reset(new LR::LRSplineSurface(otherBasis.get()));
      } else {
        projBasis.reset(new LR::LRSplineSurface(otherBasis.get()));
        refBasis = projBasis;
      }
    } else
     projBasis = refBasis = m_basis[0];
  }
  projBasis->generateIDs();
  refBasis->generateIDs();
  lrspline = m_basis[geoBasis-1];

  nb.resize(m_basis.size());
  for (size_t i=0; i < m_basis.size(); ++i)
    nb[i] = m_basis[i]->nBasisFunctions();

  if (shareFE == 'F') return true;

#ifdef SP_DEBUG
  size_t nbasis=0;
  for (auto& it : m_basis) {
    std::cout << "Basis " << ++nbasis << ":\n";
    std::cout <<"numCoefs: "<< it->nBasisFunctions();
    std::cout <<"\norder: "<< it->order(0) <<" "<< it->order(1) << std::endl;
  }
#endif

  nel = m_basis[geoBasis-1]->nElements();

  nnod = std::accumulate(nb.begin(), nb.end(), 0);

  myMLGE.resize(nel,0);
  myMLGN.resize(nnod);
  myMNPC.resize(nel);
  for (auto&& it : m_basis)
    it->generateIDs();

  std::vector<LR::Element*>::iterator el_it1 = m_basis[geoBasis-1]->elementBegin();
  for (size_t iel=0; iel<nel; iel++, ++el_it1)
  {
    double uh = ((*el_it1)->umin()+(*el_it1)->umax())/2.0;
    double vh = ((*el_it1)->vmin()+(*el_it1)->vmax())/2.0;
    size_t nfunc = 0;
    for (size_t i=0; i<m_basis.size();++i) {
      auto el_it2 = m_basis[i]->elementBegin() +
                    m_basis[i]->getElementContaining(uh, vh);
      nfunc += (*el_it2)->nBasisFunctions();
    }
    myMLGE[iel] = ++gEl; // global element number over all patches
    myMNPC[iel].resize(nfunc);

    int lnod = 0;
    size_t ofs=0;
    for (size_t i=0; i<m_basis.size();++i) {
      auto el_it2 = m_basis[i]->elementBegin() +
                    m_basis[i]->getElementContaining(uh, vh);
      for (LR::Basisfunction *b : (*el_it2)->support())
        myMNPC[iel][lnod++] = b->getId()+ofs;
      ofs += nb[i];
    }
  }

  for (size_t inod = 0; inod < nnod; ++inod)
    myMLGN[inod] = ++gNod;

#ifdef SP_DEBUG
  std::cout <<"NEL = "<< nel <<" NNOD = "<< nnod << std::endl;
#endif

  geo = m_basis[geoBasis-1].get();

  return true;
}


bool ASMu2Dmx::integrate (Integrand& integrand,
                          GlobalIntegral& glInt,
                          const TimeDomain& time)
{
  if (m_basis.empty())
    return true; // silently ignore empty patches

  PROFILE2("ASMu2Dmx::integrate(I)");

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;
  bool use2ndDer = integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES;

  // === Assembly loop over all elements in the patch ==========================
  bool ok = true;
  for (size_t t = 0; t < threadGroups[0].size() && ok; ++t)
  {
//#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < threadGroups[0][t].size(); ++e)
    {
      if (!ok)
        continue;
      int iel = threadGroups[0][t][e] + 1;
      auto el1 = threadBasis->elementBegin()+iel-1;
      double uh = ((*el1)->umin()+(*el1)->umax())/2.0;
      double vh = ((*el1)->vmin()+(*el1)->vmax())/2.0;
      std::vector<size_t> els;
      std::vector<size_t> elem_sizes;
      for (size_t i=0; i < m_basis.size(); ++i) {
        els.push_back(m_basis[i]->getElementContaining(uh, vh)+1);
        elem_sizes.push_back(m_basis[i]->getElement(els.back()-1)->nBasisFunctions());
      }

      int geoEl = els[geoBasis-1];

      MxFiniteElement fe(elem_sizes);
      fe.iel = MLGE[geoEl-1];
      std::vector<Matrix> dNxdu(m_basis.size());
      Matrix   Xnod, Jac;
      double   param[3] = { 0.0, 0.0, 0.0 };
      Vec4     X(param);
      std::vector<Matrix3D> d2Nxdu2(m_basis.size());
      Matrix3D Hess;
      double   dXidu[2];

      // Get element area in the parameter space
      double dA = this->getParametricArea(geoEl);
      if (dA < 0.0)  // topology error (probably logic error)
      {
        ok = false;
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

      // Compute parameter values of the Gauss points over this element
      std::array<RealArray,2> gpar;
      for (int d = 0; d < 2; d++)
        this->getGaussPointParameters(gpar[d],d,nGauss,geoEl,xg);

      // Initialize element quantities
      LocalIntegral* A = integrand.getLocalIntegral(elem_sizes,fe.iel,false);
      if (!integrand.initElement(MNPC[geoEl-1], elem_sizes, nb, *A))
      {
        A->destruct();
        ok = false;
        continue;
      }

      // --- Integration loop over all Gauss points in each direction ------------

      int jp = (iel-1)*nGauss*nGauss;
      fe.iGP = firstIp + jp; // Global integration point counter

      for (int j = 0; j < nGauss; j++)
        for (int i = 0; i < nGauss; i++, fe.iGP++)
        {
          // Local element coordinates of current integration point
          fe.xi  = xg[i];
          fe.eta = xg[j];

          // Parameter values of current integration point
          fe.u = param[0] = gpar[0][i];
          fe.v = param[1] = gpar[1][j];

          // Compute basis function derivatives at current integration point
          if (use2ndDer)
            for (size_t b = 0; b < m_basis.size(); ++b) {
              Go::BasisDerivsSf2 spline;
              m_basis[b]->computeBasis(fe.u,fe.v,spline,els[b]-1);
              SplineUtils::extractBasis(spline,fe.basis(b+1),dNxdu[b],d2Nxdu2[b]);
            }
          else
            for (size_t b=0; b < m_basis.size(); ++b) {
              Go::BasisDerivsSf spline;
              m_basis[b]->computeBasis(fe.u, fe.v, spline, els[b]-1);
              SplineUtils::extractBasis(spline,fe.basis(b+1),dNxdu[b]);
            }

          // Compute Jacobian inverse of coordinate mapping and derivatives
          // basis function derivatives w.r.t. Cartesian coordinates
          fe.detJxW = utl::Jacobian(Jac,fe.grad(geoBasis),Xnod,dNxdu[geoBasis-1]);
          if (fe.detJxW == 0.0) continue; // skip singular points
          for (size_t b = 0; b < m_basis.size(); ++b)
            if (b != (size_t)geoBasis-1)
              fe.grad(b+1).multiply(dNxdu[b],Jac);

          // Compute Hessian of coordinate mapping and 2nd order derivatives
          if (use2ndDer) {
            if (!utl::Hessian(Hess,fe.hess(geoBasis),Jac,Xnod,
                              d2Nxdu2[geoBasis-1],fe.grad(geoBasis),true))
              ok = false;

            for (size_t b = 0; b < m_basis.size() && ok; ++b)
              if ((int)b != geoBasis)
                utl::Hessian(Hess,fe.hess(b+1),Jac,Xnod,
                             d2Nxdu2[b],fe.grad(b+1),false);
          }

          // Compute G-matrix
          if (integrand.getIntegrandType() & Integrand::G_MATRIX)
            utl::getGmat(Jac,dXidu,fe.G);

          // Cartesian coordinates of current integration point
          X.assign(Xnod * fe.basis(geoBasis));
          X.t = time.t;

          // Evaluate the integrand and accumulate element contributions
          fe.detJxW *= 0.25*dA*wg[i]*wg[j];
          if (!integrand.evalIntMx(*A,fe,time,X))
          {
            ok = false;
            continue;
          }
        }

      // Finalize the element quantities
      if (!integrand.finalizeElement(*A,time,firstIp+jp))
      {
        ok = false;
        continue;
      }

      // Assembly of global system integral
      if (!glInt.assemble(A->ref(),fe.iel))
      {
        ok = false;
        continue;
      }

      A->destruct();
    }
  }

  return ok;
}


bool ASMu2Dmx::integrate (Integrand& integrand, int lIndex,
                          GlobalIntegral& glInt,
                          const TimeDomain& time)
{
  if (!m_basis[0] || !m_basis[1])
    return true; // silently ignore empty patches

  PROFILE2("ASMu2Dmx::integrate(B)");

  // Get Gaussian quadrature points and weights
  int nGP = integrand.getBouIntegrationPoints(nGauss);
  const double* xg = GaussQuadrature::getCoord(nGP);
  const double* wg = GaussQuadrature::getWeight(nGP);
  if (!xg || !wg) return false;

  // Find the parametric direction of the edge normal {-2,-1, 1, 2}
  const int edgeDir = (lIndex%10+1)/(lIndex%2 ? -2 : 2);

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
  double   param[3] = { 0.0, 0.0, 0.0 };
  Vec4   X(param);
  Vec3   normal;

  // === Assembly loop over all elements on the patch edge =====================

  std::vector<LR::Element*>::iterator el1 = m_basis[geoBasis-1]->elementBegin();
  for (int iel = 1; el1 != m_basis[geoBasis-1]->elementEnd(); ++el1, ++iel)
  {
    // Skip elements that are not on current boundary edge
    bool skipMe = false;
    switch (edgeDir)
    {
      case -1: if ((*el1)->umin() != m_basis[geoBasis-1]->startparam(0)) skipMe = true; break;
      case  1: if ((*el1)->umax() != m_basis[geoBasis-1]->endparam(0)  ) skipMe = true; break;
      case -2: if ((*el1)->vmin() != m_basis[geoBasis-1]->startparam(1)) skipMe = true; break;
      case  2: if ((*el1)->vmax() != m_basis[geoBasis-1]->endparam(1)  ) skipMe = true; break;
    }
    if (skipMe) continue;

    double uh = ((*el1)->umin()+(*el1)->umax())/2.0;
    double vh = ((*el1)->vmin()+(*el1)->vmax())/2.0;
    std::vector<size_t> els;
    std::vector<size_t> elem_sizes;
    for (size_t i=0; i < m_basis.size(); ++i) {
      els.push_back(m_basis[i]->getElementContaining(uh, vh)+1);
      elem_sizes.push_back((*(m_basis[i]->elementBegin()+(els.back()-1)))->nBasisFunctions());
    }
    int geoEl = els[geoBasis-1];

    // Get element edge length in the parameter space
    double dS = this->getParametricLength(geoEl,t1);
    if (dS < 0.0) return false; // topology error (probably logic error)

    // Set up control point coordinates for current element
    if (!this->getElementCoordinates(Xnod,geoEl))
      return false;

    // Initialize element quantities
    MxFiniteElement fe(elem_sizes);
    fe.iel = MLGE[geoEl-1];
    fe.xi = fe.eta = edgeDir < 0 ? -1.0 : 1.0;
    LocalIntegral* A = integrand.getLocalIntegral(elem_sizes,fe.iel,true);
    if (!integrand.initElementBou(MNPC[geoEl-1], elem_sizes, nb, *A))
    {
      A->destruct();
      return false;
    }

    if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
      fe.h = this->getElementCorners(iel,fe.XC);

    // Get integration gauss points over this element
    this->getGaussPointParameters(gpar[t2-1],t2-1,nGP,geoEl,xg);

    // --- Integration loop over all Gauss points along the edge -------------

    fe.iGP = firstp; // Global integration point counter
    firstp += nGP;

    for (int i = 0; i < nGP; i++, ++fe.iGP)
    {
      // Local element coordinates and parameter values
      // of current integration point
      fe.xi = xg[i];
      fe.eta = xg[i];
      fe.u = param[0] = gpar[0][i];
      fe.v = param[1] = gpar[1][i];

      // Evaluate basis function derivatives at current integration points
      std::vector<Go::BasisDerivsSf> splinex(m_basis.size());
      for (size_t b=0; b < m_basis.size(); ++b) {
        m_basis[b]->computeBasis(fe.u, fe.v, splinex[b], els[b]-1);
        SplineUtils::extractBasis(splinex[b],fe.basis(b+1),dNxdu[b]);
      }

      // Compute Jacobian inverse of the coordinate mapping and
      // basis function derivatives w.r.t. Cartesian coordinates
      fe.detJxW = utl::Jacobian(Jac,normal,fe.grad(geoBasis),Xnod,dNxdu[geoBasis-1],t1,t2);
      if (fe.detJxW == 0.0) continue; // skip singular points
      for (size_t b = 0; b < m_basis.size(); ++b)
        if (b != (size_t)geoBasis-1)
          fe.grad(b+1).multiply(dNxdu[b],Jac);

      if (edgeDir < 0)
        normal *= -1.0;

      // Cartesian coordinates of current integration point
      X.assign(Xnod * fe.basis(geoBasis));
      X.t = time.t;

      // Evaluate the integrand and accumulate element contributions
      fe.detJxW *= 0.5*dS*wg[i];
      if (!integrand.evalBouMx(*A,fe,time,X,normal))
        return false;
    }

    // Finalize the element quantities
    if (!integrand.finalizeElementBou(*A,fe,time))
      return false;

    // Assembly of global system integral
    if (!glInt.assemble(A,fe.iel))
      return false;

    A->destruct();
  }

  return true;
}


bool ASMu2Dmx::integrate (Integrand& integrand,
                          GlobalIntegral& glInt,
                          const TimeDomain& time,
                          const ASM::InterfaceChecker& iChkgen)
{
  if (!geo) return true; // silently ignore empty patches
  if (!(integrand.getIntegrandType() & Integrand::INTERFACE_TERMS)) return true;

  PROFILE2("ASMu2Dmx::integrate(J)");

  const ASMu2D::InterfaceChecker& iChk =
                          static_cast<const ASMu2D::InterfaceChecker&>(iChkgen);

  // Get Gaussian quadrature points and weights
  int nGP = integrand.getBouIntegrationPoints(nGauss);
  const double* xg = GaussQuadrature::getCoord(nGP);
  const double* wg = GaussQuadrature::getWeight(nGP);
  if (!xg || !wg) return false;

  Matrix   Xnod, Jac;
  Vec4     X;
  Vec3     normal;

  std::vector<LR::Element*>::iterator el1 = m_basis[0]->elementBegin();
  for (int iel = 1; el1 != m_basis[0]->elementEnd(); ++el1, ++iel) {
    short int status = iChk.hasContribution(iel);
    if (!status) continue; // no interface contributions for this element
    // first push the hosting elements
    std::vector<size_t> els(1,iel);
    std::vector<size_t> elem_sizes(1,(*el1)->nBasisFunctions());

    double uh = ((*el1)->umin()+(*el1)->umax())/2.0;
    double vh = ((*el1)->vmin()+(*el1)->vmax())/2.0;

    for (size_t i=1; i < m_basis.size(); ++i) {
      els.push_back(m_basis[i]->getElementContaining(uh, vh)+1);
      elem_sizes.push_back(m_basis[i]->getElement(els.back()-1)->nBasisFunctions());
    }

    // Set up control point coordinates for current element
    if (!this->getElementCoordinates(Xnod,els[geoBasis-1]))
      return false;

    LocalIntegral* A = integrand.getLocalIntegral(elem_sizes, iel);
    integrand.initElement(MNPC[els[geoBasis-1]-1],elem_sizes,nb,*A);
    size_t origSize = A->vec.size();

    int bit = 8;
    for (int iedge = 4; iedge > 0 && status > 0; iedge--, bit /= 2) {
      if (status & bit) {
        const int edgeDir = (iedge+1)/(iedge%2 ? -2 : 2);
        const int t1 = abs(edgeDir);   // Tangent direction normal to the edge
        const int t2 = 3-abs(edgeDir); // Tangent direction along the edge

        // Set up parameters
        double u1 = iedge != 2 ? (*el1)->umin() : (*el1)->umax();
        double v1 = iedge < 4 ? (*el1)->vmin() : (*el1)->vmax();
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
        for (size_t i = 0; i < intersections.size(); ++i) {
          std::vector<double> parval(2);
          parval[0] = u1-epsu;
          parval[1] = v1-epsv;

          if (iedge == 1 || iedge == 2)
            v2 = intersections[i];
          else
            u2 = intersections[i];

          int el_neigh = this->getBasis(1)->getElementContaining(parval)+1;
          const LR::Element* el2 = m_basis[0]->getElement(el_neigh-1);
          uh = (el2->umin()+el2->umax())/2.0;
          vh = (el2->vmin()+el2->vmax())/2.0;

          std::vector<size_t> els2(1,el_neigh);
          std::vector<size_t> elem_sizes2(1,el2->nBasisFunctions());
          for (size_t i=1; i < m_basis.size(); ++i) {
            els2.push_back(m_basis[i]->getElementContaining(uh, vh)+1);
            elem_sizes2.push_back(m_basis[i]->getElement(els2.back()-1)->nBasisFunctions());
          }

          LocalIntegral* A_neigh = integrand.getLocalIntegral(elem_sizes2, el_neigh);
          integrand.initElement(MNPC[els2[geoBasis-1]-1],elem_sizes2,nb,*A_neigh);

          // Element sizes for both elements
          std::vector<size_t> elem_sizes3(elem_sizes);
          std::copy(elem_sizes2.begin(), elem_sizes2.end(),
                    std::back_inserter(elem_sizes3));

          MxFiniteElement fe(elem_sizes3);
          fe.h = this->getElementCorners(els2[geoBasis-1], fe.XC);

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
          if (!this->getElementCoordinates(Xnod2,els2[geoBasis-1]))
            return false;

          for (int g = 0; g < nGP; g++, ++fe.iGP)
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
              m_basis[b]->computeBasis(fe.u+epsu, fe.v+epsv, spline, els[b]-1);
              SplineUtils::extractBasis(spline,fe.basis(b+1),dNxdu[b]);
              m_basis[b]->computeBasis(fe.u-epsu, fe.v-epsv, spline, els2[b]-1);
              SplineUtils::extractBasis(spline,fe.basis(b+1+m_basis.size()),
                                        dNxdu[b+m_basis.size()]);
            }

            // Compute Jacobian inverse of the coordinate mapping and
            // basis function derivatives w.r.t. Cartesian coordinates
            fe.detJxW = utl::Jacobian(Jac2, normal,
                                      fe.grad(geoBasis+m_basis.size()),
                                      Xnod2,dNxdu[geoBasis-1+m_basis.size()],t1,t2);
            fe.detJxW = utl::Jacobian(Jac, normal,
                                      fe.grad(geoBasis),Xnod,dNxdu[geoBasis-1],t1,t2);
            if (fe.detJxW == 0.0) continue; // skip singular points
            for (size_t b = 0; b < m_basis.size(); ++b)
              if (b != (size_t)geoBasis-1) {
                fe.grad(b+1).multiply(dNxdu[b],Jac);
                fe.grad(b+1+m_basis.size()).multiply(dNxdu[b+m_basis.size()],Jac);
              }

            if (edgeDir < 0)
              normal *= -1.0;

            // Cartesian coordinates of current integration point
            X = Xnod * fe.basis(geoBasis);
            X.t = time.t;

            // Evaluate the integrand and accumulate element contributions
            fe.detJxW *= 0.5*dS*wg[g];
            if (!integrand.evalIntMx(*A,fe,time,X,normal))
              return false;
          }

          if (iedge == 1 || iedge == 2)
            v1 = v2;
          else
            u1 = u2;
        }
      }
    }
    // Finalize the element quantities
    if (!integrand.finalizeElement(*A,time,0))
      return false;

    // Assembly of global system integral
    if (!glInt.assemble(A,els[geoBasis-1]))
      return false;

    A->destruct();
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
    size_t ofs=0;
    Vector Ztmp;
    for (size_t j=0; j <  m_basis.size(); ++j) {
      if (nc[j] == 0)
        continue;
      // Fetch element containing evaluation point.
      // Sadly, points are not always ordered in the same way as the elements.
      int iel = m_basis[j]->getElementContaining(gpar[0][i],gpar[1][i]);

      // Evaluate basis function values/derivatives at current parametric point
      // and multiply with control point values to get the point-wise solution
      m_basis[j]->computeBasis(gpar[0][i],gpar[1][i],splinex[j],iel);

      std::vector<LR::Element*>::iterator el_it = m_basis[j]->elementBegin()+iel;
      Matrix val1(nc[j], splinex[j].basisValues.size());
      size_t col=1;
      for (auto* b : (*el_it)->support()) {
        for (size_t n = 1; n <= nc[j]; ++n)
          val1(n, col) = locSol(b->getId()*nc[j]+n+ofs);
        ++col;
      }
      Vector Ytmp;
      val1.multiply(splinex[j].basisValues,Ytmp);
      Ztmp.insert(Ztmp.end(),Ytmp.begin(),Ytmp.end());
      ofs += nb[j]*nc[j];
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
    std::vector<size_t> els;
    std::vector<size_t> elem_sizes;
    for (size_t b = 0; b < m_basis.size(); ++b) {
      els.push_back(m_basis[b]->getElementContaining(gpar[0][i],gpar[1][i])+1);
      elem_sizes.push_back((*(m_basis[b]->elementBegin()+els.back()-1))->nBasisFunctions());
    }

    // Evaluate the basis functions at current parametric point
    MxFiniteElement fe(elem_sizes);
    std::vector<Matrix> dNxdu(m_basis.size());
    Matrix Jac, Xnod;
    std::vector<Matrix3D> d2Nxdu2(m_basis.size());
    Matrix3D Hess;
    if (use2ndDer)
      for (size_t b = 0; b < m_basis.size(); ++b) {
        Go::BasisDerivsSf2 spline;
        m_basis[b]->computeBasis(gpar[0][i],gpar[1][i],spline,els[b]-1);
        SplineUtils::extractBasis(spline,fe.basis(b+1),dNxdu[b],d2Nxdu2[b]);
      }
    else
      for (size_t b = 0; b < m_basis.size(); ++b) {
        Go::BasisDerivsSf spline;
        m_basis[b]->computeBasis(gpar[0][i],gpar[1][i],spline,els[b]-1);
        SplineUtils::extractBasis(spline,fe.basis(b+1),dNxdu[b]);
      }

    // Set up control point (nodal) coordinates for current element
    if (!this->getElementCoordinates(Xnod,els[geoBasis-1])) return false;

    // Compute the Jacobian inverse
    fe.detJxW = utl::Jacobian(Jac,fe.grad(geoBasis),Xnod,dNxdu[geoBasis-1]);
    for (size_t b = 0; b < m_basis.size(); ++b)
      if (b != (size_t)geoBasis-1)
        fe.grad(b+1).multiply(dNxdu[b],Jac);

    // Compute Hessian of coordinate mapping and 2nd order derivatives
    if (use2ndDer) {
      if (!utl::Hessian(Hess,fe.hess(geoBasis),Jac,Xnod,
                        d2Nxdu2[geoBasis-1],fe.grad(geoBasis),true))
        return false;

      for (size_t b = 0; b < m_basis.size(); ++b)
        if (b != (size_t)geoBasis)
          utl::Hessian(Hess,fe.hess(b+1),Jac,Xnod,
                        d2Nxdu2[b],fe.grad(b+1),false);
    }

    // Now evaluate the solution field
    Vector solPt;
    if (!integrand.evalSol(solPt,fe,Xnod*fe.basis(geoBasis),
                           MNPC[els[geoBasis-1]-1],elem_sizes,nb))
      return false;
    else if (sField.empty())
      sField.resize(solPt.size(),nPoints,true);

    sField.fillColumn(1+i,solPt);
  }

  return true;
}


bool ASMu2Dmx::refine (const LR::RefineData& prm,
                       Vectors& sol, const char* fName)
{
  if (shareFE) return true;

  auto&& storeMesh = [this, fName]()
                     {
                       for (size_t b = 1; b <= m_basis.size(); ++b) {
                         std::stringstream str;
                         str << "_patch" << idx << "_basis" << b << "_" << fName;
                         std::ofstream paramMeshFile("param"+str.str());
                         this->getBasis(b)->writePostscriptMesh(paramMeshFile);
                         std::ofstream physMeshFile("physical"+str.str());
                         this->getBasis(b)->writePostscriptElements(physMeshFile);
                         std::ofstream pdotFile("param_dot"+str.str());
                         this->getBasis(b)->writePostscriptElements(pdotFile);
                         std::ofstream physdotFile("physical_dot"+str.str());
                         this->getBasis(b)->writePostscriptMeshWithControlPoints(physdotFile);
                       }
                     };

  if (!prm.errors.empty() || !prm.elements.empty()) {
    for (size_t j = 0; j < sol.size(); ++j) {
      size_t ofs = 0;
      for (size_t i = 0; i< m_basis.size(); ++i) {
        LR::extendControlPoints(m_basis[i].get(), sol[j], nfx[i], ofs);
        ofs += nfx[i]*nb[i];
      }
    }
  } else {
    if (fName)
      storeMesh();
    return true; // No refinement
  }

  if (doRefine(prm, refBasis.get())) {
    for (const LR::Meshline* line : refBasis->getAllMeshlines())
      for (size_t j = 0; j < m_basis.size(); ++j)
        if (refBasis == m_basis[j])
          continue;
        else {
          int p = m_basis[j]->order(line->span_u_line_ ? 1 : 0);
          int mult = 1;
          if (ASMmxBase::Type == ASMmxBase::REDUCED_CONT_RAISE_BASIS1 ||
              ASMmxBase::Type == ASMmxBase::REDUCED_CONT_RAISE_BASIS2) {
            if (line->multiplicity_ > 1)
              mult = p;
            else
              mult = (j == 0 && ASMmxBase::Type == REDUCED_CONT_RAISE_BASIS1) ||
                     (j == 1 && ASMmxBase::Type == REDUCED_CONT_RAISE_BASIS2) ? 2 : 1;
          }

          if (line->span_u_line_)
            m_basis[j]->insert_const_v_edge(line->const_par_,
                                            line->start_, line->stop_, mult);
          else
            m_basis[j]->insert_const_u_edge(line->const_par_,
                                            line->start_, line->stop_, mult);
        }

    // Uniformly refine to find basis 1
    if (ASMmxBase::Type == ASMmxBase::SUBGRID) {
      m_basis[0].reset(refBasis->copy());
      projBasis = m_basis[0];
      size_t nFunc = refBasis->nBasisFunctions();
      IntVec elems(nFunc);
      std::iota(elems.begin(),elems.end(),0);
      m_basis[0]->refineBasisFunction(elems);
    }

    size_t len = 0;
    for (size_t j = 0; j< m_basis.size(); ++j) {
      m_basis[j]->generateIDs();
      nb[j] = m_basis[j]->nBasisFunctions();
      len += nfx[j]*nb[j];
    }

    size_t ofs = 0;
    for (int i = sol.size()-1; i > 0; i--)
      for (size_t j = 0; j < m_basis.size(); ++j) {
        sol[i].resize(len);
        LR::contractControlPoints(m_basis[j].get(), sol[i], nfx[j], ofs);
        ofs += nfx[j]*nb[j];
      }

    if (fName)
      storeMesh();

  #ifdef SP_DEBUG
    std::cout <<"Refined mesh: ";
    for (const auto& it : m_basis)
      std::cout << it->nElements() <<" ";
    std::cout <<"elements ";
    for (const auto& it : m_basis)
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
  // TODO: Support for div-compatible
  int p1 = 0;
  for (size_t i = 1; i <= m_basis.size(); ++i)
    if (this->getBasis(i)->order(0) > p1) {
      threadBasis = this->getBasis(i);
      p1 = threadBasis->order(0);
    }

  LR::generateThreadGroups(threadGroups,threadBasis);
  if (silence || threadGroups[0].size() < 2) return;

  std::cout <<"\nMultiple threads are utilized during element assembly.";
  for (size_t i = 0; i < threadGroups[0].size(); i++)
    std::cout <<"\n Color "<< i+1 <<": "
              << threadGroups[0][i].size() <<" elements";
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


void ASMu2Dmx::remapErrors(RealArray& errors,
                           const RealArray& origErr, bool elemErrors) const
{
  const LR::LRSplineSurface* geo = this->getBasis(ASMmxBase::geoBasis);
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


size_t ASMu2Dmx::getNoProjectionNodes() const
{
  return projBasis->nBasisFunctions();
}


Fields* ASMu2Dmx::getProjectedFields(const Vector& coefs, size_t nf) const
{
  return new LRSplineFields2D(projBasis.get(), coefs, nf);
}


size_t ASMu2Dmx::getNoRefineNodes() const
{
  return refBasis->nBasisFunctions();
}


size_t ASMu2Dmx::getNoRefineElms() const
{
  return refBasis->nElements();
}


bool ASMu2Dmx::evalProjSolution (Matrix& sField, const Vector& locSol,
                                 const int* npe, int nf) const
{
#ifdef SP_DEBUG
  std::cout <<"ASMu2D::evalProjSolution(Matrix&,const Vector&,const int*,int)\n";
#endif
  if (projBasis == m_basis[0])
    return this->evalSolution(sField, locSol, npe, nf);

  // Compute parameter values of the result sampling points
  std::array<RealArray,2> gpar;
  for (int dir = 0; dir < 2; dir++)
    if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
      return false;

  size_t nComp = locSol.size() / this->getNoProjectionNodes();
  if (nComp*this->getNoProjectionNodes() != locSol.size())
    return false;

  size_t nPoints = gpar[0].size();
  if (nPoints != gpar[1].size())
    return false;

  Fields* f = this->getProjectedFields(locSol, nComp);

  // Evaluate the primary solution field at each point
  sField.resize(nComp,nPoints);
  for (size_t i = 0; i < nPoints; i++)
  {
    Vector vals;
    FiniteElement fe;
    fe.u = gpar[0][i];
    fe.v = gpar[1][i];
    f->valueFE(fe, vals);
    sField.fillColumn(1+i, vals);
  }

  delete f;

  return true;
}


const LR::LRSpline* ASMu2Dmx::getRefinementBasis() const
{
  return refBasis.get();
}
