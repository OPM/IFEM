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

#include <array>
#include <numeric>


ASMu2Dmx::ASMu2Dmx (unsigned char n_s,
		    const std::vector<unsigned char>& n_f)
  : ASMu2D(n_s), ASMmxBase(n_f)
{
}


ASMu2Dmx::ASMu2Dmx (const ASMu2Dmx& patch,
                    const std::vector<unsigned char>& n_f)
  : ASMu2D(patch), ASMmxBase(n_f[0]==0?patch.nfx:n_f)
{
  m_basis = patch.m_basis;
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
  return false;
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
  if (m_basis.empty()) {
    auto vec = ASMmxBase::establishBases(tensorspline, ASMmxBase::Type);
    m_basis.resize(vec.size());
    for (size_t i=0;i<vec.size();++i)
      m_basis[i].reset(new LR::LRSplineSurface(vec[i].get()));
  }
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

  std::vector<Matrix> dNxdu(m_basis.size());
  Matrix   Xnod, Jac;
  Vec4     X;

  // === Assembly loop over all elements in the patch ==========================

  std::vector<LR::Element*>::iterator el1 = m_basis[geoBasis-1]->elementBegin();
  for (int iel = 1; el1 != m_basis[geoBasis-1]->elementEnd(); ++el1, ++iel)
  {
    double uh = ((*el1)->umin()+(*el1)->umax())/2.0;
    double vh = ((*el1)->vmin()+(*el1)->vmax())/2.0;
    std::vector<size_t> els;
    std::vector<size_t> elem_sizes;
    for (size_t i=0; i < m_basis.size(); ++i) {
      els.push_back(m_basis[i]->getElementContaining(uh, vh)+1);
      elem_sizes.push_back((*(m_basis[i]->elementBegin()+els.back()-1))->nBasisFunctions());
    }

    int geoEl = els[geoBasis-1];

    MxFiniteElement fe(elem_sizes);
    fe.iel = MLGE[iel-1];

    // Get element area in the parameter space
    double dA = this->getParametricArea(geoEl);
    if (dA < 0.0) return false; // topology error (probably logic error)

    // Set up control point (nodal) coordinates for current element
    if (!this->getElementCoordinates(Xnod,geoEl))
      return false;

    // Compute parameter values of the Gauss points over this element
    std::array<RealArray,2> gpar;
    for (int d = 0; d < 2; d++)
      this->getGaussPointParameters(gpar[d],d,nGauss,geoEl,xg);

    // Initialize element quantities
    LocalIntegral* A = integrand.getLocalIntegral(elem_sizes,fe.iel,false);
    if (!integrand.initElement(MNPC[iel-1], elem_sizes, nb, *A))
    {
      A->destruct();
      return false;
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
	fe.u = gpar[0][i];
	fe.v = gpar[1][j];

	// Compute basis function derivatives at current integration point
        std::vector<Go::BasisDerivsSf> splinex(m_basis.size());
        for (size_t i=0; i < m_basis.size(); ++i) {
          m_basis[i]->computeBasis(fe.u, fe.v, splinex[i], els[i]-1);
          SplineUtils::extractBasis(splinex[i],fe.basis(i+1),dNxdu[i]);
        }

	// Compute Jacobian inverse of coordinate mapping and derivatives
        // basis function derivatives w.r.t. Cartesian coordinates
        fe.detJxW = utl::Jacobian(Jac,fe.grad(geoBasis),Xnod,dNxdu[geoBasis-1]);
	if (fe.detJxW == 0.0) continue; // skip singular points
        for (size_t b = 0; b < m_basis.size(); ++b)
          if (b != (size_t)geoBasis-1)
            fe.grad(b+1).multiply(dNxdu[b],Jac);

	// Cartesian coordinates of current integration point
        X = Xnod * fe.basis(geoBasis);
	X.t = time.t;

	// Evaluate the integrand and accumulate element contributions
	fe.detJxW *= 0.25*dA*wg[i]*wg[j];
	PROFILE3("Integrand::evalInt");
	if (!integrand.evalIntMx(*A,fe,time,X))
	  return false;
      }

    // Finalize the element quantities
    if (!integrand.finalizeElement(*A,time,firstIp+jp))
      return false;

    // Assembly of global system integral
    if (!glInt.assemble(A->ref(),fe.iel))
      return false;

    A->destruct();
  }

  return true;
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
  const int edgeDir = (lIndex+1)/(lIndex%2 ? -2 : 2);

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

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;

  std::vector<Matrix> dNxdu(m_basis.size());
  Matrix Xnod, Jac;
  Vec4   X;
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
      elem_sizes.push_back((*(m_basis[i]->elementBegin()+els.back()-1))->nBasisFunctions());
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
    fe.iel = MLGE[iel-1];
    fe.xi = fe.eta = edgeDir < 0 ? -1.0 : 1.0;
    LocalIntegral* A = integrand.getLocalIntegral(elem_sizes,fe.iel,true);
    if (!integrand.initElementBou(MNPC[iel-1], elem_sizes, nb, *A))
    {
      A->destruct();
      return false;
    }

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
      fe.u = gpar[0][i];
      fe.v = gpar[1][i];

      // Evaluate basis function derivatives at current integration points
      std::vector<Go::BasisDerivsSf> splinex(m_basis.size());
      for (size_t i=0; i < m_basis.size(); ++i) {
        m_basis[i]->computeBasis(fe.u, fe.v, splinex[i], els[i]-1);
        SplineUtils::extractBasis(splinex[i],fe.basis(i+1),dNxdu[i]);
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
      X = Xnod * fe.basis(geoBasis);
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


bool ASMu2Dmx::evalSolution (Matrix& sField, const Vector& locSol,
                             const RealArray* gpar, bool regular,
                             int deriv) const
{
  size_t nPoints = gpar[0].size();
  if (nPoints != gpar[1].size())
    return false;

  Vector   ptSol;
  std::vector<Matrix> dNxdu(m_basis.size()), dNxdX(m_basis.size());
  Matrix   Jac, Xnod, eSol, ptDer;

  std::vector<Go::BasisPtsSf> splinex(m_basis.size());

  // Evaluate the primary solution field at each point
  sField.resize(std::accumulate(nfx.begin(), nfx.end(), 0), nPoints);
  for (size_t i = 0; i < nPoints; i++)
  {
    size_t ofs=0;
    Vector Ztmp;
    for (size_t j=0; j < m_basis.size(); ++j) {
      // Fetch element containing evaluation point.
      // Sadly, points are not always ordered in the same way as the elements.
      int iel = m_basis[j]->getElementContaining(gpar[0][i],gpar[1][i]);

      // Evaluate basis function values/derivatives at current parametric point
      // and multiply with control point values to get the point-wise solution
      m_basis[j]->computeBasis(gpar[0][i],gpar[1][i],splinex[j],iel);

      std::vector<LR::Element*>::iterator el_it = m_basis[j]->elementBegin()+iel;
      Matrix val1(nfx[j], splinex[j].basisValues.size());
      size_t col=1;
      for (auto* b : (*el_it)->support()) {
        for (size_t n=1;n<=nfx[j];++n)
          val1(n, col) = locSol(b->getId()*nfx[j]+n+ofs);
        ++col;
      }
      Vector Ytmp;
      val1.multiply(splinex[j].basisValues,Ytmp);
      Ztmp.insert(Ztmp.end(),Ytmp.begin(),Ytmp.end());
      ofs += nb[j]*nfx[j];
    }

    sField.fillColumn(i+1, Ztmp);
  }

  return true;
}


bool ASMu2Dmx::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			     const RealArray* gpar, bool regular) const
{
  return evalSolution(sField,
                     const_cast<IntegrandBase&>(integrand).getSolution(0),
                     gpar, regular, 0);
}


bool ASMu2Dmx::refine (const LR::RefineData& prm,
                       Vectors& sol, const char* fName)
{
  if (shareFE) return true;

  char doRefine = 0;
  if (!prm.errors.empty())
    doRefine = 'E'; // Refine based on error indicators
  else if (!prm.elements.empty())
    doRefine = 'I'; // Refine the specified elements
  else
    return true; // No refinement

  // which basis to refine
  size_t bas = geoBasis-1;

  // to pick up if LR splines get stuck while doing refinement
  // print entry and exit point of this function

  double beta         = prm.options.size() > 0 ? prm.options[0]/100.0 : 0.10;
  int    multiplicity = prm.options.size() > 1 ? prm.options[1]       : 1;
  if (multiplicity > 1)
    for (int d = 0; d < geo->nVariate(); d++) {
      int p = geo->order(d) - 1;
      if (multiplicity > p) multiplicity = p;
    }

  enum refinementStrategy strat = LR_FULLSPAN;
  if (prm.options.size() > 2)
    switch (prm.options[2]) {
    case 1: strat = LR_MINSPAN; break;
    case 2: strat = LR_STRUCTURED_MESH; break;
    }

  int  maxTjoints     = prm.options.size() > 4 ? prm.options[4]      : -1;
  int  maxAspectRatio = prm.options.size() > 5 ? prm.options[5]      : -1;
  bool closeGaps      = prm.options.size() > 6 ? prm.options[6] != 0 : false;

  // set refinement parameters
  if (maxTjoints > 0)
    m_basis[bas]->setMaxTjoints(maxTjoints);
  if (maxAspectRatio > 0)
    m_basis[bas]->setMaxAspectRatio((double)maxAspectRatio);
  m_basis[bas]->setCloseGaps(closeGaps);
  m_basis[bas]->setRefMultiplicity(multiplicity);
  m_basis[bas]->setRefStrat(strat);

  for (size_t j = 0; j < sol.size(); ++j) {
    size_t ofs = 0;
    for (size_t i = 0; i< m_basis.size(); ++i) {
      LR::extendControlPoints(m_basis[i].get(), sol[j], nfx[i], ofs);
      ofs += nfx[i]*nb[i];
    }
  }

  // do actual refinement
  if (doRefine == 'E')
    m_basis[bas]->refineByDimensionIncrease(prm.errors,beta);
  else if(strat == LR_STRUCTURED_MESH)
    m_basis[bas]->refineBasisFunction(prm.elements);
  else
    m_basis[bas]->refineElement(prm.elements);

  const std::vector<LR::Meshline*> lines = m_basis[bas]->getAllMeshlines();
  for (auto line : lines)
    for (size_t j = 0; j < m_basis.size(); ++j)
      if (j == bas)
        continue;
      else if (line->span_u_line_)
        m_basis[j]->insert_const_u_edge(line->const_par_,
                                        line->start_, line->stop_,
                                        line->multiplicity_);
      else
        m_basis[j]->insert_const_v_edge(line->const_par_,
                                        line->start_, line->stop_,
                                        line->multiplicity_);

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

#ifdef SP_DEBUG
  std::cout <<"Refined mesh: ";
  for (const auto& it : m_basis)
    std::cout << it->nElements() <<" ";
  std::cout <<"elements";
  for (const auto& it : m_basis)
    std::cout << it->nBasisFunctions() <<" ";
  std::cout <<"nodes."<< std::endl;
#endif

  return true;
}
