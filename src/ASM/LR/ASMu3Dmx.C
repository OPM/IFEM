// $Id$
//==============================================================================
//!
//! \file ASMu3Dmx.C
//!
//! \date Mar 8 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for assembly of unstructured 3D spline mixed FE models.
//!
//==============================================================================

#include "ASMu3Dmx.h"

#include "GoTools/trivariate/SplineVolume.h"

#include "LRSpline/LRSplineVolume.h"
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
#include <numeric>
#include <utility>


ASMu3Dmx::ASMu3Dmx (const CharVec& n_f)
  : ASMu3D(std::accumulate(n_f.begin(), n_f.end(), 0)), ASMmxBase(n_f),
    bezierExtractmx(myBezierExtractmx)
{
  threadBasis = nullptr;
  myGeoBasis = ASMmxBase::itgBasis;
}


ASMu3Dmx::ASMu3Dmx (const ASMu3Dmx& patch, const CharVec& n_f)
  : ASMu3D(patch), ASMmxBase(n_f[0]==0?patch.nfx:n_f),
    m_basis(patch.m_basis),
    bezierExtractmx(patch.myBezierExtractmx)
{
  threadBasis = patch.threadBasis;
  nfx = patch.nfx;
  nb =  patch.nb;
  myGeoBasis = ASMmxBase::itgBasis;
}


const LR::LRSplineVolume* ASMu3Dmx::getBasis (int basis) const
{
  if (basis < 1)
    return this->ASMu3D::getBasis(basis);

  if (basis > static_cast<int>(m_basis.size()))
    return nullptr;

  return m_basis[basis-1].get();
}


LR::LRSplineVolume* ASMu3Dmx::getBasis (int basis)
{
  return const_cast<LR::LRSplineVolume*>(std::as_const(*this).getBasis(basis));
}


bool ASMu3Dmx::readBasis (std::istream& is, size_t basis)
{
  if (basis < 1 || basis > nfx.size())
    return false;

  if (m_basis.empty())
  {
    m_basis.resize(nfx.size());
    nb.resize(nfx.size(), 0);
  }

  m_basis[--basis] = std::make_shared<LR::LRSplineVolume>();
  is >> *m_basis[basis];
  nb[basis] = m_basis[basis]->nBasisFunctions();
  m_basis[basis]->generateIDs();

  return true;
}


void ASMu3Dmx::clear (bool retainGeometry)
{
  if (!retainGeometry) {
    // Erase the spline data
    for (SplinePtr& patch : m_basis)
      patch.reset();

    m_basis.clear();
  }

  // Erase the FE data
  this->ASMu3D::clear(retainGeometry);
}


size_t ASMu3Dmx::getNoNodes (int basis) const
{
  if (basis > (int)nb.size() || basis < 1)
    return this->ASMbase::getNoNodes(basis);

  return nb[basis-1];
}


unsigned char ASMu3Dmx::getNoFields (int basis) const
{
  if (basis > (int)m_basis.size() || basis < 0)
    basis = 0;

  if (basis == 0)
    return std::accumulate(nfx.begin(), nfx.end(), 0);

  return nfx[basis-1];
}


unsigned char ASMu3Dmx::getNodalDOFs (size_t inod) const
{
  if (this->isLMn(inod)) return nLag;
  size_t nbc=0;
  for (size_t i=0;i<nb.size();++i)
    if (inod <= (nbc+=nb[i]))
      return nfx[i];

  return nfx[0];
}


char ASMu3Dmx::getNodeType (size_t inod) const
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


void ASMu3Dmx::initMADOF (const int* sysMadof)
{
  this->initMx(MLGN,sysMadof);
}


void ASMu3Dmx::extractNodeVec (const RealArray& globRes, RealArray& nodeVec,
                               unsigned char, int basis) const
{
  this->extractNodeVecMx(globRes,nodeVec,basis);
}


bool ASMu3Dmx::injectNodeVec (const RealArray& nodeRes, RealArray& globRes,
                              unsigned char, int basis) const
{
  this->injectNodeVecMx(globRes,nodeRes,basis);
  return true;
}


bool ASMu3Dmx::getSolution (Matrix& sField, const Vector& locSol,
                            const IntVec& nodes) const
{
  return this->getSolutionMx(sField,locSol,nodes);
}


bool ASMu3Dmx::generateFEMTopology ()
{
  if (!myMLGN.empty())
    return true;

  if (m_basis.empty()) {
    VolumeVec vvec = ASMmxBase::establishBases(tensorspline, ASMmxBase::Type);
    for (size_t b = 0; b < vvec.size(); b++)
      m_basis.push_back(std::make_shared<LR::LRSplineVolume>(vvec[b].get()));

    // make a backup as establishBases resets it
    int geoB = ASMmxBase::itgBasis;
    std::shared_ptr<Go::SplineVolume> otherBasis =
        ASMmxBase::establishBases(tensorspline,
                                  ASMmxBase::FULL_CONT_RAISE_BASIS1).front();
    itgBasis = geoB;

    // we need to project on something that is not one of our bases
    if (!projB) {
      if (ASMmxBase::Type == ASMmxBase::REDUCED_CONT_RAISE_BASIS1 ||
          ASMmxBase::Type == ASMmxBase::REDUCED_CONT_RAISE_BASIS2 ||
          ASMmxBase::Type == ASMmxBase::DIV_COMPATIBLE)
        projB = std::make_shared<LR::LRSplineVolume>(otherBasis.get());
      else if (ASMmxBase::Type == ASMmxBase::SUBGRID)
        projB = m_basis.front();
      else // FULL_CONT_RAISE_BASISx
        projB = m_basis[2-ASMmxBase::itgBasis];
    }

    if (ASMmxBase::Type == ASMmxBase::SUBGRID)
      projB2 = refB = std::make_shared<LR::LRSplineVolume>(otherBasis.get());
    else
      refB = projB;

    delete tensorspline;
    tensorspline = nullptr;
  }
  lrspline = m_basis[itgBasis-1];
  projB->generateIDs();
  projB->getElementContaining(projB->getElement(0)->midpoint()); // to force cache generation
  if (projB2) {
    projB2->generateIDs();
    projB2->getElementContaining(projB2->getElement(0)->midpoint()); // to force cache generation
  }
  myGeoBasis = ASMmxBase::itgBasis;

  nb.clear();
  nb.reserve(m_basis.size());
  for (const SplinePtr& it : m_basis) {
    nb.push_back(it->nBasisFunctions());
#ifdef SP_DEBUG
    std::cout <<"Basis "<< nb.size()
              <<":\nnumCoefs: "<< nb.back()
              <<"\norder: "<< it->order(0) <<" "<<
                              it->order(1) <<" "<< it->order(2) << std::endl;
#endif
  }

  if (shareFE == 'F') return true;

  nel = m_basis[itgBasis-1]->nElements();
  nnod = std::accumulate(nb.begin(), nb.end(), 0);

  myMLGE.resize(nel,0);
  myMLGN.resize(nnod);
  myMNPC.resize(nel);
  for (SplinePtr& it : m_basis)
    it->generateIDs();

  size_t iel = 0;
  for (const LR::Element* el1 : m_basis[itgBasis-1]->getAllElements())
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

  size_t b = 0;
  myBezierExtractmx.resize(m_basis.size());
  for (const SplinePtr& it : m_basis) {
    myBezierExtractmx[b].resize(it->nElements());
    for (int iel = 0; iel < it->nElements(); iel++)
    {
      PROFILE("Bezier extraction");
      // Get bezier extraction matrix
      RealArray extrMat;
      it->getBezierExtraction(iel,extrMat);
      myBezierExtractmx[b][iel].resize(it->getElement(iel)->nBasisFunctions(),
                                       it->order(0)*it->order(1)*it->order(2));
      myBezierExtractmx[b][iel].fill(extrMat.data(),extrMat.size());
    }
    ++b;
  }

  for (size_t inod = 0; inod < nnod; ++inod)
    myMLGN[inod] = ++gNod;

#ifdef SP_DEBUG
  std::cout <<"NEL = "<< nel <<" NNOD = "<< nnod << std::endl;
#endif

  geomB = m_basis[itgBasis-1];

  return true;
}


bool ASMu3Dmx::integrate (Integrand& integrand,
                          GlobalIntegral& glInt,
                          const TimeDomain& time)
{
  if (m_basis.empty())
    return true; // silently ignore empty patches

  PROFILE2("ASMu3Dmx::integrate(I)");

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

  for (std::unique_ptr<ASMu3D::BasisFunctionCache>& cache : myCache) {
    cache->setIntegrand(&integrand);
    cache->init(use2ndDer ? 2 : 1);
  }

  ASMu3D::BasisFunctionCache& cache = *myCache.front();

  const std::array<int,3>& ng = cache.nGauss();
  const std::array<const double*,3>& xg = cache.coord();
  const std::array<const double*,3>& wg = cache.weight();

  ThreadGroups oneGroup;
  if (glInt.threadSafe()) oneGroup.oneGroup(nel);
  const IntMat& groups = glInt.threadSafe() ? oneGroup[0] : threadGroups[0];


  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t t = 0; t < groups.size() && ok; t++)
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < groups[t].size(); e++)
    {
      if (!ok)
        continue;

      const LR::Element* el = threadBasis->getElement(groups[t][e]);

      std::vector<int>    els;
      std::vector<size_t> elem_sizes;
      this->getElementsAt(el->midpoint(),els,elem_sizes);

      MxFiniteElement fe(elem_sizes);
      Matrix   Xnod, Jac;
      Matrix3D Hess;
      double   dXidu[3];
      double   param[3] = { 0.0, 0.0, 0.0 };
      Vec4     X(param,time.t);

      int iEl = el->getId();
      fe.iel = MLGE[iEl];

      // Get element volume in the parameter space
      double dV = 0.125*el->volume();

      // Set up control point (nodal) coordinates for current element
      if (!this->getElementCoordinates(Xnod,iEl+1))
      {
        ok = false;
        continue;
      }

      if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
        fe.h = this->getElementCorners(iEl+1, fe.XC);

      if (integrand.getIntegrandType() & Integrand::G_MATRIX)
        // Element size in parametric space
        for (int i = 0; i < 3; i++)
          dXidu[i] = el->getParmax(i) - el->getParmin(i);

      else if (integrand.getIntegrandType() & Integrand::AVERAGE)
      {
        // --- Compute average value of basis functions over the element -----

        fe.Navg.resize(elem_sizes[0],true);
        double vol = 0.0;
        size_t jp = 0;
        for (int k = 0; k < ng[2]; k++)
          for (int j = 0; j < ng[1]; j++)
            for (int i = 0; i < ng[0]; i++, jp++)
            {
              // Fetch basis function derivatives at current integration point
              const BasisFunctionVals& bfs = myCache[itgBasis-1]->getVals(iEl,jp);

              // Compute Jacobian determinant of coordinate mapping
              // and multiply by weight of current integration point
              double detJac = utl::Jacobian(Jac,fe.grad(itgBasis),
                                            Xnod,bfs.dNdu,false);
              double weight = dV*wg[0][i]*wg[1][j]*wg[2][k];

              // Numerical quadrature
              fe.Navg.add(bfs.N,detJac*weight);
              vol += detJac*weight;
        }

        // Divide by element volume
        fe.Navg /= vol;
      }

      else if (integrand.getIntegrandType() & Integrand::ELEMENT_CENTER)
      {
        // Compute the element center
        Go::Point X0;
        double u0 = 0.5*(el->getParmin(0) + el->getParmax(0));
        double v0 = 0.5*(el->getParmin(1) + el->getParmax(1));
        double w0 = 0.5*(el->getParmin(2) + el->getParmax(2));
        this->getBasis(itgBasis)->point(X0,u0,v0,w0);
        X.assign(SplineUtils::toVec3(X0));
      }

      // Initialize element quantities
      LocalIntegral* A = integrand.getLocalIntegral(elem_sizes,fe.iel);
      if (!integrand.initElement(MNPC[iEl],fe,elem_sizes,nb,*A))
      {
        A->destruct();
        ok = false;
        continue;
      }

      // --- Integration loop over all Gauss points in each direction ----------

      int jp = iEl*ng[0]*ng[1]*ng[2];
      fe.iGP = firstIp + jp; // Global integration point counter

      size_t ig = 0;
      for (int k = 0; k < ng[2]; k++)
        for (int j = 0; j < ng[1]; j++)
          for (int i = 0; i < ng[0]; i++, fe.iGP++, ig++)
          {
            // Local element coordinates of current integration point
            fe.xi   = xg[0][i];
            fe.eta  = xg[1][j];
            fe.zeta = xg[2][k];

            // Parameter values of current integration point
            fe.u = param[0] = cache.getParam(0,iEl,i);
            fe.v = param[1] = cache.getParam(1,iEl,j);
            fe.w = param[2] = cache.getParam(2,iEl,k);

            std::vector<const BasisFunctionVals*> bfs(this->getNoBasis());
            for (size_t b = 0; b < m_basis.size(); ++b) {
              bfs[b] = &myCache[b]->getVals(iEl,ig);
              fe.basis(b+1) = bfs[b]->N;
            }

            // Compute Jacobian inverse of the coordinate mapping and
            // basis function derivatives w.r.t. Cartesian coordinates
            if (!fe.Jacobian(Jac,Xnod,itgBasis,&bfs))
              continue; // skip singular points

            // Compute Hessian of coordinate mapping and 2nd order derivatives
            if (use2ndDer && !fe.Hessian(Hess,Jac,Xnod,itgBasis,&bfs))
              ok = false;

            // Compute G-matrix
            if (integrand.getIntegrandType() & Integrand::G_MATRIX)
              utl::getGmat(Jac,dXidu,fe.G);

            // Cartesian coordinates of current integration point
            X.assign(Xnod * fe.basis(itgBasis));

            // Evaluate the integrand and accumulate element contributions
            fe.detJxW *= dV*wg[0][i]*wg[1][j]*wg[2][k];
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

  for (std::unique_ptr<ASMu3D::BasisFunctionCache>& cache : myCache)
    cache->finalizeAssembly();
  return ok;
}


bool ASMu3Dmx::integrate (Integrand& integrand, int lIndex,
                          GlobalIntegral& glInt,
                          const TimeDomain& time)
{
  if (m_basis.empty())
    return true; // silently ignore empty patches

  PROFILE2("ASMu3Dmx::integrate(B)");

  // Get Gaussian quadrature points and weights
  int nGP = integrand.getBouIntegrationPoints(nGauss);
  const double* xg = GaussQuadrature::getCoord(nGP);
  const double* wg = GaussQuadrature::getWeight(nGP);
  if (!xg || !wg) return false;

  // Find the parametric direction of the face normal {-3,-2,-1, 1, 2, 3}
  const int faceDir = (lIndex+1)/((lIndex%2) ? -2 : 2);

  const int t1 = 1 + abs(faceDir)%3; // first tangent direction
  const int t2 = 1 + t1%3;           // second tangent direction

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;

  // Lambda function mapping face index to LR enum value
  auto&& getFaceEnum = [](int faceIndex) -> LR::parameterEdge
  {
    switch (faceIndex) {
    case 1: return LR::WEST;
    case 2: return LR::EAST;
    case 3: return LR::SOUTH;
    case 4: return LR::NORTH;
    case 5: return LR::BOTTOM;
    case 6: return LR::TOP;
    }
    return LR::NONE;
  };

  // Fetch all elements on the chosen face
  std::vector<LR::Element*> edgeElms;
  this->getBasis(itgBasis)->getEdgeElements(edgeElms,getFaceEnum(lIndex));


  // === Assembly loop over all elements on the patch face =====================

  for (LR::Element* el : edgeElms)
  {
    int iEl = el->getId();
    if (!myElms.empty() && !glInt.threadSafe() &&
        std::find(myElms.begin(), myElms.end(), iEl) == myElms.end())
      continue;

    std::vector<int>    els;
    std::vector<size_t> elem_sizes;
    this->getElementsAt(el->midpoint(),els,elem_sizes);

    MxFiniteElement fe(elem_sizes);
    fe.iel = MLGE[iEl];

    // Compute parameter values of the Gauss points over the whole element
    std::array<Vector,3> gpar;
    for (int d = 0; d < 3; d++)
      if (-1-d == faceDir)
      {
        gpar[d].resize(1);
        gpar[d].fill(this->getBasis(itgBasis)->startparam(d));
      }
      else if (1+d == faceDir)
      {
        gpar[d].resize(1);
        gpar[d].fill(this->getBasis(itgBasis)->endparam(d));
      }
      else
        this->getGaussPointParameters(gpar[d],d,nGP,iEl+1,xg);

    fe.xi = fe.eta = fe.zeta = faceDir < 0 ? -1.0 : 1.0;
    fe.u = gpar[0](1);
    fe.v = gpar[1](1);
    fe.w = gpar[2](1);

    std::vector<Matrix> dNxdu(m_basis.size());
    Matrix Xnod, Jac;

    double param[3] = { fe.u, fe.v, fe.w };
    Vec4   X(param,time.t);
    Vec3   normal;
    double dXidu[3];

    // Get element face area in the parameter space
    double dA = 0.25*this->getParametricArea(iEl+1,abs(faceDir));
    if (dA < 0.0) return false; // topology error (probably logic error)

    // Set up control point coordinates for current element
    if (!this->getElementCoordinates(Xnod,iEl+1))
      return false;

    if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
      fe.h = this->getElementCorners(iEl+1,fe.XC);

    if (integrand.getIntegrandType() & Integrand::G_MATRIX)
      // Element size in parametric space
      for (int i = 0; i < 3; i++)
        dXidu[i] = el->getParmax(i) - el->getParmin(i);

    // Initialize element quantities
    LocalIntegral* A = integrand.getLocalIntegral(elem_sizes,fe.iel,true);
    bool ok = integrand.initElementBou(MNPC[iEl],elem_sizes,nb,*A);


    // --- Integration loop over all Gauss points in each direction ------------

    fe.iGP = firstp; // Global integration point counter
    firstp += nGP*nGP;

    for (int j = 0; j < nGP; j++)
      for (int i = 0; i < nGP && ok; i++, fe.iGP++)
      {
        // Local element coordinates and parameter values
        // of current integration point
        int k1, k2, k3;
        switch (abs(faceDir))
        {
          case 1: k2 = i; k3 = j; k1 = 0; break;
          case 2: k1 = i; k3 = j; k2 = 0; break;
          case 3: k1 = i; k2 = j; k3 = 0; break;
          default: k1 = k2 = k3 = 0;
        }
        if (gpar[0].size() > 1)
        {
          fe.xi = xg[k1];
          fe.u = param[0] = gpar[0](k1+1);
        }
        if (gpar[1].size() > 1)
        {
          fe.eta = xg[k2];
          fe.v = param[1] = gpar[1](k2+1);
        }
        if (gpar[2].size() > 1)
        {
          fe.zeta = xg[k3];
          fe.w = param[2] = gpar[2](k3+1);
        }

        // Fetch basis function derivatives at current integration point
        for (size_t b = 1; b <= m_basis.size(); ++b)
          this->evaluateBasis(iEl, fe, dNxdu[b-1], b);

        // Compute basis function derivatives and the face normal
        fe.detJxW = utl::Jacobian(Jac, normal, fe.grad(itgBasis), Xnod,
                                  dNxdu[itgBasis-1], t1, t2);
        if (fe.detJxW == 0.0) continue; // skip singular points

        for (size_t b = 1; b <= m_basis.size(); ++b)
          if ((int)b != itgBasis)
            fe.grad(b).multiply(dNxdu[b-1],Jac);

        if (faceDir < 0) normal *= -1.0;

        // Compute G-matrix
        if (integrand.getIntegrandType() & Integrand::G_MATRIX)
          utl::getGmat(Jac,dXidu,fe.G);

        // Cartesian coordinates of current integration point
        X.assign(Xnod * fe.basis(itgBasis));

        // Evaluate the integrand and accumulate element contributions
        fe.detJxW *= dA*wg[i]*wg[j];
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


bool ASMu3Dmx::evalSolution (Matrix& sField, const Vector& locSol,
                             const RealArray* gpar, bool,
                             int deriv, int nf) const
{
  size_t nPoints = gpar[0].size();
  if (nPoints != gpar[1].size() || nPoints != gpar[2].size())
    return false;

  Vector   ptSol;
  std::vector<Matrix> dNxdu(m_basis.size()), dNxdX(m_basis.size());
  Matrix   Jac, Xnod, eSol, ptDer;

  std::vector<Go::BasisPts> splinex(m_basis.size());

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
      int iel = m_basis[j]->getElementContaining(gpar[0][i],gpar[1][i],gpar[2][i]);
      const LR::Element* el = m_basis[j]->getElement(iel);

      // Evaluate basis function values/derivatives at current parametric point
      // and multiply with control point values to get the point-wise solution
      m_basis[j]->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],splinex[j],iel);

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


bool ASMu3Dmx::evalSolution (Matrix& sField, const IntegrandBase& integrand,
                             const RealArray* gpar, bool) const
{
#ifdef SP_DEBUG
  std::cout <<"ASMu3Dmx::evalSolution(Matrix&,const IntegrandBase&,const RealArray*,bool)\n";
#endif

  sField.resize(0,0);

  // TODO: investigate the possibility of doing "regular" refinement by
  //       uniform tesselation grid and ignoring LR mesh lines

  size_t nPoints = gpar[0].size();
  bool use2ndDer = integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES;
  if (nPoints != gpar[1].size() || nPoints != gpar[2].size())
    return false;

  // Evaluate the secondary solution field at each point
  for (size_t i = 0; i < nPoints; i++)
  {
    // Fetch element containing evaluation point
    // sadly, points are not always ordered in the same way as the elements
    std::vector<int>    els;
    std::vector<size_t> elem_sizes;
    this->getElementsAt({gpar[0][i],gpar[1][i],gpar[2][i]},els,elem_sizes);

    // Evaluate the basis functions at current parametric point
    MxFiniteElement       fe(elem_sizes,firstIp+i);
    std::vector<Matrix>   dNxdu(m_basis.size());
    std::vector<Matrix3D> d2Nxdu2(use2ndDer ? m_basis.size() : 0);
    Matrix Jac, Xnod;
    Matrix3D Hess;
    if (use2ndDer)
      for (size_t b = 0; b < m_basis.size(); ++b) {
        Go::BasisDerivs2 spline;
        m_basis[b]->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline,els[b]-1);
        SplineUtils::extractBasis(spline,fe.basis(b+1),dNxdu[b],d2Nxdu2[b]);
      }
    else
      for (size_t b = 0; b < m_basis.size(); ++b) {
        Go::BasisDerivs spline;
        m_basis[b]->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline,els[b]-1);
        SplineUtils::extractBasis(spline,fe.basis(b+1),dNxdu[b]);
      }

    // Set up control point (nodal) coordinates for current element
    if (!this->getElementCoordinates(Xnod,els[itgBasis-1])) return false;

    // Compute Jacobian inverse of the coordinate mapping and
    // basis function derivatives w.r.t. Cartesian coordinates
    if (!fe.Jacobian(Jac,Xnod,itgBasis,nullptr,&dNxdu))
      continue; // skip singular points

    // Compute Hessian of coordinate mapping and 2nd order derivatives
    if (use2ndDer && !fe.Hessian(Hess,Jac,Xnod,itgBasis,nullptr,&d2Nxdu2))
      return false;

    // Cartesian coordinates of current integration point
    fe.u = gpar[0][i];
    fe.v = gpar[1][i];
    fe.w = gpar[2][i];
    utl::Point X4(Xnod*fe.basis(itgBasis), {fe.u, fe.v, fe.w});

    // Now evaluate the solution field
    Vector solPt;
    if (!integrand.evalSol(solPt,fe,X4,
                           MNPC[els[itgBasis-1]-1],elem_sizes,nb))
      return false;
    else if (sField.empty())
      sField.resize(solPt.size(),nPoints,true);

    sField.fillColumn(1+i,solPt);
  }

  return true;
}


bool ASMu3Dmx::refine (const LR::RefineData& prm, Vectors& sol)
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

  if (doRefine(prm, this->getBasis(ASM::REFINEMENT_BASIS))) {
    for (size_t j = 0; j < m_basis.size(); ++j)
      if (refB != m_basis[j]) {
        if ((j == 0 && ASMmxBase::Type == REDUCED_CONT_RAISE_BASIS1) ||
            (j == 1 && ASMmxBase::Type == REDUCED_CONT_RAISE_BASIS2))
          this->copyRefinement(m_basis[j].get(), 2);
        else
          this->copyRefinement(m_basis[j].get(), 1);
      }

    // Uniformly refine to find basis 1
    if (ASMmxBase::Type == ASMmxBase::SUBGRID) {
      m_basis[0].reset(this->getBasis(ASM::REFINEMENT_BASIS)->copy());
      projB = m_basis.front();
      size_t nFunc = refB->nBasisFunctions();
      IntVec elems(nFunc);
      std::iota(elems.begin(),elems.end(),0);
      m_basis[0]->refineBasisFunction(elems);
    }

    if (projB2)
      projB2->generateIDs();

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
              << projB->nElements() << " elements "
              << projB->nBasisFunctions() << " nodes" << std::endl;
    std::cout << "Refinement basis: "
              << refB->nElements() << " elements "
              << refB->nBasisFunctions() << " nodes" << std::endl;
  #endif

    return true;
  }

  return false;
}


Vec3 ASMu3Dmx::getCoord (size_t inod) const
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


void ASMu3Dmx::generateThreadGroups (const Integrand& integrand, bool silence,
                                     bool ignoreGlobalLM)
{
  int p1 = 0;
  if (ASMmxBase::Type == ASMmxBase::DIV_COMPATIBLE)
    threadBasis = this->getBasis(4);
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
    secConstraint = {this->getBasis(1),this->getBasis(2),this->getBasis(3)};

  LR::generateThreadGroups(threadGroups,threadBasis,secConstraint);
  LR::generateThreadGroups(projThreadGroups,projB.get());
  if (projB2)
    LR::generateThreadGroups(proj2ThreadGroups,projB2.get());

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


void ASMu3Dmx::remapErrors (RealArray& errors,
                            const RealArray& origErr, bool elemErrors) const
{
  const LR::LRSplineVolume* geo = this->getBasis(ASMmxBase::itgBasis);
  for (const LR::Element* elm : geo->getAllElements()) {
    int rEl = refB->getElementContaining(elm->midpoint());
    if (elemErrors)
      errors[rEl] += origErr[elm->getId()];
    else
      for (LR::Basisfunction* b : refB->getElement(rEl)->support())
        errors[b->getId()] += origErr[elm->getId()];
  }
}


size_t ASMu3Dmx::getNoRefineNodes() const
{
  return refB->nBasisFunctions();
}


size_t ASMu3Dmx::getNoRefineElms() const
{
  return refB->nElements();
}


void ASMu3Dmx::copyRefinement (LR::LRSplineVolume* basis,
                               int multiplicity) const
{
  const LR::LRSplineVolume* ref = this->getBasis(ASM::REFINEMENT_BASIS);
  for (const LR::MeshRectangle* rect : ref->getAllMeshRectangles()) {
    int mult = rect->multiplicity_ > 1 ? basis->order(rect->constDirection())
                                       : multiplicity;
    LR::MeshRectangle* newRect = rect->copy();
    newRect->multiplicity_ = mult;

    basis->insert_line(newRect);
  }
}


void ASMu3Dmx::swapProjectionBasis ()
{
  if (projB2) {
    ASMmxBase::itgBasis = ASMmxBase::itgBasis == 1 ? 2 : 1;
    std::swap(projB, projB2);
    std::swap(projThreadGroups, proj2ThreadGroups);
  }
}


void ASMu3Dmx::getElementsAt (const RealArray& param,
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


BasisFunctionVals ASMu3Dmx::BasisFunctionCache::calculatePt (size_t el,
                                                             size_t gp,
                                                             bool reduced) const
{
  PROFILE2("Spline evaluation");
  const std::array<size_t,3> gpIdx = this->gpIndex(gp,reduced);
  FiniteElement fe;
  fe.u = this->getParam(0,el,gpIdx[0],reduced);
  fe.v = this->getParam(1,el,gpIdx[1],reduced);
  fe.w = this->getParam(2,el,gpIdx[2],reduced);

  const ASMu3Dmx& pch = static_cast<const ASMu3Dmx&>(patch);

  const LR::Element* elm = pch.lrspline->getElement(el);
  std::array<double,3> du;
  du[0] = elm->umax() - elm->umin();
  du[1] = elm->vmax() - elm->vmin();
  du[2] = elm->wmax() - elm->wmin();

  el = pch.getBasis(basis)->getElementContaining(elm->midpoint());

  return this->calculatePrm(fe,du,el,gp,reduced);
}


bool ASMu3Dmx::separateProjectionBasis () const
{
  return std::none_of(m_basis.begin(), m_basis.end(),
                      [this](const SplinePtr& entry)
                      { return entry == projB; });
}
