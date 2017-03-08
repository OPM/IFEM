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

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/trivariate/SplineVolume.h"

#include "LRSpline/LRSplineVolume.h"
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


ASMu3Dmx::ASMu3Dmx (const CharVec& n_f)
  : ASMu3D(std::accumulate(n_f.begin(), n_f.end(), 0)), ASMmxBase(n_f),
    bezierExtract(myBezierExtract)
{
  myGeoBasis = ASMmxBase::geoBasis;
}


ASMu3Dmx::ASMu3Dmx (const ASMu3Dmx& patch, const CharVec& n_f)
  : ASMu3D(patch), ASMmxBase(n_f[0]==0?patch.nfx:n_f),
    bezierExtract(patch.myBezierExtract)
{
  m_basis = patch.m_basis;
  nfx = patch.nfx;
  nb =  patch.nb;
  myGeoBasis = ASMmxBase::geoBasis;
}


const LR::LRSplineVolume* ASMu3Dmx::getBasis (int basis) const
{
  if (basis < 1 || basis > (int)m_basis.size())
    return nullptr;

  return m_basis[basis-1].get();
}


LR::LRSplineVolume* ASMu3Dmx::getBasis (int basis)
{
  if (basis < 1 || basis > (int)m_basis.size())
    return nullptr;

  return m_basis[basis-1].get();
}


bool ASMu3Dmx::write (std::ostream& os, int basis) const
{
  return false;
}


void ASMu3Dmx::clear (bool retainGeometry)
{
  if (!retainGeometry) {
    // Erase the spline data
    for (auto& patch : m_basis)
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


void ASMu3Dmx::extractNodeVec (const Vector& globRes, Vector& nodeVec,
                               unsigned char, int basis) const
{
  this->extractNodeVecMx(globRes,nodeVec,basis);
}


bool ASMu3Dmx::injectNodeVec (const Vector& nodeRes, Vector& globRes,
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
  if (m_basis.empty()) {
    auto vec = ASMmxBase::establishBases(tensorspline, ASMmxBase::Type);
    m_basis.resize(vec.size());
    for (size_t i=0;i<vec.size();++i)
      m_basis[i].reset(new LR::LRSplineVolume(vec[i].get()));
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
    std::cout <<"\norder: "<< it->order(0) <<" "<<
                              it->order(1) <<" "<< it->order(2) << std::endl;
  }
#endif

  nel = m_basis[geoBasis-1]->nElements();

  nnod = std::accumulate(nb.begin(), nb.end(), 0);

  myMLGE.resize(nel,0);
  myMLGN.resize(nnod);
  myMNPC.resize(nel);
  for (auto& it : m_basis)
    it->generateIDs();

  std::vector<LR::Element*>::iterator el_it1 = m_basis[geoBasis-1]->elementBegin();
  for (size_t iel=0; iel<nel; iel++, ++el_it1)
  {
    double uh = ((*el_it1)->umin()+(*el_it1)->umax())/2.0;
    double vh = ((*el_it1)->vmin()+(*el_it1)->vmax())/2.0;
    double wh = ((*el_it1)->wmin()+(*el_it1)->wmax())/2.0;
    size_t nfunc = 0;
    for (size_t i=0; i<m_basis.size();++i) {
      auto el_it2 = m_basis[i]->elementBegin() +
                    m_basis[i]->getElementContaining(uh, vh, wh);
      nfunc += (*el_it2)->nBasisFunctions();
    }
    myMLGE[iel] = ++gEl; // global element number over all patches
    myMNPC[iel].resize(nfunc);

    int lnod = 0;
    size_t ofs=0;
    for (size_t i=0; i<m_basis.size();++i) {
      auto el_it2 = m_basis[i]->elementBegin() +
                    m_basis[i]->getElementContaining(uh, vh, wh);
      for (LR::Basisfunction *b : (*el_it2)->support())
        myMNPC[iel][lnod++] = b->getId()+ofs;
      ofs += nb[i];
    }
  }

  myBezierExtract.resize(m_basis.size());
  for (size_t b = 1; b <= m_basis.size(); ++b) {
    myBezierExtract[b-1].resize(this->getBasis(b)->nElements());
    std::vector<LR::Element*>::const_iterator eit = this->getBasis(b)->elementBegin();
    for (int iel = 0; iel < this->getBasis(b)->nElements(); iel++, ++eit)
    {
      PROFILE("Bezier extraction");
      int p1 = this->getBasis(b)->order(0);
      int p2 = this->getBasis(b)->order(1);
      int p3 = this->getBasis(b)->order(2);

      // Get bezier extraction matrix
      RealArray extrMat;
      this->getBasis(b)->getBezierExtraction(iel,extrMat);
      myBezierExtract[b-1][iel].resize((*eit)->nBasisFunctions(),p1*p2*p3);
      myBezierExtract[b-1][iel].fill(extrMat.data(),extrMat.size());
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


bool ASMu3Dmx::integrate (Integrand& integrand,
                          GlobalIntegral& glInt,
                          const TimeDomain& time)
{
  PROFILE2("ASMu3Dmx::integrate(I)");

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // evaluate all gauss points on the bezier patch (-1, 1)
  std::vector<Matrix> BN(m_basis.size());
  std::vector<Matrix> BdNdu(m_basis.size());
  std::vector<Matrix> BdNdv(m_basis.size());
  std::vector<Matrix> BdNdw(m_basis.size());
  for (size_t b = 1; b <= m_basis.size(); ++b) {
    const LR::LRSplineVolume* lrspline = this->getBasis(b);
    int p1 = lrspline->order(0);
    int p2 = lrspline->order(1);
    int p3 = lrspline->order(2);
    Go::BsplineBasis basis1 = getBezierBasis(p1);
    Go::BsplineBasis basis2 = getBezierBasis(p2);
    Go::BsplineBasis basis3 = getBezierBasis(p3);

    BN[b-1].resize(p1*p2*p3, nGauss*nGauss*nGauss);
    BdNdu[b-1].resize(p1*p2*p3, nGauss*nGauss*nGauss);
    BdNdv[b-1].resize(p1*p2*p3, nGauss*nGauss*nGauss);
    BdNdw[b-1].resize(p1*p2*p3, nGauss*nGauss*nGauss);
    int ig=1; // gauss point iterator
    for(int zeta=0; zeta<nGauss; zeta++) {
      for(int eta=0; eta<nGauss; eta++) {
        for(int xi=0; xi<nGauss; xi++, ig++) {
          double u[2*p1];
          double v[2*p2];
          double w[2*p3];
          basis1.computeBasisValues(xg[xi],   u, 1);
          basis2.computeBasisValues(xg[eta],  v, 1);
          basis3.computeBasisValues(xg[zeta], w, 1);
          int ib=1; // basis function iterator
          double sum = 0;
          for(int k=0; k<p3; k++) {
            for(int j=0; j<p2; j++) {
              for(int i=0; i<p1; i++, ib++) {
                BN[b-1](ib,ig)    = u[2*i  ]*v[2*j  ]*w[2*k  ];
                BdNdu[b-1](ib,ig) = u[2*i+1]*v[2*j  ]*w[2*k  ];
                BdNdv[b-1](ib,ig) = u[2*i  ]*v[2*j+1]*w[2*k  ];
                BdNdw[b-1](ib,ig) = u[2*i  ]*v[2*j  ]*w[2*k+1];
                sum += BN[b-1](ib,ig);
              }
            }
          }
        }
      }
    }
  }

  // Get the reduced integration quadrature points, if needed
  const double* xr = nullptr;
  const double* wr = nullptr;
  int nRed = integrand.getReducedIntegration(nGauss);
  if (nRed > 0)
  {
    xr = GaussQuadrature::getCoord(nRed);
    wr = GaussQuadrature::getWeight(nRed);
    if (!xr || !wr) return false;
  }
  else if (nRed < 0)
    nRed = nGauss; // The integrand needs to know nGauss


  // === Assembly loop over all elements in the patch ==========================

  bool ok=true;
  for(LR::Element *el : this->getBasis(geoBasis)->getAllElements()) {
    double uh = (el->umin()+el->umax())/2.0;
    double vh = (el->vmin()+el->vmax())/2.0;
    double wh = (el->wmin()+el->wmax())/2.0;
    std::vector<size_t> els;
    std::vector<size_t> elem_sizes;
    for (size_t i=0; i < m_basis.size(); ++i) {
      els.push_back(m_basis[i]->getElementContaining(uh, vh, wh)+1);
      elem_sizes.push_back((*(m_basis[i]->elementBegin()+els.back()-1))->nBasisFunctions());
    }
    int iEl = el->getId();
    MxFiniteElement fe(elem_sizes);
    fe.iel = iEl+1;

    std::vector<Matrix> dNxdu(m_basis.size());
    Matrix   Xnod, Jac;
    std::vector<Matrix3D> d2Nxdu2(m_basis.size());
    Matrix3D Hess;
    double   dXidu[3];
    Vec4     X;
    // Get element volume in the parameter space
    double du = el->umax() - el->umin();
    double dv = el->vmax() - el->vmin();
    double dw = el->wmax() - el->wmin();
    double vol = el->volume();
    if (vol < 0.0)
    {
      ok = false; // topology error (probably logic error)
      break;
    }

    // Set up control point (nodal) coordinates for current element
    if (!this->getElementCoordinates(Xnod,iEl+1))
    {
      ok = false;
      break;
    }

    // Compute parameter values of the Gauss points over the whole element
    std::array<RealArray,3> gpar, redpar;
    for (int d = 0; d < 3; d++)
    {
      this->getGaussPointParameters(gpar[d],d,nGauss,iEl,xg);
      if (xr)
        this->getGaussPointParameters(redpar[d],d,nRed,iEl,xr);
    }


    if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
      this->getElementCorners(iEl, fe.XC);

    if (integrand.getIntegrandType() & Integrand::G_MATRIX)
    {
      // Element size in parametric space
      dXidu[0] = el->getParmin(0);
      dXidu[1] = el->getParmin(1);
      dXidu[2] = el->getParmin(2);
    }
    else if (integrand.getIntegrandType() & Integrand::AVERAGE)
    {
      // --- Compute average value of basis functions over the element -----

      fe.Navg.resize(elem_sizes[0],true);
      double vol = 0.0;
      for (int k = 0; k < nGauss; k++)
        for (int j = 0; j < nGauss; j++)
          for (int i = 0; i < nGauss; i++)
          {
            // Fetch basis function derivatives at current integration point
            for (size_t b = 1; b <= m_basis.size(); ++b)
              evaluateBasis(fe, dNxdu[b-1], b);

            // Compute Jacobian determinant of coordinate mapping
            // and multiply by weight of current integration point
            double detJac = utl::Jacobian(Jac,fe.grad(geoBasis),
                                          Xnod,dNxdu[geoBasis-1],false);
            double weight = 0.125*vol*wg[i]*wg[j]*wg[k];

            // Numerical quadrature
            fe.Navg.add(fe.N,detJac*weight);
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
      this->getBasis(geoBasis)->point(X0,u0,v0,w0);
      X = SplineUtils::toVec3(X0);
    }

    // Initialize element quantities
    LocalIntegral* A = integrand.getLocalIntegral(elem_sizes,fe.iel);
    if (!integrand.initElement(MNPC[iEl],elem_sizes,nb,*A))
    {
      A->destruct();
      ok = false;
      break;
    }

    if (xr)
    {
      std::cerr << "Haven't really figured out what this part does yet\n";
      exit(42142);
    }

    // --- Integration loop over all Gauss points in each direction --------

    fe.iGP = iEl*nGauss*nGauss*nGauss; // Global integration point counter

    std::vector<Matrix> B(m_basis.size());
    size_t ig = 1;
    for (int k = 0; k < nGauss; k++)
      for (int j = 0; j < nGauss; j++)
        for (int i = 0; i < nGauss; i++, fe.iGP++, ig++)
        {
          // Local element coordinates of current integration point
          fe.xi   = xg[i];
          fe.eta  = xg[j];
          fe.zeta = xg[k];

          // Parameter values of current integration point
          fe.u = gpar[0][i];
          fe.v = gpar[1][j];
          fe.w = gpar[2][k];

          // Extract bezier basis functions
          for (size_t b = 0; b < m_basis.size(); ++b) {
            Matrix B(BN[b].rows(), 4);
            B.fillColumn(1, BN[b].getColumn(ig));
            B.fillColumn(2, BdNdu[b].getColumn(ig)*2.0/du);
            B.fillColumn(3, BdNdv[b].getColumn(ig)*2.0/dv);
            B.fillColumn(4, BdNdw[b].getColumn(ig)*2.0/dw);

            // Fetch basis function derivatives at current integration point
            if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
              evaluateBasis(fe, dNxdu[b], d2Nxdu2[b], b+1);
            else
              evaluateBasis(fe, dNxdu[b], bezierExtract[b][els[b]-1], B, b+1) ;
          }

          // Compute Jacobian inverse of coordinate mapping and derivatives
          fe.detJxW = utl::Jacobian(Jac,fe.grad(geoBasis),Xnod,dNxdu[geoBasis-1]);
          if (fe.detJxW == 0.0) continue; // skip singular points
          for (size_t b = 0; b < m_basis.size(); ++b)
            if (b != (size_t)geoBasis-1)
              fe.grad(b+1).multiply(dNxdu[b],Jac);

          // Compute Hessian of coordinate mapping and 2nd order derivatives
          if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES) {
            if (!utl::Hessian(Hess,fe.hess(geoBasis),Jac,Xnod,
                              d2Nxdu2[geoBasis-1],dNxdu[geoBasis-1]))
              return false;

            for (size_t b = 0; b < m_basis.size() && ok; ++b)
              if ((int)b != geoBasis)
                if (!utl::Hessian(Hess,fe.hess(b+1),Jac,Xnod,
                                  d2Nxdu2[b],fe.grad(b+1),false))
                  return false;
          }

          // Compute G-matrix
          if (integrand.getIntegrandType() & Integrand::G_MATRIX)
            utl::getGmat(Jac,dXidu,fe.G);

          // Cartesian coordinates of current integration point
          X   = Xnod * fe.basis(geoBasis);
          X.t = time.t;

          // Evaluate the integrand and accumulate element contributions
          fe.detJxW *= 0.125*vol*wg[i]*wg[j]*wg[k];
          if (!integrand.evalIntMx(*A,fe,time,X))
            ok = false;

    } // end gauss integrand

    // Finalize the element quantities
    if (ok && !integrand.finalizeElement(*A,time,0))
      ok = false;

    // Assembly of global system integral
    if (ok && !glInt.assemble(A->ref(),fe.iel))
      ok = false;

    A->destruct();
  }

  return ok;
}


bool ASMu3Dmx::integrate (Integrand& integrand, int lIndex,
                          GlobalIntegral& glInt,
                          const TimeDomain& time)
{
  PROFILE2("ASMu3Dmx::integrate(B)");

  // Get Gaussian quadrature points and weights
  int nGP = integrand.getBouIntegrationPoints(nGauss);
  const double* xg = GaussQuadrature::getCoord(nGP);
  const double* wg = GaussQuadrature::getWeight(nGP);
  if (!xg || !wg) return false;

  // Find the parametric direction of the face normal {-3,-2,-1, 1, 2, 3}
  const int faceDir = (lIndex+1)/(lIndex%2 ? -2 : 2);

  const int t1 = 1 + abs(faceDir)%3; // first tangent direction
  const int t2 = 1 + t1%3;           // second tangent direction

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;

  LR::parameterEdge edge;
  switch(lIndex)
  {
  case 1: edge = LR::WEST;   break;
  case 2: edge = LR::EAST;   break;
  case 3: edge = LR::SOUTH;  break;
  case 4: edge = LR::NORTH;  break;
  case 5: edge = LR::BOTTOM; break;
  case 6: edge = LR::TOP;    break;
  default:edge = LR::NONE;
  }

  // fetch all elements along the chosen edge
  std::vector<LR::Element*> edgeElms;
  this->getBasis(geoBasis)->getEdgeElements(edgeElms, (LR::parameterEdge) edge);

  // iterate over all edge elements
  bool ok = true;
  for(LR::Element *el : edgeElms) {
    double uh = (el->umin()+el->umax())/2.0;
    double vh = (el->vmin()+el->vmax())/2.0;
    double wh = (el->wmin()+el->wmax())/2.0;
    std::vector<size_t> els;
    std::vector<size_t> elem_sizes;
    for (size_t i=0; i < m_basis.size(); ++i) {
      els.push_back(m_basis[i]->getElementContaining(uh, vh, wh)+1);
      elem_sizes.push_back((*(m_basis[i]->elementBegin()+els.back()-1))->nBasisFunctions());
    }
    int iEl = el->getId();
    MxFiniteElement fe(elem_sizes);
    fe.iel = iEl+1;

    // Compute parameter values of the Gauss points over the whole element
    std::array<Vector,3> gpar;
    for (int d = 0; d < 3; d++)
      if (-1-d == faceDir)
      {
        gpar[d].resize(1);
        gpar[d].fill(this->getBasis(geoBasis)->startparam(d));
      }
      else if (1+d == faceDir)
      {
        gpar[d].resize(1);
        gpar[d].fill(this->getBasis(geoBasis)->endparam(d));
      }
      else
        this->getGaussPointParameters(gpar[d],d,nGP,iEl,xg);

    fe.xi = fe.eta = fe.zeta = faceDir < 0 ? -1.0 : 1.0;
    fe.u = gpar[0](1);
    fe.v = gpar[1](1);
    fe.w = gpar[2](1);

    std::vector<Matrix> dNxdu(m_basis.size());
    Matrix Xnod, Jac;
    Vec4   X;
    Vec3   normal;
    double dXidu[3];

    // Get element face area in the parameter space
    double dA = this->getParametricArea(iEl+1,abs(faceDir));
    if (dA < 0.0) // topology error (probably logic error)
    {
      ok = false;
      break;
    }

    // Set up control point coordinates for current element
    if (!this->getElementCoordinates(Xnod,iEl+1))
    {
      ok = false;
      break;
    }

    if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
      this->getElementCorners(iEl,fe.XC);

    if (integrand.getIntegrandType() & Integrand::G_MATRIX)
    {
      // Element size in parametric space
      dXidu[0] = el->getParmax(0) - el->getParmin(0);
      dXidu[1] = el->getParmax(1) - el->getParmin(1);
      dXidu[2] = el->getParmax(2) - el->getParmin(2);
    }

    // Initialize element quantities
    LocalIntegral* A = integrand.getLocalIntegral(elem_sizes,fe.iel,true);
    if (!integrand.initElementBou(MNPC[iEl],elem_sizes,nb,*A))
    {
      A->destruct();
      ok = false;
      break;
    }

    // --- Integration loop over all Gauss points in each direction --------

    fe.iGP = firstp; // Global integration point counter
    int k1,k2,k3;
    for (int j = 0; j < nGP; j++)
      for (int i = 0; i < nGP; i++, fe.iGP++)
      {
        // Local element coordinates and parameter values
        // of current integration point
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
          fe.u = gpar[0](k1+1);
        }
        if (gpar[1].size() > 1)
        {
          fe.eta = xg[k2];
          fe.v = gpar[1](k2+1);
        }
        if (gpar[2].size() > 1)
        {
          fe.zeta = xg[k3];
          fe.w = gpar[2](k3+1);
        }

        // Fetch basis function derivatives at current integration point
        for (size_t b = 1; b <= m_basis.size(); ++b)
          evaluateBasis(fe, dNxdu[b-1], b);

        // Compute basis function derivatives and the face normal
        fe.detJxW = utl::Jacobian(Jac, normal, fe.grad(geoBasis), Xnod,
                                  dNxdu[geoBasis-1], t1, t2);
        if (fe.detJxW == 0.0) continue; // skip singular points
        for (size_t b = 0; b < m_basis.size(); ++b)
          if (b != (size_t)geoBasis-1)
            fe.grad(b+1).multiply(dNxdu[b],Jac);

        if (faceDir < 0) normal *= -1.0;

        // Compute G-matrix
        if (integrand.getIntegrandType() & Integrand::G_MATRIX)
          utl::getGmat(Jac,dXidu,fe.G);

        // Cartesian coordinates of current integration point
        X = Xnod * fe.basis(geoBasis);
        X.t = time.t;

        // Evaluate the integrand and accumulate element contributions
        fe.detJxW *= 0.25*dA*wg[i]*wg[j];
        if (!integrand.evalBouMx(*A,fe,time,X,normal))
          ok = false;
    }

    // Finalize the element quantities
    if (ok && !integrand.finalizeElementBou(*A,fe,time))
      ok = false;

    // Assembly of global system integral
    if (ok && !glInt.assemble(A->ref(),fe.iel))
      ok = false;
    A->destruct();

    firstp += nGP*nGP*nGP;
  }

  return ok;
}


bool ASMu3Dmx::evalSolution (Matrix& sField, const Vector& locSol,
                             const RealArray* gpar, bool regular,
                             int deriv) const
{
  return false;
}


bool ASMu3Dmx::evalSolution (Matrix& sField, const IntegrandBase& integrand,
                             const RealArray* gpar, bool regular) const
{
  return evalSolution(sField,
                     const_cast<IntegrandBase&>(integrand).getSolution(0),
                     gpar, regular, 0);
}


bool ASMu3Dmx::refine (const LR::RefineData& prm,
                       Vectors& sol, const char* fName)
{
  return false;
}


Vec3 ASMu3Dmx::getCoord (size_t inod) const
{
  size_t b = 0;
  size_t nbb = 0;
  while (nbb + nb[b] < inod && b < nb.size())
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
