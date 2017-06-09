// $Id$
//==============================================================================
//!
//! \file ASMs2Dmx.C
//!
//! \date Oct 18 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of structured 2D spline mixed FE models.
//!
//==============================================================================

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SurfaceInterpolator.h"

#include "ASMs2Dmx.h"
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
#ifdef USE_OPENMP
#include <omp.h>
#endif


ASMs2Dmx::ASMs2Dmx (unsigned char n_s,
		    const std::vector<unsigned char>& n_f)
  : ASMs2D(n_s), ASMmxBase(n_f)
{
}


ASMs2Dmx::ASMs2Dmx (const ASMs2Dmx& patch,
                    const std::vector<unsigned char>& n_f)
  : ASMs2D(patch), ASMmxBase(n_f)
{
  nb = patch.nb;
}


Go::SplineSurface* ASMs2Dmx::getBasis (int basis) const
{
  if (basis < 1 || basis > (int)m_basis.size())
    return surf;

  return m_basis[basis-1].get();
}


Go::SplineCurve* ASMs2Dmx::getBoundary (int dir, int basis)
{
  if (dir < -2 || dir == 0 || dir > 2 || basis < 1 || basis > (int)m_basis.size())
    return nullptr;

  int iedge = dir > 0 ? dir : 3*dir+6;

  return m_basis[basis-1]->edgeCurve(iedge);
}


bool ASMs2Dmx::write (std::ostream& os, int basis) const
{
  if (basis < 1 || basis > (int)m_basis.size())
    os <<"200 1 0 0\n" << *surf;
  else if (m_basis[basis-1])
    os <<"200 1 0 0\n" << *m_basis[basis-1];
  else
    return false;

  return os.good();
}


void ASMs2Dmx::clear (bool retainGeometry)
{
  // Erase the spline data
  if (!retainGeometry)
    delete surf, surf = 0;

  for (auto& it : m_basis)
    it.reset();

  // Erase the FE data
  this->ASMs2D::clear(retainGeometry);
}


size_t ASMs2Dmx::getNoNodes (int basis) const
{
  if (basis < 1 || basis > (int)nb.size())
    return this->ASMbase::getNoNodes(basis);

  return nb[basis-1];
}


unsigned char ASMs2Dmx::getNoFields (int basis) const
{
  if (basis < 1 || basis > (int)nfx.size())
    return std::accumulate(nfx.begin(), nfx.end(), 0);

  return nfx[basis-1];
}


unsigned char ASMs2Dmx::getNodalDOFs (size_t inod) const
{
  if (this->isLMn(inod))
    return nLag;

  size_t nbc=0;
  for (size_t i=0;i<nb.size();++i)
    if (inod <= (nbc+=nb[i]))
      return nfx[i];

  return nfx[0];
}


char ASMs2Dmx::getNodeType (size_t inod) const
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


void ASMs2Dmx::initMADOF (const int* sysMadof)
{
  this->initMx(MLGN,sysMadof);
}


void ASMs2Dmx::extractNodeVec (const Vector& globRes, Vector& nodeVec,
			       unsigned char, int basis) const
{
  this->extractNodeVecMx(globRes,nodeVec,basis);
}


bool ASMs2Dmx::injectNodeVec (const Vector& nodeRes, Vector& globRes,
			      unsigned char, int basis) const
{
  this->injectNodeVecMx(globRes,nodeRes,basis);
  return true;
}


bool ASMs2Dmx::getSolution (Matrix& sField, const Vector& locSol,
			    const IntVec& nodes) const
{
  return this->getSolutionMx(sField,locSol,nodes);
}


bool ASMs2Dmx::generateFEMTopology ()
{
  if (!surf) return false;

  if (m_basis.empty())
    m_basis = ASMmxBase::establishBases(surf, ASMmxBase::Type);

  nb.clear();
  for (auto it : m_basis)
    nb.push_back(it->numCoefs_u()*it->numCoefs_v());

  if (!nodeInd.empty() && !shareFE)
  {
    if (nodeInd.size() == std::accumulate(nb.begin(), nb.end(), 0u)) return true;
    std::cerr <<" *** ASMs2Dmx::generateFEMTopology: Inconsistency between the"
	      <<" number of FE nodes "<< nodeInd.size()
	      <<"\n     and the number of spline coefficients "
              << std::accumulate(nb.begin(), nb.end(), 0)
	      <<" in the patch."<< std::endl;
    return false;
  }

  if (shareFE == 'F') return true;

  int i1, i2, j1, j2;
#ifdef SP_DEBUG
  size_t nbasis=0;
  for (auto& it : m_basis) {
    int i1, i2;
    std::cout << "Basis " << ++nbasis << ":\n";
    std::cout <<"numCoefs: "<< it->numCoefs_u() <<" "<< it->numCoefs_v();
    std::cout <<"\norder: "<< it->order_u() <<" "<< it->order_v();
    std::cout <<"\ndu:";
    for (i1 = 0; i1 < it->numCoefs_u(); i1++)
      std::cout <<' '<< it->knotSpan(0,i1);
    std::cout <<"\ndv:";
    for (i2 = 0; i2 < it->numCoefs_v(); i2++)
      std::cout <<' '<< it->knotSpan(1,i2);
    std::cout << std::endl;
  }
#endif

  nel = (m_basis[geoBasis-1]->numCoefs_u()-m_basis[geoBasis-1]->order_u()+1)*
        (m_basis[geoBasis-1]->numCoefs_v()-m_basis[geoBasis-1]->order_v()+1);
  nnod = std::accumulate(nb.begin(), nb.end(), 0);

  myMLGE.resize(nel,0);
  myMLGN.resize(nnod);
  myMNPC.resize(nel);
  myNodeInd.resize(nnod);

  size_t iel, inod = 0;
  for (auto& it : m_basis) {
    for (i2 = 0; i2 < it->numCoefs_v(); i2++)
      for (i1 = 0; i1 < it->numCoefs_u(); i1++)
      {
        myNodeInd[inod].I = i1;
        myNodeInd[inod].J = i2;
        myMLGN[inod++]    = ++gNod;
      }
  }

  int lnod2 = 0;
  int lnod3 = 0;

  auto&& addBasis = [&](std::shared_ptr<Go::SplineSurface>& b, bool geo) {
    for (i2 = 1; i2 <= b->numCoefs_v(); i2++)
      for (i1 = 1; i1 <= b->numCoefs_u(); i1++, inod++)
        if (i1 >= b->order_u() && i2 >= b->order_v())
        {
          if (b->knotSpan(0,i1-1) > 0.0)
            if (b->knotSpan(1,i2-1) > 0.0)
            {
              if (geo)
                myMLGE[iel] = ++gEl; // global element number over all patches
              else
                while (iel < myMNPC.size() && myMNPC[iel].empty()) iel++;

              int lnod = lnod2;
              myMNPC[iel].resize(lnod3,0);

              for (j2 = b->order_v()-1; j2 >= 0; j2--)
                for (j1 = b->order_u()-1; j1 >= 0; j1--)
                  myMNPC[iel][lnod++] = inod - b->numCoefs_u()*j2 - j1;
              if (!geo)
                ++iel;
            }

          if (geo)
            ++iel;
        }
  };

  iel = 0, inod = std::accumulate(nb.begin(),nb.begin()+geoBasis-1,0u);

  for (i2 = 0; i2 < geoBasis-1; ++i2)
    lnod2 += m_basis[i2]->order_u()*m_basis[i2]->order_v();
  for (i2 = 0; i2 < (int)m_basis.size(); ++i2)
    lnod3 += m_basis[i2]->order_u()*m_basis[i2]->order_v();

  // Create nodal connectivities for geometry basis
  addBasis(m_basis[geoBasis-1],true);
  // Create nodal connectivities for other bases
  inod = 0;
  lnod2 = 0;
  for (size_t b = 0; b < m_basis.size(); ++b) {
    iel = 0;
    if ((int)b != geoBasis-1)
      addBasis(m_basis[b],false);
    else
      inod += nb[b];
    lnod2 += m_basis[b]->order_u()*m_basis[b]->order_v();
  }

#ifdef SP_DEBUG
  std::cout <<"NEL = "<< nel <<" NNOD = "<< nnod << std::endl;
#endif
  return true;
}


bool ASMs2Dmx::connectPatch (int edge, ASM2D& neighbor, int nedge, bool revers,
                             int basis, bool coordCheck, int thick)
{
  ASMs2Dmx* neighMx = dynamic_cast<ASMs2Dmx*>(&neighbor);
  if (!neighMx) return false;

  size_t nb1 = 0, nb2 = 0;
  for (size_t i = 1; i <= m_basis.size(); ++i) {
    if (basis == 0 || i == (size_t)basis)
      if (!this->connectBasis(edge,*neighMx,nedge,revers,i,nb1,nb2,coordCheck,thick))
        return false;

    nb1 += nb[i-1];
    nb2 += neighMx->nb[i-1];
  }

  this->addNeighbor(neighMx);
  return true;
}


void ASMs2Dmx::closeEdges (int dir, int, int)
{
  size_t nbi = 1;
  for (size_t i = 1;i <= m_basis.size(); ++i) {
    this->ASMs2D::closeEdges(dir,i,nbi);
    nbi += nb[i-1];
  }
}


bool ASMs2Dmx::getElementCoordinates (Matrix& X, int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs2Dmx::getElementCoordinates: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return false;
  }
#endif

  size_t nenod = surf->order_u()*surf->order_v();
  size_t lnod0 = 0;

  for (int i = 1; i < geoBasis; ++i)
    lnod0 += m_basis[i-1]->order_u()*m_basis[i-1]->order_v();

  X.resize(nsd,nenod);
  const IntVec& mnpc = MNPC[iel-1];

  RealArray::const_iterator cit = surf->coefs_begin();
  for (size_t n = 0; n < nenod; n++)
  {
    int iI = nodeInd[mnpc[lnod0+n]].I;
    int iJ = nodeInd[mnpc[lnod0+n]].J;
    int ip = (iJ*surf->numCoefs_u() + iI)*surf->dimension();
    for (size_t i = 0; i < nsd; i++)
      X(i+1,n+1) = *(cit+(ip+i));
  }

#if SP_DEBUG > 2
  std::cout <<"\nCoordinates for element "<< iel << X << std::endl;
#endif
  return true;
}


Vec3 ASMs2Dmx::getCoord (size_t inod) const
{
  if (inod > nodeInd.size() && inod <= MLGN.size())
  {
    // This is a node added due to constraints in local directions.
    // Find the corresponding original node (see constrainEdgeLocal)
    std::map<size_t,size_t>::const_iterator it = xnMap.find(inod);
    if (it != xnMap.end()) inod = it->second;
  }

#ifdef INDEX_CHECK
  if (inod < 1 || inod > nodeInd.size())
  {
    std::cerr <<" *** ASMs2Dmx::getCoord: Node index "<< inod
              <<" out of range [1,"<< nodeInd.size() <<"]."<< std::endl;
    return Vec3();
  }
#endif

  RealArray::const_iterator cit;
  const int I = nodeInd[inod-1].I;
  const int J = nodeInd[inod-1].J;

  size_t b = 0;
  size_t nbb = nb[0];
  while (nbb < inod)
    nbb += nb[++b];

  cit = m_basis[b]->coefs_begin()
      + (J*m_basis[b]->numCoefs_u()+I) * m_basis[b]->dimension();

  Vec3 X;
  for (size_t i = 0; i < nsd; i++, cit++)
    X[i] = *cit;

  return X;
}


bool ASMs2Dmx::getSize (int& n1, int& n2, int basis) const
{
  if (basis == 0)
    return ASMs2D::getSize(n1,n2);

  if (basis < 1 || basis > (int)m_basis.size())
    return false;

  n1 = m_basis[basis-1]->numCoefs_u();
  n2 = m_basis[basis-1]->numCoefs_v();
  return true;
}


#define DERR -999.99

double ASMs2Dmx::getParametricArea (int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs2Dmx::getParametricArea: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return DERR;
  }
#endif
  if (MNPC[iel-1].empty())
    return 0.0;

  std::vector<size_t> elem_sizes;
  for (auto& it : m_basis)
    elem_sizes.push_back(it->order_u()*it->order_v());

  int inod1 = MNPC[iel-1][std::accumulate(elem_sizes.begin(),
                                          elem_sizes.begin()+geoBasis, -1)];
#ifdef INDEX_CHECK
  if (inod1 < 0 || (size_t)inod1 >= nnod)
  {
    std::cerr <<" *** ASMs2Dmx::getParametricArea: Node index "<< inod1
	      <<" out of range [0,"<< nnod <<">."<< std::endl;
    return DERR;
  }
#endif

  double du = m_basis[geoBasis-1]->knotSpan(0,nodeInd[inod1].I);
  double dv = m_basis[geoBasis-1]->knotSpan(1,nodeInd[inod1].J);

  return du*dv;
}


double ASMs2Dmx::getParametricLength (int iel, int dir) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs2Dmx::getParametricLength: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return DERR;
  }
#endif
  if (MNPC[iel-1].empty())
    return 0.0;

  std::vector<size_t> elem_sizes;
  for (auto& it : m_basis)
    elem_sizes.push_back(it->order_u()*it->order_v());

  int inod1 = MNPC[iel-1][std::accumulate(elem_sizes.begin(),
                                          elem_sizes.begin()+geoBasis, -1)];
#ifdef INDEX_CHECK
  if (inod1 < 0 || (size_t)inod1 >= nnod)
  {
    std::cerr <<" *** ASMs2Dmx::getParametricLength: Node index "<< inod1
	      <<" out of range [0,"<< nnod <<">."<< std::endl;
    return DERR;
  }
#endif

  switch (dir)
    {
    case 1: return m_basis[geoBasis-1]->knotSpan(0,nodeInd[inod1].I);
    case 2: return m_basis[geoBasis-1]->knotSpan(1,nodeInd[inod1].J);
    }

  std::cerr <<" *** ASMs2Dmx::getParametricLength: Invalid edge direction "
	    << dir << std::endl;
  return DERR;
}


bool ASMs2Dmx::integrate (Integrand& integrand,
			  GlobalIntegral& glInt,
			  const TimeDomain& time)
{
  if (!surf) return true; // silently ignore empty patches

  PROFILE2("ASMs2Dmx::integrate(I)");

  bool useElmVtx = integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS;

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  std::array<Matrix,2> gpar;
  for (int d = 0; d < 2; d++)
    this->getGaussPointParameters(gpar[d],d,nGauss,xg);

  // Evaluate basis function derivatives at all integration points
  std::vector<std::vector<Go::BasisDerivsSf>> splinex(m_basis.size());
  std::vector<std::vector<Go::BasisDerivsSf2>> splinex2(m_basis.size());
#pragma omp parallel for schedule(static)
  for (size_t i=0;i<m_basis.size();++i)
    if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
      m_basis[i]->computeBasisGrid(gpar[0],gpar[1],splinex2[i]);
    else
      m_basis[i]->computeBasisGrid(gpar[0],gpar[1],splinex[i]);

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int n1 = surf->numCoefs_u();
  const int nel1 = n1 - p1 + 1;

  std::vector<size_t> elem_sizes;
  for (auto& it : m_basis)
    elem_sizes.push_back(it->order_u()*it->order_v());

  // === Assembly loop over all elements in the patch ==========================

  bool ok=true;
  for (size_t g=0;g<threadGroups.size() && ok;++g) {
#pragma omp parallel for schedule(static)
    for (size_t t=0;t<threadGroups[g].size();++t) {
      MxFiniteElement fe(elem_sizes);
      std::vector<Matrix> dNxdu(m_basis.size());
      std::vector<Matrix3D> d2Nxdu2(m_basis.size());
      Matrix3D Hess;
      double dXidu[2];
      Matrix Xnod, Jac;
      Vec4   X;
      for (size_t i = 0; i < threadGroups[g][t].size() && ok; ++i)
      {
        int iel = threadGroups[g][t][i];
        fe.iel = MLGE[iel];
        if (fe.iel < 1) continue; // zero-area element

        int i1 = p1 + iel % nel1;
        int i2 = p2 + iel / nel1;

        // Get element area in the parameter space
        double dA = this->getParametricArea(++iel);
        if (dA < 0.0) // topology error (probably logic error)
        {
          ok = false;
          break;
        }

        // Set up control point (nodal) coordinates for current element
        if (!this->getElementCoordinates(Xnod,iel))
        {
          ok = false;
          break;
        }

        if (useElmVtx)
          this->getElementCorners(i1-1,i2-1,fe.XC);

        if (integrand.getIntegrandType() & Integrand::G_MATRIX)
        {
          // Element size in parametric space
          dXidu[0] = surf->knotSpan(0,i1-1);
          dXidu[1] = surf->knotSpan(1,i2-1);
        }

        // Initialize element quantities
        LocalIntegral* A = integrand.getLocalIntegral(elem_sizes,fe.iel,false);
        if (!integrand.initElement(MNPC[iel-1],elem_sizes,nb,*A))
        {
          A->destruct();
          ok = false;
          break;
        }


        // --- Integration loop over all Gauss points in each direction --------

        int ip = ((i2-p2)*nGauss*nel1 + i1-p1)*nGauss;
        int jp = ((i2-p2)*nel1 + i1-p1)*nGauss*nGauss;
        fe.iGP = firstIp + jp; // Global integration point counter

        for (int j = 0; j < nGauss; j++, ip += nGauss*(nel1-1))
          for (int i = 0; i < nGauss; i++, ip++, fe.iGP++)
          {
            // Local element coordinates of current integration point
            fe.xi  = xg[i];
            fe.eta = xg[j];

            // Parameter values of current integration point
            fe.u = gpar[0](i+1,i1-p1+1);
            fe.v = gpar[1](j+1,i2-p2+1);

            // Fetch basis function derivatives at current integration point
            for (size_t b = 0; b < m_basis.size(); ++b)
              if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
                SplineUtils::extractBasis(splinex2[b][ip],fe.basis(b+1),dNxdu[b],d2Nxdu2[b]);
              else
                SplineUtils::extractBasis(splinex[b][ip],fe.basis(b+1),dNxdu[b]);

            // Compute Jacobian inverse of the coordinate mapping and
            // basis function derivatives w.r.t. Cartesian coordinates
            fe.detJxW = utl::Jacobian(Jac,fe.grad(geoBasis),Xnod,
                                      dNxdu[geoBasis-1]);
            if (fe.detJxW == 0.0) continue; // skip singular points
            for (size_t b = 0; b < m_basis.size(); ++b)
              if (b != (size_t)geoBasis-1)
                fe.grad(b+1).multiply(dNxdu[b],Jac);

            // Compute Hessian of coordinate mapping and 2nd order derivatives
            if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES) {
              if (!utl::Hessian(Hess,fe.hess(geoBasis),Jac,Xnod,
                                d2Nxdu2[geoBasis-1],fe.grad(geoBasis),true))
                ok = false;
              for (size_t b = 0; b < m_basis.size() && ok; ++b)
                if ((int)b != geoBasis)
                  if (!utl::Hessian(Hess,fe.hess(b+1),Jac,Xnod,
                                    d2Nxdu2[b],fe.grad(b+1),false))
                    ok = false;
            }

            // Compute G-matrix
            if (integrand.getIntegrandType() & Integrand::G_MATRIX)
              utl::getGmat(Jac,dXidu,fe.G);

            // Cartesian coordinates of current integration point
            X = Xnod * fe.basis(geoBasis);
            X.t = time.t;

            // Evaluate the integrand and accumulate element contributions
            fe.detJxW *= 0.25*dA*wg[i]*wg[j];
            if (!integrand.evalIntMx(*A,fe,time,X))
              ok = false;
          }

        // Finalize the element quantities
        if (ok && !integrand.finalizeElement(*A,time,firstIp+jp))
          ok = false;

        // Assembly of global system integral
        if (ok && !glInt.assemble(A->ref(),fe.iel))
          ok = false;

        A->destruct();
      }
    }
  }

  return ok;
}


bool ASMs2Dmx::integrate (Integrand& integrand, int lIndex,
			  GlobalIntegral& glInt,
			  const TimeDomain& time)
{
  if (!surf) return true; // silently ignore empty patches

  PROFILE2("ASMs2Dmx::integrate(B)");

  bool useElmVtx = integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS;

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Find the parametric direction of the edge normal {-2,-1, 1, 2}
  const int edgeDir = (lIndex%10+1)/(lIndex%2 ? -2 : 2);

  const int t1 = abs(edgeDir);   // Tangent direction normal to the patch edge
  const int t2 = 3-abs(edgeDir); // Tangent direction along the patch edge

  // Compute parameter values of the Gauss points along the whole patch edge
  std::array<Matrix,2> gpar;
  for (short int d = 0; d < 2; d++)
    if (-1-d == edgeDir)
    {
      gpar[d].resize(1,1);
      gpar[d].fill(d == 0 ? surf->startparam_u() : surf->startparam_v());
    }
    else if (1+d == edgeDir)
    {
      gpar[d].resize(1,1);
      gpar[d].fill(d == 0 ? surf->endparam_u() : surf->endparam_v());
    }
    else
      this->getGaussPointParameters(gpar[d],d,nGauss,xg);

  // Extract the Neumann order flag (1 or higher) for the integrand
  integrand.setNeumannOrder(1 + lIndex/10);

  // Evaluate basis function derivatives at all integration points
  std::vector<std::vector<Go::BasisDerivsSf>> splinex(m_basis.size());
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < m_basis.size(); ++i)
    m_basis[i]->computeBasisGrid(gpar[0],gpar[1],splinex[i]);

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();

  std::vector<size_t> elem_sizes;
  for (auto& it : m_basis)
    elem_sizes.push_back(it->order_u()*it->order_v());

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex%10);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;

  MxFiniteElement fe(elem_sizes);
  fe.xi = fe.eta = edgeDir < 0 ? -1.0 : 1.0;
  fe.u = gpar[0](1,1);
  fe.v = gpar[1](1,1);

  std::vector<Matrix> dNxdu(m_basis.size());
  Matrix Xnod, Jac;
  Vec4   X;
  Vec3   normal;


  // === Assembly loop over all elements on the patch edge =====================

  int iel = 1;
  for (int i2 = p2; i2 <= n2; i2++)
    for (int i1 = p1; i1 <= n1; i1++, iel++)
    {
      fe.iel = MLGE[iel-1];
      if (fe.iel < 1) continue; // zero-area element

      // Skip elements that are not on current boundary edge
      bool skipMe = false;
      switch (edgeDir)
	{
	case -1: if (i1 > p1) skipMe = true; break;
	case  1: if (i1 < n1) skipMe = true; break;
	case -2: if (i2 > p2) skipMe = true; break;
	case  2: if (i2 < n2) skipMe = true; break;
	}
      if (skipMe) continue;

      // Get element edge length in the parameter space
      double dS = this->getParametricLength(iel,t2);
      if (dS < 0.0) return false; // topology error (probably logic error)

      // Set up control point coordinates for current element
      if (!this->getElementCoordinates(Xnod,iel)) return false;

      if (useElmVtx)
        this->getElementCorners(i1-1,i2-1,fe.XC);

      // Initialize element quantities
      LocalIntegral* A = integrand.getLocalIntegral(elem_sizes,fe.iel,true);
      bool ok = integrand.initElementBou(MNPC[iel-1],elem_sizes,nb,*A);


      // --- Integration loop over all Gauss points along the edge -------------

      int ip = (t1 == 1 ? i2-p2 : i1-p1)*nGauss;
      fe.iGP = firstp + ip; // Global integration point counter

      for (int i = 0; i < nGauss && ok; i++, ip++, fe.iGP++)
      {
	// Parameter values of current integration point
	if (gpar[0].size() > 1)
	{
          fe.xi = xg[i];
	  fe.u = gpar[0](i+1,i1-p1+1);
	}
	if (gpar[1].size() > 1)
	{
          fe.eta = xg[i];
	  fe.v = gpar[1](i+1,i2-p2+1);
	}

	// Fetch basis function derivatives at current integration point
        for (size_t b = 0; b < m_basis.size(); ++b)
          SplineUtils::extractBasis(splinex[b][ip],fe.basis(b+1),dNxdu[b]);

	// Compute Jacobian inverse of the coordinate mapping and
	// basis function derivatives w.r.t. Cartesian coordinates
        fe.detJxW = utl::Jacobian(Jac,normal,fe.grad(geoBasis),Xnod,dNxdu[geoBasis-1],t1,t2);
	if (fe.detJxW == 0.0) continue; // skip singular points
        for (size_t b = 0; b < m_basis.size(); ++b)
          if (b != (size_t)geoBasis-1)
            fe.grad(b+1).multiply(dNxdu[b],Jac);

	if (edgeDir < 0) normal *= -1.0;

	// Cartesian coordinates of current integration point
	X = Xnod * fe.basis(geoBasis);
	X.t = time.t;

	// Evaluate the integrand and accumulate element contributions
	fe.detJxW *= 0.5*dS*wg[i];

	ok = integrand.evalBouMx(*A,fe,time,X,normal);
      }

      // Finalize the element quantities
      if (ok && !integrand.finalizeElementBou(*A,fe,time))
        ok = false;

      // Assembly of global system integral
      if (ok && !glInt.assemble(A->ref(),fe.iel))
	return false;

      A->destruct();

      if (!ok) return false;
    }

  return true;
}


int ASMs2Dmx::evalPoint (const double* xi, double* param, Vec3& X) const
{
  if (!surf) return -2;

  param[0] = (1.0-xi[0])*surf->startparam_u() + xi[0]*surf->endparam_u();
  param[1] = (1.0-xi[1])*surf->startparam_v() + xi[1]*surf->endparam_v();

  Go::Point X0;
  surf->point(X0,param[0],param[1]);
  for (unsigned char d = 0; d < nsd; d++)
    X[d] = X0[d];

  // Check if this point matches any of the control points (nodes)
  return this->searchCtrlPt(m_basis[0]->coefs_begin(),m_basis[0]->coefs_end(),
                            X,m_basis[0]->dimension());
}


bool ASMs2Dmx::evalSolution (Matrix& sField, const Vector& locSol,
                             const RealArray* gpar, bool regular, int) const
{
  // Evaluate the basis functions at all points
  std::vector<std::vector<Go::BasisPtsSf>> splinex(m_basis.size());
  if (regular)
  {
    for (size_t b = 0; b < m_basis.size(); ++b)
      m_basis[b]->computeBasisGrid(gpar[0],gpar[1],splinex[b]);
  }
  else if (gpar[0].size() == gpar[1].size())
  {
    for (size_t b = 0; b < m_basis.size(); ++b) {
      splinex[b].resize(gpar[0].size());
      for (size_t i = 0; i < splinex[b].size(); i++)
        m_basis[b]->computeBasis(gpar[0][i],gpar[1][i],splinex[b][i]);
    }
  }
  else
    return false;

  std::vector<size_t> nc(nfx.size());
  std::copy(nfx.begin(), nfx.end(), nc.begin());

  // assume first basis only
  if (locSol.size() < std::inner_product(nb.begin(), nb.end(), nfx.begin(), 0u)) {
    std::fill(nc.begin(), nc.end(), 0);
    nc[0] = nfx[0];
  }

  if (std::inner_product(nb.begin(), nb.end(), nc.begin(), 0u) != locSol.size())
    return false;

  Matrix Xtmp;
  Vector Ytmp, Ztmp;

  // Evaluate the primary solution field at each point
  size_t nPoints = splinex[0].size();
  sField.resize(std::accumulate(nc.begin(), nc.end(), 0),nPoints);
  for (size_t i = 0; i < nPoints; i++)
  {
    size_t comp=0;
    for (size_t b = 0; b < m_basis.size(); ++b) {
      if (nc[b] == 0)
        continue;
      IntVec ip;
      scatterInd(m_basis[b]->numCoefs_u(),m_basis[b]->numCoefs_v(),
                 m_basis[b]->order_u(),m_basis[b]->order_v(),splinex[b][i].left_idx,ip);

      utl::gather(ip,nc[b],locSol,Xtmp,comp);
      if (b == 0)
        Xtmp.multiply(splinex[b][i].basisValues,Ytmp);
      else {
        Xtmp.multiply(splinex[b][i].basisValues,Ztmp);
        Ytmp.insert(Ytmp.end(),Ztmp.begin(),Ztmp.end());
      }
      comp += nc[b]*nb[b];
    }
    sField.fillColumn(1+i,Ytmp);
  }

  return true;
}


bool ASMs2Dmx::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			     const RealArray* gpar, bool regular) const
{
  sField.resize(0,0);

  // Evaluate the basis functions and their derivatives at all points
  std::vector<std::vector<Go::BasisDerivsSf>> splinex(m_basis.size());
  if (regular)
  {
    for (size_t b = 0; b < m_basis.size(); ++b)
      m_basis[b]->computeBasisGrid(gpar[0],gpar[1],splinex[b]);
  }
  else if (gpar[0].size() == gpar[1].size())
  {
    for (size_t b = 0; b < m_basis.size(); ++b) {
      splinex[b].resize(gpar[0].size());
      for (size_t i = 0; i < splinex[b].size(); i++)
        m_basis[b]->computeBasis(gpar[0][i],gpar[1][i],splinex[b][i]);
    }
  }

  std::vector<size_t> elem_sizes;
  for (auto& it : m_basis)
    elem_sizes.push_back(it->order_u()*it->order_v());

  // Fetch nodal (control point) coordinates
  Matrix Xnod, Xtmp;
  this->getNodalCoordinates(Xnod);

  MxFiniteElement fe(elem_sizes,firstIp);
  Vector          solPt;
  std::vector<Matrix> dNxdu(m_basis.size());
  Matrix          Jac;
  Vec3            X;

  // Evaluate the secondary solution field at each point
  size_t nPoints = splinex[0].size();
  for (size_t i = 0; i < nPoints; i++, fe.iGP++)
  {
    // Fetch indices of the non-zero basis functions at this point
    std::vector<IntVec> ip(m_basis.size());
    IntVec ipa;
    size_t ofs = 0;
    for (size_t b = 0; b < m_basis.size(); ++b) {
      scatterInd(m_basis[b]->numCoefs_u(),m_basis[b]->numCoefs_v(),
                 m_basis[b]->order_u(),m_basis[b]->order_v(),splinex[b][i].left_idx,ip[b]);

      // Fetch associated control point coordinates
      if (b == (size_t)geoBasis-1)
        utl::gather(ip[geoBasis-1], nsd, Xnod, Xtmp);

      for (auto& it : ip[b])
        it += ofs;
      ipa.insert(ipa.end(), ip[b].begin(), ip[b].end());
      ofs += nb[b];
    }

    fe.u = splinex[0][i].param[0];
    fe.v = splinex[0][i].param[1];

    // Fetch basis function derivatives at current integration point
    for (size_t b = 0; b < m_basis.size(); ++b)
      SplineUtils::extractBasis(splinex[b][i],fe.basis(b+1),dNxdu[b]);

    // Compute Jacobian inverse of the coordinate mapping and
    // basis function derivatives w.r.t. Cartesian coordinates
    fe.detJxW = utl::Jacobian(Jac,fe.grad(geoBasis),Xtmp,dNxdu[geoBasis-1]);

    for (size_t b = 1; b <= m_basis.size(); b++)
      if (b != (size_t)geoBasis)
      {
        if (fe.detJxW == 0.0)
          fe.grad(b).clear();
        else
          fe.grad(b).multiply(dNxdu[b-1],Jac);
      }

    // Cartesian coordinates of current integration point
    X = Xtmp * fe.basis(geoBasis);

    // Now evaluate the solution field
    if (!integrand.evalSol(solPt,fe,X,ipa,elem_sizes,nb))
      return false;
    else if (sField.empty())
      sField.resize(solPt.size(),nPoints,true);

    sField.fillColumn(1+i,solPt);
  }

  return true;
}


void ASMs2Dmx::generateThreadGroups (const Integrand& integrand, bool silence,
                                     bool ignoreGlobalLM)
{
  int p1 = 0, p2 = 0;
  for (const auto& it : m_basis) {
    if (it->order_u() > p1)
      p1 = it->order_u();
    if (it->order_v() > p2)
      p2 = it->order_v();
  }

  ASMs2D::generateThreadGroups(p1-1, p2-1, silence, ignoreGlobalLM);
}


void ASMs2Dmx::getBoundaryNodes (int lIndex, IntVec& nodes, int basis,
                                 int thick, int, bool local) const
{
  if (basis > 0)
    this->ASMs2D::getBoundaryNodes(lIndex, nodes, basis, thick, 0, local);
  else
    for (size_t b = 1; b <= this->getNoBasis(); ++b)
      this->ASMs2D::getBoundaryNodes(lIndex, nodes, b, thick, 0, local);
}
