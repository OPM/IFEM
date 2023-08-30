// $Id$
//==============================================================================
//!
//! \file ASMs3Dmx.C
//!
//! \date Dec 28 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of structured 3D spline mixed FE models.
//!
//==============================================================================

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/trivariate/SplineVolume.h"

#include "ASMs3Dmx.h"
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Point.h"
#include "Profiler.h"
#include <array>
#include <numeric>
#include <utility>
#ifdef USE_OPENMP
#include <omp.h>
#endif


ASMs3Dmx::ASMs3Dmx (const CharVec& n_f)
  : ASMs3D(std::accumulate(n_f.begin(),n_f.end(),0)), ASMmxBase(n_f)
{
}


ASMs3Dmx::ASMs3Dmx (const ASMs3Dmx& patch, const CharVec& n_f)
  : ASMs3D(patch,std::accumulate(n_f.begin(),n_f.end(),0)), ASMmxBase(n_f),
    m_basis(patch.m_basis)
{
  nb = patch.nb;
}


ASMs3Dmx::~ASMs3Dmx ()
{
  // these are managed by shared ptrs, make sure base class do not delete them.
  if (!this->separateProjectionBasis())
    projB = nullptr;
  geomB = svol = nullptr;
}


const Go::SplineVolume* ASMs3Dmx::getBasis (int basis) const
{
  if (basis < 1)
    return this->ASMs3D::getBasis(basis);

  if (basis > static_cast<int>(m_basis.size()))
    return nullptr;

  return m_basis[basis-1].get();
}


Go::SplineVolume* ASMs3Dmx::getBasis (int basis)
{
  return const_cast<Go::SplineVolume*>(std::as_const(*this).getBasis(basis));
}


Go::SplineSurface* ASMs3Dmx::getBoundary (int dir, int basis)
{
  if (dir < -3 || dir == 0 || dir > 3)
    return nullptr;
  else if (basis < 1 || basis > (int)m_basis.size())
    return nullptr;

  // The boundary surfaces are stored internally in the SplineVolume object
  int iface = dir > 0 ? 2*dir-1 : -2*dir-2;
  return m_basis[basis-1]->getBoundarySurface(iface).get();
}


bool ASMs3Dmx::readBasis (std::istream& is, size_t basis)
{
  if (basis < 1 || basis > nfx.size())
    return false;

  if (m_basis.empty())
  {
    m_basis.resize(nfx.size());
    nb.resize(nfx.size(), 0);
  }

  Go::ObjectHeader head;
  m_basis[--basis] = std::make_shared<Go::SplineVolume>();
  is >> head >> *m_basis[basis];
  nb[basis] = m_basis[basis]->numCoefs(0) *
              m_basis[basis]->numCoefs(1) *
              m_basis[basis]->numCoefs(2);

  return true;
}


void ASMs3Dmx::clear (bool retainGeometry)
{
  // these are managed by shared ptrs, make sure base class do not delete them.
  if (!this->separateProjectionBasis())
    projB = nullptr;
  geomB = svol = nullptr;

  // Erase the solution field bases
  for (auto& it : m_basis)
    it.reset();

  // Erase the FE data and the geometry basis
  this->ASMs3D::clear(retainGeometry);
}


size_t ASMs3Dmx::getNoNodes (int basis) const
{
  if (basis < 1 || basis > (int)nb.size())
    return this->ASMbase::getNoNodes(basis);

  return nb[basis-1];
}


unsigned char ASMs3Dmx::getNoFields (int basis) const
{
  return basis < 1 || basis > (int)nfx.size() ? nf : nfx[basis-1];
}


unsigned char ASMs3Dmx::getNodalDOFs (size_t inod) const
{
  if (this->isLMn(inod))
    return nLag;

  size_t nbc = 0;
  for (size_t i = 0; i < nb.size(); i++)
    if (inod <= (nbc += nb[i]))
      return nfx[i];

  return nfx.front();
}


char ASMs3Dmx::getNodeType (size_t inod) const
{
  if (this->isLMn(inod))
    return 'L';

  size_t nbc = nb.front();
  if (inod <= nbc)
    return 'D';
  else for (size_t i = 1; i < nb.size(); i++)
    if (inod <= (nbc += nb[i]))
      return 'O'+i;

  return 'X';
}


void ASMs3Dmx::initMADOF (const int* sysMadof)
{
  this->initMx(MLGN,sysMadof);
}


void ASMs3Dmx::extractNodeVec (const RealArray& globRes, RealArray& nodeVec,
                               unsigned char, int basis) const
{
  this->extractNodeVecMx(globRes,nodeVec,basis);
}


bool ASMs3Dmx::injectNodeVec (const RealArray& nodeRes, RealArray& globRes,
                              unsigned char, int basis) const
{
  this->injectNodeVecMx(globRes,nodeRes,basis);
  return true;
}


bool ASMs3Dmx::getSolution (Matrix& sField, const Vector& locSol,
			    const IntVec& nodes) const
{
  return this->getSolutionMx(sField,locSol,nodes);
}


bool ASMs3Dmx::generateFEMTopology ()
{
  if (!svol) return false;

  if (m_basis.empty()) {
    m_basis = ASMmxBase::establishBases(svol, ASMmxBase::Type);

    // we need to project on something that is not one of our bases
    if (!projB) {
      if (ASMmxBase::Type == ASMmxBase::REDUCED_CONT_RAISE_BASIS1 ||
          ASMmxBase::Type == ASMmxBase::REDUCED_CONT_RAISE_BASIS2 ||
          ASMmxBase::Type == ASMmxBase::DIV_COMPATIBLE)
        projB = ASMmxBase::raiseBasis(svol);
      else if (ASMmxBase::Type == ASMmxBase::SUBGRID)
        projB = m_basis.front().get();
      else // FULL_CONT_RAISE_BASISx
        projB = m_basis[2-geoBasis].get();
    }

    if (ASMmxBase::Type == ASMmxBase::SUBGRID)
      projB2 = ASMmxBase::raiseBasis(svol);

    delete svol;
  }
  geomB = svol = m_basis[geoBasis-1].get();

  nb.clear();
  nb.reserve(m_basis.size());
  elem_size.clear();
  elem_size.reserve(m_basis.size());
  for (auto& it : m_basis) {
    nb.push_back(it->numCoefs(0)*it->numCoefs(1)*it->numCoefs(2));
    elem_size.push_back(it->order(0)*it->order(1)*it->order(2));
  }

  nnod = std::accumulate(nb.begin(), nb.end(), 0u);
  if (!nodeInd.empty() && !shareFE)
  {
    if (nodeInd.size() == nnod)
      return true;

    std::cerr <<" *** ASMs3Dmx::generateFEMTopology: Inconsistency between the"
	      <<" number of FE nodes "<< nodeInd.size()
	      <<"\n     and the number of spline coefficients "<< nnod
	      <<" in the patch."<< std::endl;
    return false;
  }
  else if (shareFE == 'F')
    return true;

#ifdef SP_DEBUG
  size_t nbasis = 0;
  for (auto it : m_basis) {
    std::cout <<"Basis "<< ++nbasis;
    std::cout <<":\nnumCoefs: "<< it->numCoefs(0) <<" "
                               << it->numCoefs(1) <<" "
                               << it->numCoefs(2);
    std::cout <<"\norder: "<< it->order(0) <<" "
                           << it->order(1) <<" "
                           << it->order(2);
    for (int d = 0; d < 3; d++)
    {
      std::cout <<"\nd"<< char('u'+d) <<':';
      for (int i = 0; i < it->numCoefs(d); i++)
        std::cout <<' '<< it->knotSpan(d,i);
    }
    std::cout << std::endl;
  }
#endif

  const Go::SplineVolume* itg = this->getBasis(ASM::INTEGRATION_BASIS);

  nel = (itg->numCoefs(0) - itg->order(0) + 1) *
        (itg->numCoefs(1) - itg->order(1) + 1) *
        (itg->numCoefs(2) - itg->order(2) + 1);

  myMLGE.resize(nel,0);
  myMLGN.resize(nnod);
  myMNPC.resize(nel);
  myNodeInd.resize(nnod);

  size_t iel = 0, inod = 0;
  for (auto& it : m_basis)
    for (int k = 0; k < it->numCoefs(2); k++)
      for (int j = 0; j < it->numCoefs(1); j++)
        for (int i = 0; i < it->numCoefs(0); i++)
        {
          myNodeInd[inod].I = i;
          myNodeInd[inod].J = j;
          myNodeInd[inod].K = k;
          myMLGN[inod++]    = ++gNod;
        }

  const int firstItgElmNode = this->getFirstItgElmNode();
  const size_t nElmNode = std::accumulate(elem_size.begin(), elem_size.end(), 0);

  // Create nodal connectivities for bases
  auto&& enumerate_elem = [this,&iel](const Go::SplineVolume& spline,
                                      size_t elm_node, int last_node)
  {
    const int nu = spline.numCoefs(0);
    const int nv = spline.numCoefs(1);
    for (int k = spline.order(2) - 1; k >= 0; k--)
      for (int j = spline.order(1) - 1; j >= 0; j--)
        for (int i = spline.order(0) - 1; i >= 0; i--)
          myMNPC[iel][elm_node++] = last_node - (k*nv + j)*nu - i;
  };

  inod = std::accumulate(nb.begin(),nb.begin()+geoBasis-1,0u);
  auto knotw = itg->basis(2).begin();
  for (int k = 1; k <= itg->numCoefs(2); k++, ++knotw) {
    auto knotv = itg->basis(1).begin();
    for (int j = 1; j <= itg->numCoefs(1); j++, ++knotv) {
      auto knotu = itg->basis(0).begin();
      for (int i = 1; i <= itg->numCoefs(0); i++, inod++, ++knotu)
        if (i >= itg->order(0) && j >= itg->order(1) && k >= itg->order(2)) {
          if (itg->knotSpan(0,i-1) > 0.0 &&
              itg->knotSpan(1,j-1) > 0.0 &&
              itg->knotSpan(2,k-1) > 0.0) {
            myMLGE[iel] = ++gEl; // global element number over all patches

            myMNPC[iel].resize(nElmNode, 0);
            enumerate_elem(*itg, firstItgElmNode, inod);

            // find knot span for other bases
            size_t elem_ofs = 0;
            size_t basis_ofs = 0;
            int basis = 0;
            for (const std::shared_ptr<Go::SplineVolume>& spline : m_basis) {
              if (basis != geoBasis-1) {
                double ku = *knotu;
                double kv = *knotv;
                double kw = *knotw;
                int iu = spline->basis(0).knotIntervalFuzzy(ku);
                int iv = spline->basis(1).knotIntervalFuzzy(kv);
                int iw = spline->basis(2).knotIntervalFuzzy(kw);
                int iinod = basis_ofs + (iw*spline->numCoefs(1) + iv)*spline->numCoefs(0) + iu;
                enumerate_elem(*spline, elem_ofs, iinod);
              }
              elem_ofs += elem_size[basis];
              basis_ofs += nb[basis++];
            }
          }

          ++iel;
        }
    }
  }

#ifdef SP_DEBUG
  std::cout <<"NEL = "<< nel <<" NNOD = "<< nnod << std::endl;
#endif
  return true;
}


bool ASMs3Dmx::connectPatch (int face, ASM3D& neighbor, int nface, int norient,
                             int basis, bool coordCheck, int thick)
{
  ASMs3Dmx* neighMx = dynamic_cast<ASMs3Dmx*>(&neighbor);
  if (!neighMx) return false;

  if (swapW && face > 4) // Account for swapped parameter direction
    face = 11-face;

  if (neighMx->swapW && face > 4) // Account for swapped parameter direction
    nface = 11-nface;

  size_t nb1 = 0, nb2 = 0;
  for (size_t i = 1; i <= nb.size(); i++) {
    if (basis == 0 || i == (size_t)basis)
      if (!this->connectBasis(face,*neighMx,nface,norient,i,nb1,nb2,
                              coordCheck,thick))
        return false;

    nb1 += nb[i-1];
    nb2 += neighMx->nb[i-1];
  }

  this->addNeighbor(neighMx);
  return true;
}


void ASMs3Dmx::closeBoundaries (int dir, int, int)
{
  size_t nbi = 1;
  for (size_t i = 0; i < nb.size(); nbi += nb[i++])
    this->ASMs3D::closeBoundaries(dir,1+i,nbi);
}


bool ASMs3Dmx::getElementCoordinates (Matrix& X, int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs3Dmx::getElementCoordinates: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return false;
  }
#endif

  size_t nenod = svol->order(0)*svol->order(1)*svol->order(2);
  size_t lnod0 = 0;
  for (int i = 1; i < geoBasis; ++i)
    lnod0 += m_basis[i-1]->order(0)*m_basis[i-1]->order(1)*m_basis[i-1]->order(2);

  X.resize(3,nenod);
  const IntVec& mnpc = MNPC[iel-1];

  RealArray::const_iterator cit = svol->coefs_begin();
  for (size_t n = 0; n < nenod; n++)
  {
    int iI = nodeInd[mnpc[lnod0+n]].I;
    int iJ = nodeInd[mnpc[lnod0+n]].J;
    int iK = nodeInd[mnpc[lnod0+n]].K;
    int ip = ((iK*svol->numCoefs(1)+ iJ)*svol->numCoefs(0)+ iI)*svol->dimension();
    for (size_t i = 0; i < 3; i++)
      X(i+1,n+1) = *(cit+(ip+i));
  }

#if SP_DEBUG > 2
  std::cout <<"\nCoordinates for element "<< iel << X << std::endl;
#endif
  return true;
}


Vec3 ASMs3Dmx::getCoord (size_t inod) const
{
  if (inod > nodeInd.size() && inod <= MLGN.size())
  {
    // This is a node added due to constraints in local directions.
    // Find the corresponding original node (see constrainFaceLocal)
    std::map<size_t,size_t>::const_iterator it = xnMap.find(inod);
    if (it != xnMap.end()) inod = it->second;
  }

#ifdef INDEX_CHECK
  if (inod < 1 || inod > nodeInd.size())
  {
    std::cerr <<" *** ASMs3Dmx::getCoord: Node index "<< inod
              <<" out of range [1,"<< nodeInd.size() <<"]."<< std::endl;
    return Vec3();
  }
#endif

  RealArray::const_iterator cit;
  const int I = nodeInd[inod-1].I;
  const int J = nodeInd[inod-1].J;
  const int K = nodeInd[inod-1].K;

  size_t b = 0;
  size_t nbb = nb.front();
  while (nbb < inod)
    nbb += nb[++b];

  cit = m_basis[b]->coefs_begin()
      + ((K*m_basis[b]->numCoefs(1)+J)*m_basis[b]->numCoefs(0)+I) * m_basis[b]->dimension();

  return RealArray(cit,cit+3);
}


bool ASMs3Dmx::getSize (int& n1, int& n2, int& n3, int basis) const
{
  if (basis == 0)
    return this->ASMs3D::getSize(n1,n2,n3);

  if (basis < 1 || basis > (int)m_basis.size())
    return false;

  n1 = m_basis[basis-1]->numCoefs(0);
  n2 = m_basis[basis-1]->numCoefs(1);
  n3 = m_basis[basis-1]->numCoefs(2);
  return true;
}


bool ASMs3Dmx::integrate (Integrand& integrand,
			  GlobalIntegral& glInt,
			  const TimeDomain& time)
{
  if (!svol) return true; // silently ignore empty patches

  PROFILE2("ASMs3Dmx::integrate(I)");

  bool use2ndDer = integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES;
  bool useElmVtx = integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS;

  if (myCache.empty()) {
    myCache.emplace_back(std::make_unique<BasisFunctionCache>(*this, cachePolicy, 1));
    for (size_t b = 2; b <= this->getNoBasis(); ++b)
      myCache.emplace_back(std::make_unique<BasisFunctionCache>(*myCache.front(), b));
  }

  for (std::unique_ptr<BasisFunctionCache>& cache : myCache) {
    cache->init(use2ndDer ? 2 : 1);
  }

  BasisFunctionCache& cache = *myCache.front();

  const std::array<int,3>& ng = cache.nGauss();
  const std::array<const double*,3>& xg = cache.coord();
  const std::array<const double*,3>& wg = cache.weight();

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;

  ThreadGroups oneGroup;
  if (glInt.threadSafe()) oneGroup.oneStripe(nel);
  const ThreadGroups& groups = glInt.threadSafe() ? oneGroup : threadGroupsVol;


  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t g = 0; g < groups.size() && ok; g++)
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < groups[g].size(); t++)
    {
      MxFiniteElement fe(elem_size);
      Matrix3D Hess;
      double dXidu[3];
      Matrix Xnod, Jac;
      double param[3] = { 0.0, 0.0, 0.0 };
      Vec4   X(param,time.t);
      for (size_t l = 0; l < groups[g][t].size() && ok; l++)
      {
        int iel = groups[g][t][l];
        fe.iel = MLGE[iel];
        if (fe.iel < 1) continue; // zero-volume element

        int i1 = p1 + iel % nel1;
        int i2 = p2 + (iel / nel1) % nel2;
        int i3 = p3 + iel / (nel1*nel2);

        // Get element volume in the parameter space
        double dV = 0.125*this->getParametricVolume(++iel);
        if (dV < 0.0) // topology error (probably logic error)
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
          fe.h = this->getElementCorners(i1-1,i2-1,i3-1,fe.XC);

        if (integrand.getIntegrandType() & Integrand::G_MATRIX)
        {
          // Element size in parametric space
          dXidu[0] = svol->knotSpan(0,i1-1);
          dXidu[1] = svol->knotSpan(1,i2-1);
          dXidu[2] = svol->knotSpan(2,i3-1);
        }

        // Initialize element quantities
        LocalIntegral* A = integrand.getLocalIntegral(elem_size,fe.iel,false);
        if (!integrand.initElement(MNPC[iel-1],fe,elem_size,nb,*A))
        {
          A->destruct();
          ok = false;
          break;
        }


        // --- Integration loop over all Gauss points in each direction --------

        size_t ip = 0;
        int jp = (((i3-p3)*nel2 + i2-p2)*nel1 + i1-p1)*ng[0]*ng[1]*ng[2];
        fe.iGP = firstIp + jp; // Global integration point counter

        for (int k = 0; k < ng[2]; k++)
          for (int j = 0; j < ng[1]; j++)
            for (int i = 0; i < ng[0]; i++, ip++, fe.iGP++)
            {
              // Local element coordinates of current integration point
              fe.xi   = xg[0][i];
              fe.eta  = xg[1][j];
              fe.zeta = xg[2][k];

              // Parameter values of current integration point
              fe.u = param[0] = cache.getParam(0,i1-p1,i);
              fe.v = param[1] = cache.getParam(1,i2-p2,j);
              fe.w = param[2] = cache.getParam(2,i3-p3,k);

              // Fetch basis function derivatives at current integration point
              std::vector<const BasisFunctionVals*> bfs(this->getNoBasis());
              for (size_t b = 0; b < m_basis.size(); ++b) {
                bfs[b] = &myCache[b]->getVals(iel-1, ip);
                fe.basis(b+1) = bfs[b]->N;
              }

              // Compute Jacobian inverse of the coordinate mapping and
              // basis function derivatives w.r.t. Cartesian coordinates
              if (!fe.Jacobian(Jac,Xnod,geoBasis,&bfs))
                continue; // skip singular points

              // Compute Hessian of coordinate mapping and 2nd order derivatives
              if (use2ndDer && !fe.Hessian(Hess,Jac,Xnod,geoBasis,&bfs))
                ok = false;

              // Compute G-matrix
              if (integrand.getIntegrandType() & Integrand::G_MATRIX)
                utl::getGmat(Jac,dXidu,fe.G);

              // Cartesian coordinates of current integration point
              X.assign(Xnod * fe.basis(geoBasis));

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
    }

  for (std::unique_ptr<BasisFunctionCache>& cache : myCache)
    cache->finalizeAssembly();
  return ok;
}


bool ASMs3Dmx::integrate (Integrand& integrand, int lIndex,
			  GlobalIntegral& glInt,
			  const TimeDomain& time)
{
  if (!svol) return true; // silently ignore empty patches

  PROFILE2("ASMs3Dmx::integrate(B)");

  bool useElmVtx = integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS;

  std::map<char,ThreadGroups>::const_iterator tit;
  if ((tit = threadGroupsFace.find(lIndex%10)) == threadGroupsFace.end())
  {
    std::cerr <<" *** ASMs3Dmx::integrate: No thread groups for face "
              << lIndex%10 << std::endl;
    return false;
  }
  const ThreadGroups& threadGrp = tit->second;

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Find the parametric direction of the face normal {-3,-2,-1, 1, 2, 3}
  const int faceDir = (lIndex%10+1) / ((lIndex%2) ? -2 : 2);

  const int t1 = 1 + abs(faceDir)%3; // first tangent direction
  const int t2 = 1 + t1%3;           // second tangent direction

  // Compute parameter values of the Gauss points over the whole patch face
  std::array<Matrix,3> gpar;
  for (int d = 0; d < 3; d++)
    if (-1-d == faceDir)
    {
      gpar[d].resize(1,1);
      gpar[d].fill(svol->startparam(d));
    }
    else if (1+d == faceDir)
    {
      gpar[d].resize(1,1);
      gpar[d].fill(svol->endparam(d));
    }
    else
      this->getGaussPointParameters(gpar[d],d,nGauss,xg);

  // Extract the Neumann order flag (1 or higher) for the integrand
  integrand.setNeumannOrder(1 + lIndex/10);

  // Evaluate basis function derivatives at all integration points
  std::vector<std::vector<Go::BasisDerivs>> splinex(m_basis.size());
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < m_basis.size(); ++i)
    m_basis[i]->computeBasisGrid(gpar[0],gpar[1],gpar[2],splinex[i]);

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex%10);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;


  // === Assembly loop over all elements on the patch face =====================

  bool ok = true;
  for (size_t g = 0; g < threadGrp.size() && ok; g++)
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < threadGrp[g].size(); t++)
    {
      MxFiniteElement fe(elem_size);
      fe.xi = fe.eta = fe.zeta = faceDir < 0 ? -1.0 : 1.0;
      fe.u = gpar[0](1,1);
      fe.v = gpar[1](1,1);
      fe.w = gpar[2](1,1);

      Matrices dNxdu(m_basis.size());
      Matrix Xnod, Jac;
      double param[3] = { fe.u, fe.v, fe.w };
      Vec4   X(param,time.t);
      Vec3   normal;
      for (size_t l = 0; l < threadGrp[g][t].size() && ok; ++l)
      {
        int iel = threadGrp[g][t][l];
        fe.iel = MLGE[iel];
        if (fe.iel < 1) continue; // zero-volume element

        int i1 = p1 + iel % nel1;
        int i2 = p2 + (iel / nel1) % nel2;
        int i3 = p3 + iel / (nel1*nel2);

	// Get element face area in the parameter space
	double dA = 0.25*this->getParametricArea(++iel,abs(faceDir));
	if (dA < 0.0) // topology error (probably logic error)
	{
          ok = false;
          break;
        }

	// Set up control point coordinates for current element
	if (!this->getElementCoordinates(Xnod,iel))
	{
          ok = false;
          break;
        }

        if (useElmVtx)
          fe.h = this->getElementCorners(i1-1,i2-1,i3-1,fe.XC);

	// Initialize element quantities
        LocalIntegral* A = integrand.getLocalIntegral(elem_size,fe.iel,true);
	if (!integrand.initElementBou(MNPC[iel-1],elem_size,nb,*A))
        {
          A->destruct();
          ok = false;
          break;
        }

        // Define some loop control variables depending on which face we are on
        int nf1, j1, j2;
        switch (abs(faceDir))
        {
          case 1: nf1 = nel2; j2 = i3-p3; j1 = i2-p2; break;
          case 2: nf1 = nel1; j2 = i3-p3; j1 = i1-p1; break;
          case 3: nf1 = nel1; j2 = i2-p2; j1 = i1-p1; break;
          default: nf1 = j1 = j2 = 0;
        }


	// --- Integration loop over all Gauss points in each direction --------

        int k1, k2, k3;
        int ip = (j2*nGauss*nf1 + j1)*nGauss;
        int jp = (j2*nf1 + j1)*nGauss*nGauss;
        fe.iGP = firstp + jp; // Global integration point counter

        for (int j = 0; j < nGauss; j++, ip += nGauss*(nf1-1))
          for (int i = 0; i < nGauss; i++, ip++, fe.iGP++)
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
              fe.u = param[0] = gpar[0](k1+1,i1-p1+1);
            }
            if (gpar[1].size() > 1)
            {
              fe.eta = xg[k2];
              fe.v = param[1] = gpar[1](k2+1,i2-p2+1);
            }
            if (gpar[2].size() > 1)
            {
              fe.zeta = xg[k3];
              fe.w = param[2] = gpar[2](k3+1,i3-p3+1);
            }

            // Fetch basis function derivatives at current integration point
            for (size_t b = 0; b < m_basis.size(); ++b)
              SplineUtils::extractBasis(splinex[b][ip],fe.basis(b+1),dNxdu[b]);

            // Compute Jacobian inverse of the coordinate mapping and
            // basis function derivatives w.r.t. Cartesian coordinates
            fe.detJxW = utl::Jacobian(Jac,normal,fe.grad(geoBasis),Xnod,
                                      dNxdu[geoBasis-1],t1,t2);
            if (fe.detJxW == 0.0) continue; // skip singular points

            for (size_t b = 0; b < m_basis.size(); ++b)
              if (b != (size_t)geoBasis-1)
                fe.grad(b+1).multiply(dNxdu[b],Jac);

            if (faceDir < 0) normal *= -1.0;

            // Cartesian coordinates of current integration point
            X.assign(Xnod * fe.basis(geoBasis));

            // Evaluate the integrand and accumulate element contributions
            fe.detJxW *= dA*wg[i]*wg[j];
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
      }
    }

  return ok;
}


bool ASMs3Dmx::integrate (Integrand& integrand,
                          GlobalIntegral& glInt,
                          const TimeDomain& time,
                          const ASM::InterfaceChecker& iChk)
{
  if (!svol) return true; // silently ignore empty patches
  if (!(integrand.getIntegrandType() & Integrand::INTERFACE_TERMS))
    return true; // silently ignore if no interface terms
  else if (integrand.getIntegrandType() & Integrand::NORMAL_DERIVS)
  {
    std::cerr <<" *** Normal derivatives not implemented for mixed integrands."
              << std::endl;
    return false;
  }
  else if (MLGE.size() > nel && MLGE.size() != 2*nel)
  {
    std::cerr <<" *** Interface elements not implemented for mixed integrands."
              << std::endl;
    return false;
  }

  PROFILE2("ASMs3Dmx::integrate(J)");

  // Get Gaussian quadrature points and weights
  int nGP = integrand.getBouIntegrationPoints(nGauss);
  const double* xg = GaussQuadrature::getCoord(nGP);
  const double* wg = GaussQuadrature::getWeight(nGP);
  if (!xg || !wg) return false;

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);
  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);
  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;

  std::vector<size_t> elem_sizes2(elem_size);
  std::copy(elem_size.begin(), elem_size.end(), std::back_inserter(elem_sizes2));

  MxFiniteElement fe(elem_sizes2);
  Matrix        dNdu, Xnod, Jac;
  Vector        dN;
  Vec4          X(nullptr,time.t);
  Vec3          normal;
  double        u[2], v[2], w[2];

  // === Assembly loop over all elements in the patch ==========================

  int iel = 0;
  for (int i3 = p3; i3 <= n3; i3++)
    for (int i2 = p2; i2 <= n2; i2++)
      for (int i1 = p1; i1 <= n1; i1++, iel++)
      {
        fe.iel = abs(MLGE[iel]);
        if (fe.iel < 1) continue; // zero-area element

        short int status = iChk.hasContribution(iel,i1,i2,i3);
        if (!status) continue; // no interface contributions for this element

  #if SP_DEBUG > 3
        std::cout <<"\n\nIntegrating interface terms for element "<< fe.iel
                  << std::endl;
  #endif

        // Set up control point (nodal) coordinates for current element
        if (!this->getElementCoordinates(Xnod,1+iel)) return false;

        // Compute parameter values of the element edges
        this->getElementBorders(i1-1,i2-1,i3-1,u,v,w);

        if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
          fe.h = this->getElementCorners(i1-1,i2-1,i3-1,fe.XC);

        // Initialize element quantities
        LocalIntegral* A = integrand.getLocalIntegral(elem_size,fe.iel);
        bool ok = integrand.initElement(MNPC[iel],fe,elem_size,nb,*A);
        size_t origSize = A->vec.size();

        // Loop over the element edges with contributions
        int bit = 32;
        for (int iface = 6; iface > 0 && status > 0 && ok; iface--, bit /= 2)
          if (status & bit)
          {
            // Find the parametric direction of the edge normal {-2,-1, 1, 2}
            int faceDir;
            switch(iface) {
              case 1: faceDir = -1; break;
              case 2: faceDir = 1; break;
              case 3: faceDir = -2; break;
              case 4: faceDir = 2; break;
              case 5: faceDir = -3; break;
              case 6: faceDir = 3; break;
            }
            const int t1 = 1 + abs(faceDir)%3; // first tangent direction
            const int t2 = 1 + t1%3;           // second tangent direction

            int kel = iel;
            if (t1 == 2)
              kel += iface > 1  ? 1 : -1;
            else if (t1 == 3)
              kel += iface == 4 ? n1-p1+1 : -(n1-p1+1);
            else
              kel += iface == 6 ? (n2-p2+1)*(n1-p1+1) : -(n2-p2+1)*(n1-p1+1);

            // initialize neighbor element
            LocalIntegral* A_neigh = integrand.getLocalIntegral(elem_size,kel+1);
            ok &= integrand.initElement(MNPC[kel],fe,elem_size,nb,*A_neigh);
            if (!A_neigh->vec.empty()) {
              A->vec.resize(origSize+A_neigh->vec.size());
              std::copy(A_neigh->vec.begin(), A_neigh->vec.end(), A->vec.begin()+origSize);
            }
            A_neigh->destruct();

            // Get element area in the parameter space
            double dA = this->getParametricArea(1+iel,abs(faceDir));
            if (dA < 0.0) // topology error (probably logic error)
              ok = false;

            // Define some loop control variables depending on which face we are on
            int nf1, j1, j2;
            switch (abs(faceDir))
            {
              case 1: nf1 = nel2; j2 = i3-p3; j1 = i2-p2; break;
              case 2: nf1 = nel1; j2 = i3-p3; j1 = i1-p1; break;
              case 3: nf1 = nel1; j2 = i2-p2; j1 = i1-p1; break;
              default: nf1 = j1 = j2 = 0;
            }

            int ip = (j2*nGP*nf1 + j1)*nGP;

            // --- Integration loop over all Gauss points along the face ---------

            for (int j = 0; j < nGP && ok; j++, ip += nGP*(nf1-1))
              for (int i = 0; i < nGP && ok; i++, ip++, fe.iGP++)
              {
                // Local element coordinates and parameter values
                // of current integration point
                if (abs(faceDir) == 1)
                {
                  fe.xi = faceDir > 0 ? 1.0 : 0.0;
                  fe.eta = xg[i];
                  fe.zeta = xg[j];
                  fe.u = faceDir > 0 ? u[1] : u[0];
                  fe.v = 0.5*((v[1]-v[0])*xg[i] + v[1]+v[0]);
                  fe.w = 0.5*((w[1]-w[0])*xg[j] + w[1]+w[0]);
                  fe.p = p1 - 1;
                }
                else if (abs(faceDir) == 2)
                {
                  fe.xi = xg[i];
                  fe.eta = faceDir > 0 ? 1.0 : 0.0;
                  fe.zeta = xg[j];
                  fe.u = 0.5*((u[1]-u[0])*xg[i] + u[1]+u[0]);
                  fe.v = faceDir > 0 ? v[1] : v[0];
                  fe.w = 0.5*((w[1]-w[0])*xg[j] + w[1]+w[0]);
                  fe.p = p2 - 1;
                }
                else
                {
                  fe.xi = xg[i];
                  fe.eta = xg[j];
                  fe.zeta = faceDir > 0 ? 1.0 : 0.0;
                  fe.u = 0.5*((u[1]-u[0])*xg[i] + u[1]+u[0]);
                  fe.v = 0.5*((v[1]-v[0])*xg[j] + v[1]+v[0]);
                  fe.w = faceDir > 0 ? w[1] : w[0];
                  fe.p = p3 - 1;
                }

                // Fetch basis function derivatives at current integration point
                Matrices dNxdu(m_basis.size()*2);
                for (size_t b = 0; b < m_basis.size(); ++b) {
                  Go::BasisDerivs spline;
                  this->getBasis(b+1)->computeBasis(fe.u, fe.v, fe.w, spline, faceDir < 0);
                  SplineUtils::extractBasis(spline, fe.basis(b+1), dNxdu[b]);
                  this->getBasis(b+1)->computeBasis(fe.u, fe.v, fe.w, spline, faceDir > 0);
                  SplineUtils::extractBasis(spline, fe.basis(b+1+m_basis.size()),
                                            dNxdu[b+m_basis.size()]);
                }

              // Compute basis function derivatives and the edge normal
              fe.detJxW = utl::Jacobian(Jac,normal,fe.grad(geoBasis+m_basis.size()),
                                        Xnod,dNxdu[geoBasis-1+m_basis.size()],t1,t2);
              fe.detJxW = utl::Jacobian(Jac,normal,fe.grad(geoBasis),Xnod,
                                        dNxdu[geoBasis-1],t1,t2);
              if (fe.detJxW == 0.0) continue; // skip singular points
              for (size_t b = 0; b < m_basis.size(); ++b)
                if (b != (size_t)geoBasis-1) {
                  fe.grad(b+1).multiply(dNxdu[b],Jac);
                  fe.grad(b+1+m_basis.size()).multiply(dNxdu[b+m_basis.size()],Jac);
                }

              if (faceDir < 0) normal *= -1.0;

              // Cartesian coordinates of current integration point
              X.assign(Xnod * fe.basis(geoBasis));

              // Evaluate the integrand and accumulate element contributions
              fe.detJxW *= 0.25*dA*wg[i]*wg[j];
              ok = integrand.evalIntMx(*A,fe,time,X,normal);
            }
          }

        // Finalize the element quantities
        if (ok && !integrand.finalizeElement(*A,fe,time))
          ok = false;

        // Assembly of global system integral
        if (ok && !glInt.assemble(A->ref(),fe.iel))
          ok = false;

        A->destruct();

        if (!ok) return false;
      }

  return true;
}


int ASMs3Dmx::evalPoint (const double* xi, double* param, Vec3& X) const
{
  if (!svol) return -3;

  int i;
  for (i = 0; i < 3; i++)
    param[i] = (1.0-xi[i])*svol->startparam(i) + xi[i]*svol->endparam(i);

  Go::Point X0;
  svol->point(X0,param[0],param[1],param[2]);
  for (i = 0; i < 3 && i < svol->dimension(); i++)
    X[i] = X0[i];

  // Check if this point matches any of the control points (nodes)
  return this->searchCtrlPt(m_basis.front()->coefs_begin(),
                            m_basis.front()->coefs_end(), X,
                            m_basis.front()->dimension());
}


bool ASMs3Dmx::evalSolution (Matrix& sField, const Vector& locSol,
                             const RealArray* gpar,
                             bool regular, int, int nf) const
{
  // Evaluate the basis functions at all points
  std::vector<std::vector<Go::BasisPts>> splinex(m_basis.size());
  if (regular)
  {
    for (size_t b = 0; b < m_basis.size(); ++b)
      m_basis[b]->computeBasisGrid(gpar[0],gpar[1],gpar[2],splinex[b]);
  }
  else if (gpar[0].size() == gpar[1].size() && gpar[0].size() == gpar[2].size())
  {
    for (size_t b = 0; b < m_basis.size(); ++b) {
      splinex[b].resize(gpar[0].size());
      for (size_t i = 0; i < splinex[b].size(); i++)
        m_basis[b]->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],splinex[b][i]);
    }
  }
  else
    return false;

  std::vector<size_t> nc(nfx.size(), 0);
  if (nf)
    nc.front() = nf;
  else
    std::copy(nfx.begin(), nfx.end(), nc.begin());

  if (std::inner_product(nb.begin(), nb.end(), nc.begin(), 0u) != locSol.size())
    return false;

  Matrix Xtmp;
  Vector Ytmp, Ztmp;

  // Evaluate the primary solution field at each point
  size_t nPoints = splinex.front().size();
  sField.resize(std::accumulate(nc.begin(), nc.end(), 0),nPoints);
  for (size_t i = 0; i < nPoints; i++)
  {
    size_t comp = 0;
    for (size_t b = 0; b < m_basis.size(); ++b) {
      if (nc[b] == 0)
        continue;
      IntVec ip;
      scatterInd(m_basis[b]->numCoefs(0),m_basis[b]->numCoefs(1),m_basis[b]->numCoefs(2),
                 m_basis[b]->order(0), m_basis[b]->order(1), m_basis[b]->order(2),
                 splinex[b][i].left_idx,ip);

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


bool ASMs3Dmx::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			     const RealArray* gpar, bool regular) const
{
  sField.resize(0,0);

  // Evaluate the basis functions and their derivatives at all points
  std::vector<std::vector<Go::BasisDerivs>> splinex(m_basis.size());
  if (regular)
  {
    for (size_t b = 0; b < m_basis.size(); ++b)
      m_basis[b]->computeBasisGrid(gpar[0],gpar[1],gpar[2],splinex[b]);
  }
  else if (gpar[0].size() == gpar[1].size() && gpar[0].size() == gpar[2].size())
  {
    for (size_t b = 0; b < m_basis.size(); ++b) {
      splinex[b].resize(gpar[0].size());
      for (size_t i = 0; i < splinex[b].size(); i++)
        m_basis[b]->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],splinex[b][i]);
    }
  }

  // Fetch nodal (control point) coordinates
  Matrix Xnod, Xtmp;
  this->getNodalCoordinates(Xnod);

  MxFiniteElement fe(elem_size,firstIp);
  Vector          solPt;
  Matrices        dNxdu(m_basis.size());
  Matrix          Jac;

  // Evaluate the secondary solution field at each point
  size_t nPoints = splinex.front().size();
  for (size_t i = 0; i < nPoints; i++, fe.iGP++)
  {
    // Fetch indices of the non-zero basis functions at this point
    IntMat ip(m_basis.size());
    IntVec ipa;
    size_t ofs = 0;
    for (size_t b = 0; b < m_basis.size(); ++b) {
      scatterInd(m_basis[b]->numCoefs(0),m_basis[b]->numCoefs(1),m_basis[b]->numCoefs(2),
                 m_basis[b]->order(0),m_basis[b]->order(1),m_basis[b]->order(2),
                 splinex[b][i].left_idx,ip[b]);

      // Fetch associated control point coordinates
      if (b == (size_t)geoBasis-1)
        utl::gather(ip[geoBasis-1], 3, Xnod, Xtmp);

      for (int& c : ip[b]) c += ofs;
      ipa.insert(ipa.end(), ip[b].begin(), ip[b].end());
      ofs += nb[b];
    }

    fe.u = splinex[0][i].param[0];
    fe.v = splinex[0][i].param[1];
    fe.w = splinex[0][i].param[2];

    // Fetch basis function derivatives at current integration point
    for (size_t b = 0; b < m_basis.size(); ++b)
      SplineUtils::extractBasis(splinex[b][i],fe.basis(b+1),dNxdu[b]);

    // Compute Jacobian inverse of the coordinate mapping and
    // basis function derivatives w.r.t. Cartesian coordinate
    if (!fe.Jacobian(Jac,Xtmp,geoBasis,nullptr,&dNxdu))
      continue; // skip singular points

    // Cartesian coordinates of current integration point
    utl::Point X4(Xtmp * fe.basis(geoBasis),{fe.u,fe.v,fe.w});

    // Now evaluate the solution field
    if (!integrand.evalSol(solPt,fe,X4,ipa,elem_size,nb))
      return false;
    else if (sField.empty())
      sField.resize(solPt.size(),nPoints,true);

    sField.fillColumn(1+i,solPt);
  }

  return true;
}


void ASMs3Dmx::generateThreadGroups (const Integrand& integrand, bool silence,
                                     bool ignoreGlobalLM)
{
  int p[3] = { 0, 0, 0 };
  for (const auto& it : m_basis)
    for (size_t d = 0; d < 3; d++)
      if (it->order(d) > p[d])
        p[d] = it->order(d);

  this->ASMs3D::generateThreadGroups(p[0]-1, p[1]-1, p[2]-1,
                                     silence, ignoreGlobalLM);
}


void ASMs3Dmx::getBoundaryNodes (int lIndex, IntVec& nodes, int basis,
                                 int thick, int, bool local) const
{
  if (basis > 0)
    this->ASMs3D::getBoundaryNodes(lIndex, nodes, basis, thick, 0, local);
  else
    for (size_t b = 1; b <= m_basis.size(); b++)
      this->ASMs3D::getBoundaryNodes(lIndex, nodes, b, thick, 0, local);
}


void ASMs3Dmx::swapProjectionBasis ()
{
  if (projB2) {
    ASMmxBase::geoBasis = ASMmxBase::geoBasis == 1 ? 2 : 1;
    std::swap(projB, projB2);
    svol = this->getBasis(ASMmxBase::geoBasis);
  }
}


int ASMs3Dmx::getFirstItgElmNode () const
{
  return std::accumulate(elem_size.begin(), elem_size.begin() + geoBasis-1, 0);
}


int ASMs3Dmx::getLastItgElmNode () const
{
  return std::accumulate(elem_size.begin(), elem_size.begin() + geoBasis, -1);
}


bool ASMs3Dmx::separateProjectionBasis () const
{
  return std::none_of(m_basis.begin(), m_basis.end(),
                      [this](const std::shared_ptr<Go::SplineVolume>& entry)
                      { return entry.get() == projB; });
}
