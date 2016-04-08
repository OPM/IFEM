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

#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/VolumeInterpolator.h"

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
#include "Profiler.h"
#include "Vec3Oper.h"
#include "Vec3.h"
#include <array>
#include <numeric>
#ifdef USE_OPENMP
#include <omp.h>
#endif


ASMs3Dmx::ASMs3Dmx (const std::vector<unsigned char>& n_f)
  : ASMs3D(std::accumulate(n_f.begin(), n_f.end(), 0)), ASMmxBase(n_f)
{
}


ASMs3Dmx::ASMs3Dmx (const ASMs3Dmx& patch,
                    const std::vector<unsigned char>& n_f)
  : ASMs3D(patch), ASMmxBase(n_f)
{
  m_basis = patch.m_basis;
}


Go::SplineVolume* ASMs3Dmx::getBasis (int basis) const
{
  if (basis < 1 || basis > (int)m_basis.size())
    basis = 1;

  return m_basis[basis-1].get();
}


Go::SplineSurface* ASMs3Dmx::getBoundary (int dir, int basis)
{
  if (dir < -3 || dir == 0 || dir > 3)
    return nullptr;

  // The boundary surfaces are stored internally in the SplineVolume object
  int iface = dir > 0 ? 2*dir-1 : -2*dir-2;
  return m_basis[basis-1]->getBoundarySurface(iface).get();
}


bool ASMs3Dmx::write (std::ostream& os, int basis) const
{
  if (m_basis[0] && basis == 1)
    os <<"700 1 0 0\n" << *m_basis[0];
  else if (m_basis[1] && basis == 2)
    os <<"700 1 0 0\n" << *m_basis[1];
  else if (svol)
    os <<"700 1 0 0\n" << *svol;
  else
    return false;

  return os.good();
}


void ASMs3Dmx::clear (bool retainGeometry)
{
  // Erase the spline data
  if (!retainGeometry)
    svol = 0;

  if (!shareFE)
    m_basis[0].reset();
  if (!shareFE)
    m_basis[1].reset();

  // Erase the FE data
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
  if (basis < 1 || basis > (int)nfx.size())
    return std::accumulate(nfx.begin(), nfx.end(), 0);

  return nfx[basis-1];
}


unsigned char ASMs3Dmx::getNodalDOFs (size_t inod) const
{
  if (this->isLMn(inod)) return nLag;
  size_t nbc=0;
  for (size_t i=0;i<nb.size();++i)
    if (inod <= (nbc+=nb[i]))
      return nfx[i];

  return nfx[0];
}


char ASMs3Dmx::getNodeType (size_t inod) const
{
  if (this->isLMn(inod)) return 'L';
  size_t nbc=nb.front();
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


void ASMs3Dmx::extractNodeVec (const Vector& globRes, Vector& nodeVec,
			       unsigned char, int basis) const
{
  this->extractNodeVecMx(globRes,nodeVec,basis);
}


bool ASMs3Dmx::injectNodeVec (const Vector& nodeRes, Vector& globRes,
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

  if (m_basis.empty())
    m_basis = ASMmxBase::establishBases(svol, ASMmxBase::Type);

  nb.clear();
  for (auto it : m_basis)
    nb.push_back(it->numCoefs(0)*it->numCoefs(1)*it->numCoefs(2));

  if (!nodeInd.empty() && !shareFE)
  {
    if (nodeInd.size() == std::accumulate(nb.begin(), nb.end(), 0u)) return true;
    std::cerr <<" *** ASMs3Dmx::generateFEMTopology: Inconsistency between the"
	      <<" number of FE nodes "<< nodeInd.size()
	      <<"\n     and the number of spline coefficients "
              << std::accumulate(nb.begin(), nb.end(), 0)
	      <<" in the patch."<< std::endl;
    return false;
  }

  if (shareFE == 'F') return true;

#ifdef SP_DEBUG
  for (auto it : m_basis) {
    std::cout <<"numCoefs: "<< it->numCoefs(0) <<" "
                            << it->numCoefs(1) <<" "<< it->numCoefs(2);
    std::cout <<"\norder: "<< it->order(0) <<" "
                           << it->order(1) <<" "<< it->order(2);
    for (int d = 0; d < 3; d++)
    {
      std::cout <<"\nd"<< char('u'+d) <<':';
      for (int i = 0; i < it->numCoefs(d); i++)
        std::cout <<' '<< it->knotSpan(d,i);
    }
    std::cout << std::endl;
  }
#endif

  nel = (m_basis[geoBasis-1]->numCoefs(0)-m_basis[geoBasis-1]->order(0)+1)*
        (m_basis[geoBasis-1]->numCoefs(1)-m_basis[geoBasis-1]->order(1)+1)*
        (m_basis[geoBasis-1]->numCoefs(2)-m_basis[geoBasis-1]->order(2)+1);

  nnod = std::accumulate(nb.begin(), nb.end(), 0);

  myMLGE.resize(nel,0);
  myMLGN.resize(nnod);
  myMNPC.resize(nel);
  myNodeInd.resize(nnod);

  int i1, i2, i3, j1, j2, j3;
  size_t iel, inod = 0;
  for (auto& it : m_basis) {
    for (i3 = 0; i3 < it->numCoefs(2); i3++)
      for (i2 = 0; i2 <  it->numCoefs(1); i2++)
        for (i1 = 0; i1 < it->numCoefs(0); i1++)
        {
          myNodeInd[inod].I = i1;
          myNodeInd[inod].J = i2;
          myNodeInd[inod].K = i3;
          myMLGN[inod++]    = ++gNod;
        }
  }

  int lnod2 = 0;
  int lnod3 = 0;

  auto&& addBasis = [&](std::shared_ptr<Go::SplineVolume>& b, bool geo) {
    for (i3 = 1; i3 <= b->numCoefs(2); i3++)
      for (i2 = 1; i2 <= b->numCoefs(1); i2++)
	for (i1 = 1; i1 <= b->numCoefs(0); i1++, inod++)
	  if (i1 >= b->order(0) && i2 >= b->order(1) && i3 >= b->order(2))
	  {
	    if (b->knotSpan(0,i1-1) > 0.0)
	      if (b->knotSpan(1,i2-1) > 0.0)
		if (b->knotSpan(2,i3-1) > 0.0)
		{
                  if (geo)
                    myMLGE[iel] = ++gEl; // global element number over all patches
                  else
                    while (iel < myMNPC.size() && myMNPC[iel].empty()) iel++;

		  int lnod = lnod2;
		  myMNPC[iel].resize(lnod3,0);
		  for (j3 = b->order(2)-1; j3 >= 0; j3--)
		    for (j2 = b->order(1)-1; j2 >= 0; j2--)
		      for (j1 = b->order(0)-1; j1 >= 0; j1--)
			myMNPC[iel][lnod++] = inod - b->numCoefs(0)*b->numCoefs(1)*j3 - b->numCoefs(0)*j2 - j1;

                  if (!geo)
                    ++iel;
		}

            if (geo)
              ++iel;
	  }
  };

  iel = 0, inod = std::accumulate(nb.begin(),nb.begin()+geoBasis-1,0u);

  for (i2 = 0; i2 < geoBasis-1; ++i2)
    lnod2 += m_basis[i2]->order(0)*m_basis[i2]->order(1)*m_basis[i2]->order(2);
  for (i2 = 0; i2 < (int)m_basis.size(); ++i2)
    lnod3 += m_basis[i2]->order(0)*m_basis[i2]->order(1)*m_basis[i2]->order(2);

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
    lnod2 += m_basis[b]->order(0)*m_basis[b]->order(1)*m_basis[b]->order(2);
  }

#ifdef SP_DEBUG
  std::cout <<"NEL = "<< nel <<" NNOD = "<< nnod << std::endl;
#endif
  return true;
}


bool ASMs3Dmx::connectPatch (int face, ASMs3D& neighbor, int nface, int norient)
{
  ASMs3Dmx* neighMx = dynamic_cast<ASMs3Dmx*>(&neighbor);
  if (!neighMx) return false;

  if (swapW && face > 4) // Account for swapped parameter direction
    face = 11-face;

  if (neighMx->swapW && face > 4) // Account for swapped parameter direction
    nface = 11-nface;

  size_t nbi=0;
  for (size_t i = 1;i <= m_basis.size(); ++i) {
    if (!this->connectBasis(face,neighbor,nface,norient,i,nbi,nbi))
      return false;
    nbi += nb[i-1];
  }

  this->addNeighbor(neighMx);
  return true;
}


void ASMs3Dmx::closeFaces (int dir, int, int)
{
  size_t nbi = 0;
  for (size_t i = 1;i <= m_basis.size(); ++i) {
    this->ASMs3D::closeFaces(dir,i,nbi+1);
    nbi += nb[i-1];
  }
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
  size_t nbb = nb[0];
  while (nbb < inod)
    nbb += nb[++b];

  cit = m_basis[b]->coefs_begin()
      + ((K*m_basis[b]->numCoefs(1)+J)*m_basis[b]->numCoefs(0)+I) * m_basis[b]->dimension();

  return Vec3(*cit,*(cit+1),*(cit+2));
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
  if (m_basis.empty()) return false;

  PROFILE2("ASMs3Dmx::integrate(I)");

  bool useElmVtx = integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS;

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  std::array<Matrix,3> gpar;
  for (int d = 0; d < 3; d++)
    this->getGaussPointParameters(gpar[d],d,nGauss,xg);

  // Evaluate basis function derivatives at all integration points
  std::vector<std::vector<Go::BasisDerivs>> splinex(m_basis.size());
#pragma omp parallel for schedule(static)
  for (size_t i=0;i<m_basis.size();++i)
    m_basis[i]->computeBasisGrid(gpar[0],gpar[1],gpar[2],splinex[i]);

  std::vector<size_t> elem_sizes;
  for (auto& it : m_basis)
    elem_sizes.push_back(it->order(0)*it->order(1)*it->order(2));

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);
  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;
  const int nel3 = n3 - p3 + 1;

  // === Assembly loop over all elements in the patch ==========================

  bool ok=true;
  for (size_t g=0;g<threadGroupsVol.size() && ok;++g) {
#pragma omp parallel for schedule(static)
    for (size_t t=0;t<threadGroupsVol[g].size();++t) {
      MxFiniteElement fe(elem_sizes);
      std::vector<Matrix> dNxdu(m_basis.size());
      Matrix Xnod, Jac;
      Vec4   X;
      for (size_t l = 0; l < threadGroupsVol[g][t].size() && ok; ++l)
      {
        int iel = threadGroupsVol[g][t][l];
        fe.iel = MLGE[iel];
        if (fe.iel < 1) continue; // zero-volume element

        int i1 = p1 + iel % nel1;
        int i2 = p2 + (iel / nel1) % nel3;
        int i3 = p3 + iel / (nel1*nel2);

        // Get element volume in the parameter space
        double dV = this->getParametricVolume(++iel);
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
          this->getElementCorners(i1-1,i2-1,i3-1,fe.XC);

        // Initialize element quantities
        LocalIntegral* A = integrand.getLocalIntegral(elem_sizes,fe.iel,false);
        if (!integrand.initElement(MNPC[iel-1],elem_sizes,nb,*A))
        {
          A->destruct();
          ok = false;
          break;
        }


        // --- Integration loop over all Gauss points in each direction --------

        int ip = (((i3-p3)*nGauss*nel2 + i2-p2)*nGauss*nel1 + i1-p1)*nGauss;
        int jp = (((i3-p3)*nel2 + i2-p2*nel1 + i1-p1))*nGauss*nGauss*nGauss;
        fe.iGP = firstIp + jp; // Global integration point counter

        for (int k = 0; k < nGauss; k++, ip += nGauss*(nel2-1)*nGauss*nel1)
          for (int j = 0; j < nGauss; j++, ip += nGauss*(nel1-1))
            for (int i = 0; i < nGauss; i++, ip++, fe.iGP++)
            {
              // Local element coordinates of current integration point
              fe.xi   = xg[i];
              fe.eta  = xg[j];
              fe.zeta = xg[k];

              // Parameter values of current integration point
              fe.u = gpar[0](i+1,i1-p1+1);
              fe.v = gpar[1](j+1,i2-p2+1);
              fe.w = gpar[2](k+1,i3-p3+1);

              // Fetch basis function derivatives at current integration point
              for (size_t b = 0; b < m_basis.size(); ++b)
                SplineUtils::extractBasis(splinex[b][ip],fe.basis(b+1),dNxdu[b]);

              // Compute Jacobian inverse of the coordinate mapping and
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
              fe.detJxW *= 0.125*dV*wg[i]*wg[j]*wg[k];
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


bool ASMs3Dmx::integrate (Integrand& integrand, int lIndex,
			  GlobalIntegral& glInt,
			  const TimeDomain& time)
{
  if (!svol) return true; // silently ignore empty patches
  if (m_basis.empty()) return false;

  PROFILE2("ASMs3Dmx::integrate(B)");

  bool useElmVtx = integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS;

  std::map<char,ThreadGroups>::const_iterator tit;
  if ((tit = threadGroupsFace.find(lIndex)) == threadGroupsFace.end())
  {
    std::cerr <<" *** ASMs3D::integrate: No thread groups for face "<< lIndex
	      << std::endl;
    return false;
  }
  const ThreadGroups& threadGrp = tit->second;

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Find the parametric direction of the face normal {-3,-2,-1, 1, 2, 3}
  const int faceDir = (lIndex+1)/(lIndex%2 ? -2 : 2);

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

  std::vector<size_t> elem_sizes;
  for (auto& it : m_basis)
    elem_sizes.push_back(it->order(0)*it->order(1)*it->order(2));

  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;


  // === Assembly loop over all elements on the patch face =====================

  bool ok = true;
  for (size_t g = 0; g < threadGrp.size() && ok; ++g) {
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < threadGrp[g].size(); ++t) {
      MxFiniteElement fe(elem_sizes);
      fe.xi = fe.eta = fe.zeta = faceDir < 0 ? -1.0 : 1.0;
      fe.u = gpar[0](1,1);
      fe.v = gpar[1](1,1);
      fe.w = gpar[2](1,1);

      std::vector<Matrix> dNxdu(m_basis.size());
      Matrix Xnod, Jac;
      Vec4   X;
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
	double dA = this->getParametricArea(++iel,abs(faceDir));
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
          this->getElementCorners(i1-1,i2-1,i3-1,fe.XC);

	// Initialize element quantities
        LocalIntegral* A = integrand.getLocalIntegral(elem_sizes,fe.iel,true);
	if (!integrand.initElementBou(MNPC[iel-1],elem_sizes,nb,*A))
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
              case 1: k2 = i+1; k3 = j+1; k1 = 0; break;
              case 2: k1 = i+1; k3 = j+1; k2 = 0; break;
              case 3: k1 = i+1; k2 = j+1; k3 = 0; break;
              default: k1 = k2 = k3 = 0;
            }
            if (gpar[0].size() > 1)
            {
              fe.xi = xg[k1];
              fe.u = gpar[0](k1,i1-p1+1);
            }
            if (gpar[1].size() > 1)
            {
              fe.eta = xg[k2];
              fe.v = gpar[1](k2,i2-p2+1);
            }
            if (gpar[2].size() > 1)
            {
              fe.zeta = xg[k3];
              fe.w = gpar[2](k3,i3-p3+1);
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

            if (faceDir < 0) normal *= -1.0;

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
      }
    }
  }

  return ok;
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
  return this->searchCtrlPt(m_basis[0]->coefs_begin(),m_basis[0]->coefs_end(),
                            X,m_basis[0]->dimension());
}


bool ASMs3Dmx::evalSolution (Matrix& sField, const Vector& locSol,
                             const RealArray* gpar, bool regular, int) const
{
  if (m_basis.empty()) return false;

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

  std::vector<size_t> nc(nfx.size());
  std::copy(nfx.begin(), nfx.end(), nc.begin());

  // assume geobasis only
  if (locSol.size() < std::inner_product(nb.begin(), nb.end(), nfx.begin(), 0u)) {
    std::fill(nc.begin(), nc.end(), 0);
    nc[geoBasis-1] = nfx[geoBasis-1];
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

  if (m_basis.empty()) return false;

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

  std::vector<size_t> elem_sizes;
  for (auto& it : m_basis)
    elem_sizes.push_back(it->order(0)*it->order(1)*it->order(2));

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
      scatterInd(m_basis[b]->numCoefs(0),m_basis[b]->numCoefs(1),m_basis[b]->numCoefs(2),
                 m_basis[b]->order(0),m_basis[b]->order(1),m_basis[b]->order(2),
                 splinex[b][i].left_idx,ip[b]);

      // Fetch associated control point coordinates
      if (b == (size_t)geoBasis-1)
        utl::gather(ip[geoBasis-1], 3, Xnod, Xtmp);

      for (auto& it : ip[b])
        it += ofs;
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


void ASMs3Dmx::generateThreadGroups (const Integrand& integrand, bool silence)
{
  int p[3] = { 0, 0, 0 };
  for (const auto& it : m_basis)
    for (size_t d = 0; d < 3; d++)
      if (it->order(d) > p[d])
        p[d] = it->order(d);

  ASMs3D::generateThreadGroups(p[0]-1,p[1]-1,p[2]-1,silence);
}


void ASMs3Dmx::generateThreadGroups (char lIndex, bool silence)
{
#ifdef USE_OPENMP
  omp_set_num_threads(1);
#endif
  ASMs3D::generateThreadGroups(lIndex,silence);
}


#define DERR -999.99

double ASMs3Dmx::getParametricVolume (int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs3D::getParametricVolume: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return DERR;
  }
#endif
  if (MNPC[iel-1].empty())
    return 0.0;

  std::vector<size_t> elem_sizes;
  for (auto& it : m_basis)
    elem_sizes.push_back(it->order(0)*it->order(1)*it->order(2));

  int inod1 = MNPC[iel-1][std::accumulate(elem_sizes.begin(),
                                          elem_sizes.begin()+geoBasis, -1)];
#ifdef INDEX_CHECK
  if (inod1 < 0 || (size_t)inod1 >= nnod)
  {
    std::cerr <<" *** ASMs3Dmx::getParametricVolume: Node index "<< inod1
	      <<" out of range [0,"<< nnod <<">."<< std::endl;
    return DERR;
  }
#endif

  double du = svol->knotSpan(0,nodeInd[inod1].I);
  double dv = svol->knotSpan(1,nodeInd[inod1].J);
  double dw = svol->knotSpan(2,nodeInd[inod1].K);
  return du*dv*dw;
}


double ASMs3Dmx::getParametricArea (int iel, int dir) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs3D::getParametricArea: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return DERR;
  }
#endif
  if (MNPC[iel-1].empty())
    return 0.0;

  std::vector<size_t> elem_sizes;
  for (auto& it : m_basis)
    elem_sizes.push_back(it->order(0)*it->order(1)*it->order(2));

  int inod1 = MNPC[iel-1][std::accumulate(elem_sizes.begin(),
                                          elem_sizes.begin()+geoBasis, -1)];
#ifdef INDEX_CHECK
  if (inod1 < 0 || (size_t)inod1 >= nnod)
  {
    std::cerr <<" *** ASMs3D::getParametricArea: Node index "<< inod1
	      <<" out of range [0,"<< nnod <<">."<< std::endl;
    return DERR;
  }
#endif

  const int ni = nodeInd[inod1].I;
  const int nj = nodeInd[inod1].J;
  const int nk = nodeInd[inod1].K;
  switch (dir)
    {
    case 1: return svol->knotSpan(1,nj)*svol->knotSpan(2,nk);
    case 2: return svol->knotSpan(0,ni)*svol->knotSpan(2,nk);
    case 3: return svol->knotSpan(0,ni)*svol->knotSpan(1,nj);
    }

  std::cerr <<" *** ASMs3D::getParametricArea: Invalid face direction "
	    << dir << std::endl;
  return DERR;
}


void ASMs3Dmx::getBoundaryNodes (int lIndex, IntVec& nodes, int basis) const
{
  for (int b = (basis==0?1:basis); b <= (basis==0?(int)getNoBasis():basis); ++b)
    ASMs3D::getBoundaryNodes(lIndex, nodes, b);
}
