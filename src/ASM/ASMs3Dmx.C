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
#include "SplineFields3D.h"
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


ASMs3Dmx::ASMs3Dmx (const CharVec& n_f)
  : ASMs3D(std::accumulate(n_f.begin(), n_f.end(), 0)), ASMmxBase(n_f)
{
}


ASMs3Dmx::ASMs3Dmx (const ASMs3Dmx& patch, const CharVec& n_f)
  : ASMs3D(patch), ASMmxBase(n_f), m_basis(patch.m_basis)
{
}


Go::SplineVolume* ASMs3Dmx::getBasis (int basis) const
{
  if (basis < 1 || basis > (int)m_basis.size())
    return svol;

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
  if (basis == -1)
    os <<"700 1 0 0\n" << *projBasis;
  else if (basis < 1 || basis > (int)m_basis.size())
    os <<"700 1 0 0\n" << *svol;
  else if (m_basis[basis-1])
    os <<"700 1 0 0\n" << *m_basis[basis-1];
  else
    return false;

  return os.good();
}


void ASMs3Dmx::clear (bool retainGeometry)
{
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
  if (basis < 1 || basis > (int)nfx.size())
    return std::accumulate(nfx.begin(), nfx.end(), 0);

  return nfx[basis-1];
}


unsigned char ASMs3Dmx::getNodalDOFs (size_t inod) const
{
  if (this->isLMn(inod))
    return nLag;

  size_t nbc=0;
  for (size_t i=0;i<nb.size();++i)
    if (inod <= (nbc+=nb[i]))
      return nfx[i];

  return nfx[0];
}


char ASMs3Dmx::getNodeType (size_t inod) const
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

  if (m_basis.empty()) {
    m_basis = ASMmxBase::establishBases(svol, ASMmxBase::Type);

    // we need to project on something that is not one of our bases
    if (ASMmxBase::Type == ASMmxBase::REDUCED_CONT_RAISE_BASIS1 ||
        ASMmxBase::Type == ASMmxBase::DIV_COMPATIBLE)
      projBasis = ASMmxBase::establishBases(svol,
                                            ASMmxBase::FULL_CONT_RAISE_BASIS1).front();
    else
      projBasis = m_basis.front();
  }

  delete svol;
  geomB = svol = m_basis[geoBasis-1]->clone();

  nb.clear();
  elem_size.clear();
  nb.reserve(m_basis.size());
  elem_size.reserve(m_basis.size());
  for (auto& it : m_basis) {
    nb.push_back(it->numCoefs(0)*it->numCoefs(1)*it->numCoefs(2));
    elem_size.push_back(it->order(0)*it->order(1)*it->order(2));
  }

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

  iel = 0, inod = std::accumulate(nb.begin(),nb.begin()+geoBasis-1,0u);

  for (i2 = 0; i2 < geoBasis-1; ++i2)
    lnod2 += m_basis[i2]->order(0)*m_basis[i2]->order(1)*m_basis[i2]->order(2);
  for (i2 = 0; i2 < (int)m_basis.size(); ++i2)
    lnod3 += m_basis[i2]->order(0)*m_basis[i2]->order(1)*m_basis[i2]->order(2);

  // Create nodal connectivities for bases
  Go::SplineVolume* b = m_basis[geoBasis-1].get();
  auto knotw = b->basis(2).begin();
  for (i3 = 1; i3 <= b->numCoefs(2); i3++, ++knotw) {
    auto knotv = b->basis(1).begin();
    for (i2 = 1; i2 <= b->numCoefs(1); i2++, ++knotv) {
      auto knotu = b->basis(1).begin();
      for (i1 = 1; i1 <= b->numCoefs(0); i1++, inod++, ++knotu)
        if (i1 >= b->order(0) && i2 >= b->order(1) && i3 >= b->order(2))
        {
          if (b->knotSpan(0,i1-1) > 0.0)
            if (b->knotSpan(1,i2-1) > 0.0)
              if (b->knotSpan(2,i3-1) > 0.0)
              {
                myMLGE[iel] = ++gEl; // global element number over all patches

                int lnod = lnod2;
                myMNPC[iel].resize(lnod3,0);
                for (j3 = b->order(2)-1; j3 >= 0; j3--)
                  for (j2 = b->order(1)-1; j2 >= 0; j2--)
                    for (j1 = b->order(0)-1; j1 >= 0; j1--)
                      myMNPC[iel][lnod++] = inod - b->numCoefs(0)*b->numCoefs(1)*j3 - b->numCoefs(0)*j2 - j1;
                // find knotspan spanning element for other bases
                lnod = 0;
                size_t lnod4 = 0;
                for (size_t bas = 0; bas <  m_basis.size(); ++bas) {
                  if (bas != (size_t)geoBasis-1) {
                    double ku = *knotu;
                    double kv = *knotv;
                    double kw = *knotw;
                    int bknotu = m_basis[bas]->basis(0).knotIntervalFuzzy(ku);
                    int bknotv = m_basis[bas]->basis(1).knotIntervalFuzzy(kv);
                    int bknotw = m_basis[bas]->basis(2).knotIntervalFuzzy(kw);
                    size_t iinod = lnod4+(bknotw*m_basis[bas]->numCoefs(1)+bknotv)*m_basis[bas]->numCoefs(0) + bknotu;
                    for (j3 = m_basis[bas]->order(2)-1; j3 >= 0; j3--)
                      for (j2 = m_basis[bas]->order(1)-1; j2 >= 0; j2--)
                        for (j1 = m_basis[bas]->order(0)-1; j1 >= 0; j1--)
                          myMNPC[iel][lnod++] = iinod - (j3*m_basis[bas]->numCoefs(1)+j2)*m_basis[bas]->numCoefs(0) - j1;
                  } else
                    lnod += m_basis[bas]->order(0)*m_basis[bas]->order(1)*m_basis[bas]->order(2);
                  lnod4 += nb[bas];
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
  for (size_t i = 1; i <= m_basis.size(); ++i) {
    if (basis == 0 || i == (size_t)basis)
      if (!this->connectBasis(face,*neighMx,nface,norient,i,nb1,nb2,coordCheck,thick))
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
  for (size_t i = 0; i < m_basis.size(); nbi += nb[i++])
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

  bool use2ndDer = integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES;
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
  std::vector<std::vector<Go::BasisDerivs2>> splinex2(m_basis.size());
  if (use2ndDer)
#pragma omp parallel for schedule(static)
    for (size_t i=0;i<m_basis.size();++i)
      m_basis[i]->computeBasisGrid(gpar[0],gpar[1],gpar[2],splinex2[i]);
  else
#pragma omp parallel for schedule(static)
    for (size_t i=0;i<m_basis.size();++i)
      m_basis[i]->computeBasisGrid(gpar[0],gpar[1],gpar[2],splinex[i]);

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
  for (size_t g = 0; g < groups.size() && ok; ++g) {
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < groups[g].size(); ++t) {
      MxFiniteElement fe(elem_size);
      std::vector<Matrix> dNxdu(m_basis.size());
      std::vector<Matrix3D> d2Nxdu2(m_basis.size());
      Matrix3D Hess;
      double dXidu[3];
      Matrix Xnod, Jac;
      double   param[3];
      Vec4   X(param);
      for (size_t l = 0; l < groups[g][t].size() && ok; ++l)
      {
        int iel = groups[g][t][l];
        fe.iel = MLGE[iel];
        if (fe.iel < 1) continue; // zero-volume element

        int i1 = p1 + iel % nel1;
        int i2 = p2 + (iel / nel1) % nel2;
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
        if (!integrand.initElement(MNPC[iel-1],elem_size,nb,*A))
        {
          A->destruct();
          ok = false;
          break;
        }


        // --- Integration loop over all Gauss points in each direction --------

        int ip = (((i3-p3)*nGauss*nel2 + i2-p2)*nGauss*nel1 + i1-p1)*nGauss;
        int jp = (((i3-p3)*nel2 + i2-p2)*nel1 + i1-p1)*nGauss*nGauss*nGauss;
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
              fe.u = param[0] = gpar[0](i+1,i1-p1+1);
              fe.v = param[1] = gpar[1](j+1,i2-p2+1);
              fe.w = param[2] = gpar[2](k+1,i3-p3+1);

              // Fetch basis function derivatives at current integration point
              if (use2ndDer)
                for (size_t b = 0; b < m_basis.size(); ++b)
                  SplineUtils::extractBasis(splinex2[b][ip],fe.basis(b+1),dNxdu[b], d2Nxdu2[b]);
              else
                for (size_t b = 0; b < m_basis.size(); ++b)
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
              if (use2ndDer) {
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
              X.assign(Xnod * fe.basis(geoBasis));
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
  const int faceDir = (lIndex%10+1)/((lIndex%2) ? -2 : 2);

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
  for (size_t g = 0; g < threadGrp.size() && ok; ++g) {
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < threadGrp[g].size(); ++t) {
      MxFiniteElement fe(elem_size);
      fe.xi = fe.eta = fe.zeta = faceDir < 0 ? -1.0 : 1.0;
      fe.u = gpar[0](1,1);
      fe.v = gpar[1](1,1);
      fe.w = gpar[2](1,1);

      Matrices dNxdu(m_basis.size());
      Matrix Xnod, Jac;
      double param[3] = { fe.u, fe.v, fe.w };
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


bool ASMs3Dmx::integrate (Integrand& integrand,
                          GlobalIntegral& glInt,
                          const TimeDomain& time,
                          const ASM::InterfaceChecker& iChk)
{
  if (!svol) return true; // silently ignore empty patches
  if (!(integrand.getIntegrandType() & Integrand::INTERFACE_TERMS)) return true;

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
  Vec4          X;
  Vec3          normal;
  double        u[2], v[2], w[2];
  if (MLGE.size() > nel && MLGE.size() != 2*nel) {
    std::cerr << "Interface elements not implemented for mixed integrands." << std::endl;
    return false;
  }

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
        bool ok = integrand.initElement(MNPC[iel],elem_size,nb,*A);
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
            ok &= integrand.initElement(MNPC[kel],elem_size,nb,*A_neigh);
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
              X = Xnod * fe.basis(geoBasis);
              X.t = time.t;

              if (integrand.getIntegrandType() & Integrand::NORMAL_DERIVS)
              {
                std::cerr << "Normal derivs not implemented for mixed integrands." << std::endl;
                return false;
              }

              // Evaluate the integrand and accumulate element contributions
              fe.detJxW *= 0.25*dA*wg[i]*wg[j];
              ok = integrand.evalIntMx(*A,fe,time,X,normal);
            }
          }

        // Finalize the element quantities
        if (ok && !integrand.finalizeElement(*A,time,0))
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
  return this->searchCtrlPt(m_basis[0]->coefs_begin(),m_basis[0]->coefs_end(),
                            X,m_basis[0]->dimension());
}


bool ASMs3Dmx::evalSolution (Matrix& sField, const Vector& locSol,
                             const RealArray* gpar,
                             bool regular, int, int nf) const
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

  std::vector<size_t> nc(nfx.size(), 0);
  if (nf)
    nc[0] = nf;
  else
    std::copy(nfx.begin(), nfx.end(), nc.begin());

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

  // Fetch nodal (control point) coordinates
  Matrix Xnod, Xtmp;
  this->getNodalCoordinates(Xnod);

  MxFiniteElement fe(elem_size,firstIp);
  Vector          solPt;
  Matrices        dNxdu(m_basis.size());
  Matrix          Jac;
  Vec3            X;

  // Evaluate the secondary solution field at each point
  size_t nPoints = splinex[0].size();
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
    if (!integrand.evalSol(solPt,fe,X,ipa,elem_size,nb))
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

  ASMs3D::generateThreadGroups(p[0]-1,p[1]-1,p[2]-1,silence,ignoreGlobalLM);
}


void ASMs3Dmx::generateThreadGroups (char lIndex, bool silence, bool)
{
#ifdef USE_OPENMP
  omp_set_num_threads(1);
#endif
  ASMs3D::generateThreadGroups(lIndex,silence,false);
}


#define DERR -999.99

double ASMs3Dmx::getParametricVolume (int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs3Dmx::getParametricVolume: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return DERR;
  }
#endif
  if (MNPC[iel-1].empty())
    return 0.0;

  int inod1 = MNPC[iel-1][std::accumulate(elem_size.begin(),
                                          elem_size.begin()+geoBasis, -1)];
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
    std::cerr <<" *** ASMs3Dmx::getParametricArea: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return DERR;
  }
#endif
  if (MNPC[iel-1].empty())
    return 0.0;

  int inod1 = MNPC[iel-1][std::accumulate(elem_size.begin(),
                                          elem_size.begin()+geoBasis, -1)];
#ifdef INDEX_CHECK
  if (inod1 < 0 || (size_t)inod1 >= nnod)
  {
    std::cerr <<" *** ASMs3Dmx::getParametricArea: Node index "<< inod1
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

  std::cerr <<" *** ASMs3Dmx::getParametricArea: Invalid face direction "
	    << dir << std::endl;
  return DERR;
}


void ASMs3Dmx::getBoundaryNodes (int lIndex, IntVec& nodes, int basis,
                                 int thick, int, bool local) const
{
  if (basis > 0)
    this->ASMs3D::getBoundaryNodes(lIndex, nodes, basis, thick, 0, local);
  else
    for (size_t b = 1; b <= this->getNoBasis(); ++b)
      this->ASMs3D::getBoundaryNodes(lIndex, nodes, b, thick, 0, local);
}


Fields* ASMs3Dmx::getProjectedFields(const Vector& coefs, size_t nf) const
{
  if (projBasis != m_basis[0])
    return new SplineFields3D(projBasis.get(), coefs, nf);

  return nullptr;
}


size_t ASMs3Dmx::getNoProjectionNodes() const
{
  return projBasis->numCoefs(0) *
         projBasis->numCoefs(1) *
         projBasis->numCoefs(2);
}


bool ASMs3Dmx::evalProjSolution (Matrix& sField, const Vector& locSol,
                                 const int* npe, int nf) const
{
#ifdef SP_DEBUG
  std::cout <<"ASMu3Dmx::evalProjSolution(Matrix&,const Vector&,const int*,int)\n";
#endif
  if (projBasis == m_basis[0])
    return this->evalSolution(sField, locSol, npe, nf);

  // Compute parameter values of the nodal points
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
      return false;

  size_t nComp = locSol.size() / this->getNoProjectionNodes();
  if (nComp*this->getNoProjectionNodes() != locSol.size())
    return false;

  Fields* f = this->getProjectedFields(locSol, nComp);

  // Evaluate the primary solution field at each point
  sField.resize(nComp,gpar[0].size()*gpar[1].size()*gpar[2].size());
  size_t l = 1;
  for (size_t k = 0; k < gpar[2].size(); k++)
    for (size_t j = 0; j < gpar[1].size(); j++)
      for (size_t i = 0; i < gpar[0].size(); i++)
      {
        Vector vals;
        FiniteElement fe;
        fe.u = gpar[0][i];
        fe.v = gpar[1][j];
        fe.w = gpar[2][k];
        f->valueFE(fe, vals);
        sField.fillColumn(l++, vals);
      }

  delete f;

  return true;
}
