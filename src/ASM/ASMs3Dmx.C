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
    return NULL;

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
  if (basis > (int)nb.size())
    basis = 0;

  if (basis == 0)
    return this->ASMbase::getNoNodes(basis);

  return nb[basis-1];
}


unsigned char ASMs3Dmx::getNoFields (int basis) const
{
  if (basis == 0)
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
  size_t nbc=0;
  for (size_t i=0;i<nb.size();++i)
    if (inod <= (nbc+=nb[i]))
      return i == 0 ? 'D' : 'P';

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
  {
    m_basis.resize(2);
    // With mixed methods we need two separate spline spaces
    if (ASMmxBase::Type == FULL_CONT_RAISE_BASIS1 ||
        ASMmxBase::Type == FULL_CONT_RAISE_BASIS2)
    {
      // basis1 should be one degree higher than basis2 and C^p-1 continuous
      int ndim = svol->dimension();
      Go::BsplineBasis b1 = svol->basis(0).extendedBasis(svol->order(0)+1);
      Go::BsplineBasis b2 = svol->basis(1).extendedBasis(svol->order(1)+1);
      Go::BsplineBasis b3 = svol->basis(2).extendedBasis(svol->order(2)+1);
      /* To lower order and regularity this can be used instead
      std::vector<double>::const_iterator first =  ++svol->basis(0).begin();
      std::vector<double>::const_iterator last  =  --svol->basis(0).end();
      Go::BsplineBasis b1 = Go::BsplineBasis(svol->order(0)-1,first,last);
      first =  ++svol->basis(1).begin();
      last  =  --svol->basis(1).end();
      Go::BsplineBasis b2 = Go::BsplineBasis(svol->order(1)-1,first,last);
      first =  ++svol->basis(2).begin();
      last  =  --svol->basis(2).end();
      Go::BsplineBasis b3 = Go::BsplineBasis(svol->order(2)-1,first,last);
      */

      // Note: Currently this is implemented for non-rational splines only.
      // TODO: Ask the splines people how to fix this properly, that is, how
      // may we obtain the correct weights for basis1 when *svol is a NURBS?
      if (svol->rational())
	std::cout <<"WARNING: The geometry basis is rational (using NURBS).\n"
		  <<"         The basis for the unknown fields of one degree"
		  <<" higher will however be non-rational.\n"
		  <<"         This may affect accuracy.\n"<< std::endl;

      // Compute parameter values of the Greville points of the new basis
      size_t i;
      RealArray ug(b1.numCoefs()), vg(b2.numCoefs()), wg(b3.numCoefs());
      for (i = 0; i < ug.size(); i++)
	ug[i] = b1.grevilleParameter(i);
      for (i = 0; i < vg.size(); i++)
	vg[i] = b2.grevilleParameter(i);
      for (i = 0; i < wg.size(); i++)
	wg[i] = b3.grevilleParameter(i);

      // Evaluate the spline volume at all points
      RealArray XYZ(ndim*ug.size()*vg.size()*wg.size());
      svol->gridEvaluator(ug,vg,wg,XYZ);

      // Project the coordinates onto the new basis (the 2nd XYZ is dummy here)
      m_basis[0].reset(Go::VolumeInterpolator::regularInterpolation(b1,b2,b3,
							            ug,vg,wg,XYZ,ndim,
                                                                    false,XYZ));
    }
    else if (ASMmxBase::Type == REDUCED_CONT_RAISE_BASIS1 ||
             ASMmxBase::Type == REDUCED_CONT_RAISE_BASIS2)
    {
      // Order-elevate basis1 such that it is of one degree higher than basis2
      // but only C^p-2 continuous
      m_basis[0].reset(new Go::SplineVolume(*svol));
      m_basis[0]->raiseOrder(1,1,1);
    }
    m_basis[1].reset(new Go::SplineVolume(*svol));

    if (ASMmxBase::Type == FULL_CONT_RAISE_BASIS2 ||
        ASMmxBase::Type == REDUCED_CONT_RAISE_BASIS2) {
      std::swap(m_basis[0], m_basis[1]);
      svol = m_basis[0].get();
    }
  }

  const int n1 = m_basis[0]->numCoefs(0);
  const int n2 = m_basis[0]->numCoefs(1);
  const int n3 = m_basis[0]->numCoefs(2);
  const int m1 = m_basis[1]->numCoefs(0);
  const int m2 = m_basis[1]->numCoefs(1);
  const int m3 = m_basis[1]->numCoefs(2);

  nb.resize(2);
  nb[0] = n1*n2*n3; // Number of functions in first basis
  nb[1] = m1*m2*m3; // Number of functions in second basis

  if (!nodeInd.empty() && !shareFE)
  {
    if (nodeInd.size() == nb[0] + nb[1]) return true;
    std::cerr <<" *** ASMs3Dmx::generateFEMTopology: Inconsistency between the"
	      <<" number of FE nodes "<< nodeInd.size()
	      <<"\n     and the number of spline coefficients "<< nb[0] + nb[1]
	      <<" in the patch."<< std::endl;
    return false;
  }

  if (shareFE == 'F') return true;

  const int p1 = m_basis[0]->order(0);
  const int p2 = m_basis[0]->order(1);
  const int p3 = m_basis[0]->order(2);
  const int q1 = m_basis[1]->order(0);
  const int q2 = m_basis[1]->order(1);
  const int q3 = m_basis[1]->order(2);

#ifdef SP_DEBUG
  std::cout <<"numCoefs: "<< n1 <<" "<< n2 <<" "<< n3
	    <<", "<< m1 <<" "<< m2 <<" "<< m3;
  std::cout <<"\norder: "<< p1 <<" "<< p2 <<" "<< p3
	    <<", "<< q1 <<" "<< q2 <<" "<< q3;
  for (int d = 0; d < 3; d++)
  {
    std::cout <<"\nd"<< char('u'+d) <<':';
    for (int i = 0; i < m_basis[0]->numCoefs(d); i++)
      std::cout <<' '<< m_basis[0]->knotSpan(d,i);
    for (int j = 0; j < m_basis[1]->numCoefs(d); j++)
      std::cout <<' '<< m_basis[1]->knotSpan(d,j);
  }
  std::cout << std::endl;
#endif
  // Consistency checks, just to be fool-proof
  if (m1 <  2 || m2 <  2 || m3 <  2) return false;
  if (q1 <  1 || q2 <  1 || q3 <  1) return false;
  if (p1 > n1 || p2 > n2 || p3 > n3) return false;
  if (q1 > m1 || q2 > m2 || q3 > m3) return false;

  nel = geoBasis == 1 ? (n1-p1+1)*(n2-p2+1)*(n3-p3+1):
                        (m1-q1+1)*(m2-q2+1)*(m3-q3+1);
  nnod = nb[0] + nb[1];

  myMLGE.resize(nel,0);
  myMLGN.resize(nnod);
  myMNPC.resize(nel);
  myNodeInd.resize(nnod);

  int i1, i2, i3, j1, j2, j3;
  size_t iel, inod = 0;
  for (i3 = 0; i3 < n3; i3++)
    for (i2 = 0; i2 < n2; i2++)
      for (i1 = 0; i1 < n1; i1++)
      {
	myNodeInd[inod].I = i1;
	myNodeInd[inod].J = i2;
	myNodeInd[inod].K = i3;
	myMLGN[inod++]    = ++gNod;
      }

  for (i3 = 0; i3 < m3; i3++)
    for (i2 = 0; i2 < m2; i2++)
      for (i1 = 0; i1 < m1; i1++)
      {
	myNodeInd[inod].I = i1;
	myNodeInd[inod].J = i2;
	myNodeInd[inod].K = i3;
	myMLGN[inod++]    = ++gNod;
      }

  if (geoBasis == 1)
  {
    // Create nodal connectivities for basis 1
    iel = inod = 0;
    for (i3 = 1; i3 <= n3; i3++)
      for (i2 = 1; i2 <= n2; i2++)
	for (i1 = 1; i1 <= n1; i1++, inod++)
	  if (i1 >= p1 && i2 >= p2 && i3 >= p3)
	  {
	    if (m_basis[0]->knotSpan(0,i1-1) > 0.0)
	      if (m_basis[0]->knotSpan(1,i2-1) > 0.0)
		if (m_basis[0]->knotSpan(2,i3-1) > 0.0)
		{
		  myMLGE[iel] = ++gEl; // global element number over all patches
		  myMNPC[iel].resize(p1*p2*p3+q1*q2*q3,0);

		  int lnod = 0;
		  for (j3 = p3-1; j3 >= 0; j3--)
		    for (j2 = p2-1; j2 >= 0; j2--)
		      for (j1 = p1-1; j1 >= 0; j1--)
			myMNPC[iel][lnod++] = inod - n1*n2*j3 - n1*j2 - j1;
		}

	    iel++;
	  }

    // Create nodal connectivities for basis 2
    iel = 0;
    for (i3 = 1; i3 <= m3; i3++)
      for (i2 = 1; i2 <= m2; i2++)
	for (i1 = 1; i1 <= m1; i1++, inod++)
	  if (i1 >= q1 && i2 >= q2 && i3 >= q3)
	    if (m_basis[1]->knotSpan(0,i1-1) > 0.0)
	      if (m_basis[1]->knotSpan(1,i2-1) > 0.0)
		if (m_basis[1]->knotSpan(2,i3-1) > 0.0)
		{
		  while (iel < myMNPC.size() && myMNPC[iel].empty()) iel++;

		  int lnod = p1*p2*p3;
		  for (j3 = q3-1; j3 >= 0; j3--)
		    for (j2 = q2-1; j2 >= 0; j2--)
		      for (j1 = q1-1; j1 >= 0; j1--)
			myMNPC[iel][lnod++] = inod - m1*m2*j3 - m1*j2 - j1;

		  iel++;
		}
  }
  else
  {
    // Create nodal connectivities for basis 2
    iel = 0;
    inod = n1*n2*n3;
    for (i3 = 1; i3 <= m3; i3++)
      for (i2 = 1; i2 <= m2; i2++)
	for (i1 = 1; i1 <= m1; i1++, inod++)
	  if (i1 >= q1 && i2 >= q2 && i3 >= q3)
	  {
	    if (m_basis[1]->knotSpan(0,i1-1) > 0.0)
	      if (m_basis[1]->knotSpan(1,i2-1) > 0.0)
		if (m_basis[1]->knotSpan(2,i3-1) > 0.0)
		{
		  myMLGE[iel] = ++gEl; // global element number over all patches
		  myMNPC[iel].resize(p1*p2*p3+q1*q2*q3,0);

		  int lnod = p1*p2*p3;
		  for (j3 = q3-1; j3 >= 0; j3--)
		    for (j2 = q2-1; j2 >= 0; j2--)
		      for (j1 = q1-1; j1 >= 0; j1--)
			myMNPC[iel][lnod++] = inod - m1*m2*j3 - m1*j2 - j1;
		}

	    iel++;
	  }

    // Create nodal connectivities for basis 1
    iel = inod = 0;
    for (i3 = 1; i3 <= n3; i3++)
      for (i2 = 1; i2 <= n2; i2++)
	for (i1 = 1; i1 <= n1; i1++, inod++)
	  if (i1 >= p1 && i2 >= p2 && i3 >= p3)
	    if (m_basis[0]->knotSpan(0,i1-1) > 0.0)
	      if (m_basis[0]->knotSpan(1,i2-1) > 0.0)
		if (m_basis[0]->knotSpan(2,i3-1) > 0.0)
		{
		  while (iel < myMNPC.size() && myMNPC[iel].empty()) iel++;

		  int lnod = 0;
		  for (j3 = p3-1; j3 >= 0; j3--)
		    for (j2 = p2-1; j2 >= 0; j2--)
		      for (j1 = p1-1; j1 >= 0; j1--)
			myMNPC[iel][lnod++] = inod - n1*n2*j3 - n1*j2 - j1;

		  iel++;
		}
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

  if (!this->connectBasis(face,neighbor,nface,norient,1,0,0) ||
      !this->connectBasis(face,neighbor,nface,norient,2,nb[0],neighMx->nb[0]))
    return false;

  this->addNeighbor(neighMx);
  return true;
}


void ASMs3Dmx::closeFaces (int dir, int, int)
{
  this->ASMs3D::closeFaces(dir,1,1);
  this->ASMs3D::closeFaces(dir,2,nb[0]+1);
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
  size_t lnod0 = m_basis[0]->order(0)*m_basis[0]->order(1)*m_basis[0]->order(2);
  if (geoBasis == 1) lnod0 = 0;

  X.resize(3,nenod);
  const IntVec& mnpc = MNPC[iel-1];

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int nd = svol->dimension();
  RealArray::const_iterator cit = svol->coefs_begin();
  for (size_t n = 0; n < nenod; n++)
  {
    int iI = nodeInd[mnpc[lnod0+n]].I;
    int iJ = nodeInd[mnpc[lnod0+n]].J;
    int iK = nodeInd[mnpc[lnod0+n]].K;
    int ip = ((iK*n2 + iJ)*n1 + iI)*nd;
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
  if (inod <= nb[0])
    cit = m_basis[0]->coefs_begin()
      + ((K*m_basis[0]->numCoefs(1)+J)*m_basis[0]->numCoefs(0)+I) * m_basis[0]->dimension();
  else
    cit = m_basis[1]->coefs_begin()
      + ((K*m_basis[1]->numCoefs(1)+J)*m_basis[1]->numCoefs(0)+I) * m_basis[1]->dimension();

  return Vec3(*cit,*(cit+1),*(cit+2));
}


bool ASMs3Dmx::getSize (int& n1, int& n2, int& n3, int basis) const
{
  switch (basis)
    {
    case 1:
      if (m_basis.empty() || !m_basis[0]) return false;
      n1 = m_basis[0]->numCoefs(0);
      n2 = m_basis[0]->numCoefs(1);
      n3 = m_basis[0]->numCoefs(2);
      return true;

    case 2:
      if (m_basis.size() < 2 || ! m_basis[1]) return false;
      n1 = m_basis[1]->numCoefs(0);
      n2 = m_basis[1]->numCoefs(1);
      n3 = m_basis[1]->numCoefs(2);
      return true;
    }

  return this->ASMs3D::getSize(n1,n2,n3);
}


bool ASMs3Dmx::integrate (Integrand& integrand,
			  GlobalIntegral& glInt,
			  const TimeDomain& time)
{
  if (!svol) return true; // silently ignore empty patches
  if (m_basis.empty()) return false;

  PROFILE2("ASMs3Dmx::integrate(I)");

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  std::array<Matrix,3> gpar;
  for (int d = 0; d < 3; d++)
    this->getGaussPointParameters(gpar[d],d,nGauss,xg);

  // Evaluate basis function derivatives at all integration points
  std::vector<Go::BasisDerivs> spline1, spline2;
  m_basis[0]->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline1);
  m_basis[1]->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline2);

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);

  std::vector<size_t> elem_sizes = {(size_t)m_basis[0]->order(0)*m_basis[0]->order(1)*m_basis[0]->order(2),
                                    (size_t)m_basis[1]->order(0)*m_basis[1]->order(1)*m_basis[1]->order(2)};

  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;
  const int nel3 = n3 - p3 + 1;

  // === Assembly loop over all elements in the patch ==========================

  bool ok=true;
  for (size_t g=0;g<threadGroupsVol.size() && ok;++g) {
#pragma omp parallel for schedule(static)
    for (size_t t=0;t<threadGroupsVol[g].size();++t) {
      MxFiniteElement fe(elem_sizes);
      Matrix dN1du, dN2du, Xnod, Jac;
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
              SplineUtils::extractBasis(spline1[ip],fe.basis(1),dN1du);
              SplineUtils::extractBasis(spline2[ip],fe.basis(2),dN2du);

              // Compute Jacobian inverse of the coordinate mapping and
              // basis function derivatives w.r.t. Cartesian coordinates
              if (geoBasis == 1)
              {
                fe.detJxW = utl::Jacobian(Jac,fe.grad(1),Xnod,dN1du);
                fe.grad(2).multiply(dN2du,Jac); // dN2dX = dN2du * J^-1
              }
              else
              {
                fe.detJxW = utl::Jacobian(Jac,fe.grad(2),Xnod,dN2du);
                fe.grad(1).multiply(dN1du,Jac); // dN1dX = dN1du * J^-1
              }
              if (fe.detJxW == 0.0) continue; // skip singular points

              // Cartesian coordinates of current integration point
              X = Xnod * (geoBasis == 1 ? fe.basis(1) : fe.basis(2));
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
  std::vector<Go::BasisDerivs> spline1, spline2;
  m_basis[0]->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline1);
  m_basis[1]->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline2);

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  std::vector<size_t> elem_sizes = {(size_t)m_basis[0]->order(0)*m_basis[0]->order(1)*m_basis[0]->order(2),
                                    (size_t)m_basis[1]->order(0)*m_basis[1]->order(1)*m_basis[1]->order(2)};

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

      Matrix dN1du, dN2du, Xnod, Jac;
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
            SplineUtils::extractBasis(spline1[ip],fe.basis(1),dN1du);
            SplineUtils::extractBasis(spline2[ip],fe.basis(2),dN2du);

            // Compute Jacobian inverse of the coordinate mapping and
            // basis function derivatives w.r.t. Cartesian coordinates
            if (geoBasis == 1)
            {
              fe.detJxW = utl::Jacobian(Jac,normal,fe.grad(1),Xnod,dN1du,t1,t2);
              fe.grad(2).multiply(dN2du,Jac); // dN2dX = dN2du * J^-1
            }
            else
            {
              fe.detJxW = utl::Jacobian(Jac,normal,fe.grad(2),Xnod,dN2du,t1,t2);
              fe.grad(1).multiply(dN1du,Jac); // dN1dX = dN1du * J^-1
            }
            if (fe.detJxW == 0.0) continue; // skip singular points

            if (faceDir < 0) normal *= -1.0;

            // Cartesian coordinates of current integration point
            X = Xnod * (geoBasis == 1 ? fe.basis(1) : fe.basis(2));
            X.t = time.t;

            // Evaluate the integrand and accumulate element contributions
            fe.detJxW *= 0.25*dA*wg[i]*wg[j];
            if (!integrand.evalBouMx(*A,fe,time,X,normal))
              ok = false;
          }

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
  std::vector<Go::BasisPts> spline1, spline2;
  if (regular)
  {
    m_basis[0]->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline1);
    m_basis[1]->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline2);
  }
  else if (gpar[0].size() == gpar[1].size() && gpar[0].size() == gpar[2].size())
  {
    spline1.resize(gpar[0].size());
    spline2.resize(gpar[0].size());
    for (size_t i = 0; i < spline1.size(); i++)
    {
      m_basis[0]->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline1[i]);
      m_basis[1]->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline2[i]);
    }
  }
  else
    return false;

  const int p1 = m_basis[0]->order(0);
  const int p2 = m_basis[0]->order(1);
  const int p3 = m_basis[0]->order(2);
  const int n1 = m_basis[0]->numCoefs(0);
  const int n2 = m_basis[0]->numCoefs(1);
  const int n3 = m_basis[0]->numCoefs(2);

  const int q1 = m_basis[1]->order(0);
  const int q2 = m_basis[1]->order(1);
  const int q3 = m_basis[1]->order(2);
  const int m1 = m_basis[1]->numCoefs(0);
  const int m2 = m_basis[1]->numCoefs(1);
  const int m3 = m_basis[1]->numCoefs(2);

  size_t nc1 = nfx[0];
  size_t nc2 = 0;
  if (nc1*nb[0] < locSol.size())
    nc2 = (locSol.size() - nc1*nb[0])/nb[1];
  else
    nc1 = locSol.size()/nb[0];

  if (nc1*nb[0] + nc2*nb[1] != locSol.size())
    return false;

  Matrix Xtmp;
  Vector Ytmp, Ztmp;

  // Evaluate the primary solution field at each point
  size_t nPoints = spline1.size();
  sField.resize(nc1+nc2,nPoints);
  for (size_t i = 0; i < nPoints; i++)
  {
    IntVec ip;
    scatterInd(n1,n2,n3,p1,p2,p3,spline1[i].left_idx,ip);

    utl::gather(ip,nc1,locSol,Xtmp);
    Xtmp.multiply(spline1[i].basisValues,Ytmp);

    if (nc2 > 0)
    {
      ip.clear();
      scatterInd(m1,m2,m3,q1,q2,q3,spline2[i].left_idx,ip);

      utl::gather(ip,nc2,locSol,Xtmp,nc1*nb[0]);
      Xtmp.multiply(spline2[i].basisValues,Ztmp);

      Ytmp.insert(Ytmp.end(),Ztmp.begin(),Ztmp.end());
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
  std::vector<Go::BasisDerivs> spline1, spline2;
  if (regular)
  {
    m_basis[0]->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline1);
    m_basis[1]->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline2);
  }
  else if (gpar[0].size() == gpar[1].size() && gpar[0].size() == gpar[2].size())
  {
    spline1.resize(gpar[0].size());
    spline2.resize(gpar[0].size());
    for (size_t i = 0; i < spline1.size(); i++)
    {
      m_basis[0]->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline1[i]);
      m_basis[1]->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline2[i]);
    }
  }

  const size_t p1 = m_basis[0]->order(0);
  const size_t p2 = m_basis[0]->order(1);
  const size_t p3 = m_basis[0]->order(2);
  const size_t n1 = m_basis[0]->numCoefs(0);
  const size_t n2 = m_basis[0]->numCoefs(1);
  const size_t n3 = m_basis[0]->numCoefs(2);

  const size_t q1 = m_basis[1]->order(0);
  const size_t q2 = m_basis[1]->order(1);
  const size_t q3 = m_basis[1]->order(2);
  const size_t m1 = m_basis[1]->numCoefs(0);
  const size_t m2 = m_basis[1]->numCoefs(1);
  const size_t m3 = m_basis[1]->numCoefs(2);
  std::vector<size_t> elem_sizes = {p1*p2*p3, q1*q2*q3};

  // Fetch nodal (control point) coordinates
  Matrix Xnod, Xtmp;
  this->getNodalCoordinates(Xnod);

  MxFiniteElement fe(elem_sizes);
  Vector          solPt;
  Matrix          dN1du, dN2du, Jac;
  Vec3            X;

  // Evaluate the secondary solution field at each point
  size_t nPoints = spline1.size();
  for (size_t i = 0; i < nPoints; i++)
  {
    // Fetch indices of the non-zero basis functions at this point
    IntVec ip1, ip2;
    scatterInd(n1,n2,n3,p1,p2,p3,spline1[i].left_idx,ip1);
    scatterInd(m1,m2,m3,q1,q2,q3,spline2[i].left_idx,ip2);
    ip1.insert(ip1.end(), ip2.begin(), ip2.end());
    fe.u = spline1[i].param[0];
    fe.v = spline1[i].param[1];
    fe.w = spline1[i].param[2];
    fe.iGP = firstIp + i;

    // Fetch associated control point coordinates
    utl::gather(geoBasis == 1 ? ip1 : ip2, 3, Xnod, Xtmp);

    // Fetch basis function derivatives at current integration point
    SplineUtils::extractBasis(spline1[i],fe.basis(1),dN1du);
    SplineUtils::extractBasis(spline2[i],fe.basis(2),dN2du);

    // Compute Jacobian inverse of the coordinate mapping and
    // basis function derivatives w.r.t. Cartesian coordinates
    if (geoBasis == 1)
      if (utl::Jacobian(Jac,fe.grad(1),Xtmp,dN1du) == 0.0) // Jac = (Xtmp*dN1du)^-1
	continue; // skip singular points
      else
	fe.grad(2).multiply(dN2du,Jac); // dN2dX = dN2du * J^-1
    else
      if (utl::Jacobian(Jac,fe.grad(2),Xtmp,dN2du) == 0.0) // Jac = (Xtmp*dN2du)^-1
	continue; // skip singular points
      else
	fe.grad(1).multiply(dN1du,Jac); // dN1dX = dN1du * J^-1

    // Cartesian coordinates of current integration point
    X = Xtmp * (geoBasis == 1 ? fe.basis(1) : fe.basis(2));

    // Now evaluate the solution field
    if (!integrand.evalSol(solPt,fe,X,ip1,elem_sizes))
      return false;
    else if (sField.empty())
      sField.resize(solPt.size(),nPoints,true);

    sField.fillColumn(1+i,solPt);
  }

  return true;
}


void ASMs3Dmx::generateThreadGroups (const Integrand& integrand, bool silence)
{
#ifdef USE_OPENMP
  omp_set_num_threads(1);
#endif
  ASMs3D::generateThreadGroups(integrand,silence);
}


void ASMs3Dmx::generateThreadGroups (char lIndex, bool silence)
{
#ifdef USE_OPENMP
  omp_set_num_threads(1);
#endif
  ASMs3D::generateThreadGroups(lIndex,silence);
}
