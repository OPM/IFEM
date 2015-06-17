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
  : ASMs2D(patch), ASMmxBase(n_f[0]==0?patch.nfx:n_f)
{
  nb = patch.nb;
}


Go::SplineSurface* ASMs2Dmx::getBasis (int basis) const
{
  return basis == 2 ? m_basis[1].get() : m_basis[0].get();
}


Go::SplineCurve* ASMs2Dmx::getBoundary (int dir, int basis)
{
  if (dir < -2 || dir == 0 || dir > 2 || basis > (int)m_basis.size())
    return NULL;

  int iedge = dir > 0 ? dir : 3*dir+6;

  return m_basis[basis-1]->edgeCurve(iedge);
}


bool ASMs2Dmx::write (std::ostream& os, int basis) const
{
  if (m_basis[0] && basis == 1)
    os <<"200 1 0 0\n" << *m_basis[0];
  else if (m_basis[1] && basis == 2)
    os <<"200 1 0 0\n" << *m_basis[1];
  else if (surf)
    os <<"200 1 0 0\n" << *surf;
  else
    return false;

  return os.good();
}


void ASMs2Dmx::clear (bool retainGeometry)
{
  // Erase the spline data
  if (!retainGeometry)
    surf = 0;

  if (!shareFE)
    m_basis[0].reset();
  if (!shareFE)
    m_basis[1].reset();

  // Erase the FE data
  this->ASMs2D::clear(retainGeometry);
}


size_t ASMs2Dmx::getNoNodes (int basis) const
{
  if (basis > (int)nb.size())
    basis = 0;

  if (basis == 0)
    return this->ASMbase::getNoNodes(basis);

  return nb[basis-1];
}


unsigned char ASMs2Dmx::getNoFields (int basis) const
{
  if (basis == 0)
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
  size_t nbc=0;
  for (size_t i=0;i<nb.size();++i)
    if (inod <= (nbc+=nb[i]))
      return i == 0 ? 'D' : 'P';

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
  {
    m_basis.resize(2);
    // With mixed methods we need two separate spline spaces
    if (ASMmxBase::Type == ASMmxBase::FULL_CONT_RAISE_BASIS1 ||
        ASMmxBase::Type == ASMmxBase::FULL_CONT_RAISE_BASIS2)
    {
      // basis1 should be one degree higher than basis2 and C^p-1 continuous
      int ndim = surf->dimension();
      Go::BsplineBasis b1 = surf->basis(0).extendedBasis(surf->order_u()+1);
      Go::BsplineBasis b2 = surf->basis(1).extendedBasis(surf->order_v()+1);
      /* To lower order and regularity this can be used instead
      std::vector<double>::const_iterator first = ++surf->basis(0).begin();
      std::vector<double>::const_iterator last  = --surf->basis(0).end();
      Go::BsplineBasis b1 = Go::BsplineBasis(surf->order_u()-1,first,last);
      first =  ++surf->basis(1).begin();
      last  =  --surf->basis(1).end();
      Go::BsplineBasis b2 = Go::BsplineBasis(surf->order_v()-1,first,last);
      */

      // Note: Currently this is implemented for non-rational splines only.
      // TODO: Ask the splines people how to fix this properly, that is, how
      // may we obtain the correct weights for basis1 when *surf is a NURBS?
      if (surf->rational())
	std::cout <<"WARNING: The geometry basis is rational (using NURBS).\n"
		  <<"         The basis for the unknown fields of one degree"
		  <<" higher will however be non-rational.\n"
		  <<"         This may affect accuracy.\n"<< std::endl;

      // Compute parameter values of the Greville points
      size_t i;
      RealArray ug(b1.numCoefs()), vg(b2.numCoefs());
      for (i = 0; i < ug.size(); i++)
	ug[i] = b1.grevilleParameter(i);
      for (i = 0; i < vg.size(); i++)
	vg[i] = b2.grevilleParameter(i);

      // Evaluate the spline surface at all points
      RealArray XYZ(ndim*ug.size()*vg.size());
      surf->gridEvaluator(XYZ,ug,vg);

      // Project the coordinates onto the new basis (the 2nd XYZ is dummy here)
      m_basis[0].reset(Go::SurfaceInterpolator::regularInterpolation(b1,b2,
                                                                     ug,vg,XYZ,ndim,
                                                                     false,XYZ));
    }
    else if (ASMmxBase::Type == REDUCED_CONT_RAISE_BASIS1 ||
             ASMmxBase::Type == REDUCED_CONT_RAISE_BASIS2)
    {
      // Order-elevate basis1 such that it is of one degree higher than basis2
      // but only C^p-2 continuous
      m_basis[0].reset(new Go::SplineSurface(*surf));
      m_basis[0]->raiseOrder(1,1);
    }
    m_basis[1].reset(new Go::SplineSurface(*surf));

    if (ASMmxBase::Type == FULL_CONT_RAISE_BASIS2 ||
        ASMmxBase::Type == REDUCED_CONT_RAISE_BASIS2) {
      std::swap(m_basis[0], m_basis[1]);
      surf = m_basis[0].get();
    }
  }

  const int n1 = m_basis[0]->numCoefs_u();
  const int n2 = m_basis[0]->numCoefs_v();
  const int m1 = m_basis[1]->numCoefs_u();
  const int m2 = m_basis[1]->numCoefs_v();

  nb.resize(2);
  nb[0] = n1*n2; // Number of functions in first basis
  nb[1] = m1*m2; // Number of functions in second basis

  if (!nodeInd.empty() && !shareFE)
  {
    if (nodeInd.size() == nb[0] + nb[1]) return true;
    std::cerr <<" *** ASMs2Dmx::generateFEMTopology: Inconsistency between the"
	      <<" number of FE nodes "<< nodeInd.size()
	      <<"\n     and the number of spline coefficients "<< nb[0] + nb[1]
	      <<" in the patch."<< std::endl;
    return false;
  }

  if (shareFE == 'F') return true;

  const int p1 = m_basis[0]->order_u();
  const int p2 = m_basis[0]->order_v();
  const int q1 = m_basis[1]->order_u();
  const int q2 = m_basis[1]->order_v();

  int i1, i2, j1, j2;
#ifdef SP_DEBUG
  std::cout <<"numCoefs: "<< n1 <<" "<< n2 <<", "<< m1 <<" "<< m2;
  std::cout <<"\norder: "<< p1 <<" "<< p2 <<", "<< q1 <<" "<< q2;
  std::cout <<"\ndu:";
  for (i1 = 0; i1 < n1; i1++)
    std::cout <<' '<< basis1->knotSpan(0,i1);
  for (i1 = 0; i1 < m1; i1++)
    std::cout <<' '<< basis2->knotSpan(0,i1);
  std::cout <<"\ndv:";
  for (i2 = 0; i2 < n2; i2++)
    std::cout <<' '<< basis1->knotSpan(1,i2);
  for (i2 = 0; i2 < m2; i2++)
    std::cout <<' '<< basis2->knotSpan(1,i2);
  std::cout << std::endl;
#endif
  // Consistency checks, just to be fool-proof
  if (m1 <  2 || m2 <  2) return false;
  if (q1 <  1 || q2 <  1) return false;
  if (p1 > n1 || p2 > n2) return false;
  if (q1 > m1 || q2 > m2) return false;

  nel = geoBasis == 1 ? (n1-p1+1)*(n2-p2+1) : (m1-q1+1)*(m2-q2+1);
  nnod = nb[0] + nb[1];

  myMLGE.resize(nel,0);
  myMLGN.resize(nnod);
  myMNPC.resize(nel);
  myNodeInd.resize(nnod);

  size_t iel, inod = 0;
  for (i2 = 0; i2 < n2; i2++)
    for (i1 = 0; i1 < n1; i1++)
    {
      myNodeInd[inod].I = i1;
      myNodeInd[inod].J = i2;
      myMLGN[inod++]    = ++gNod;
    }

  for (i2 = 0; i2 < m2; i2++)
    for (i1 = 0; i1 < m1; i1++)
    {
      myNodeInd[inod].I = i1;
      myNodeInd[inod].J = i2;
      myMLGN[inod++]    = ++gNod;
    }

  if (geoBasis == 1)
  {
    // Create nodal connectivities for basis 1
    iel = inod = 0;
    for (i2 = 1; i2 <= n2; i2++)
      for (i1 = 1; i1 <= n1; i1++, inod++)
	if (i1 >= p1 && i2 >= p2)
	{
	  if (m_basis[0]->knotSpan(0,i1-1) > 0.0)
	    if (m_basis[0]->knotSpan(1,i2-1) > 0.0)
	    {
	      myMLGE[iel] = ++gEl; // global element number over all patches
	      myMNPC[iel].resize(p1*p2+q1*q2,0);

	      int lnod = 0;
	      for (j2 = p2-1; j2 >= 0; j2--)
		for (j1 = p1-1; j1 >= 0; j1--)
		  myMNPC[iel][lnod++] = inod - n1*j2 - j1;
	    }

	  iel++;
	}

    // Create nodal connectivities for basis 2
    iel = 0;
    for (i2 = 1; i2 <= m2; i2++)
      for (i1 = 1; i1 <= m1; i1++, inod++)
	if (i1 >= q1 && i2 >= q2)
	  if (m_basis[1]->knotSpan(0,i1-1) > 0.0)
	    if (m_basis[1]->knotSpan(1,i2-1) > 0.0)
	    {
	      while (iel < myMNPC.size() && myMNPC[iel].empty()) iel++;

	      int lnod = p1*p2;
	      for (j2 = q2-1; j2 >= 0; j2--)
		for (j1 = q1-1; j1 >= 0; j1--)
		  myMNPC[iel][lnod++] = inod - m1*j2 - j1;

	      iel++;
	    }
  }
  else
  {
    // Create nodal connectivities for basis 2
    iel = 0;
    inod = n1*n2;
    for (i2 = 1; i2 <= m2; i2++)
      for (i1 = 1; i1 <= m1; i1++, inod++)
	if (i1 >= q1 && i2 >= q2)
	{
	  if (m_basis[1]->knotSpan(0,i1-1) > 0.0)
	    if (m_basis[1]->knotSpan(1,i2-1) > 0.0)
	    {
	      myMLGE[iel] = ++gEl; // global element number over all patches
	      myMNPC[iel].resize(p1*p2+q1*q2,0);

	      int lnod = p1*p2;
	      for (j2 = q2-1; j2 >= 0; j2--)
		for (j1 = q1-1; j1 >= 0; j1--)
		  myMNPC[iel][lnod++] = inod - m1*j2 - j1;
	    }

	  iel++;
	}

    // Create nodal connectivities for basis 1
    iel = inod = 0;
    for (i2 = 1; i2 <= n2; i2++)
      for (i1 = 1; i1 <= n1; i1++, inod++)
	if (i1 >= p1 && i2 >= p2)
	{
	  if (m_basis[0]->knotSpan(0,i1-1) > 0.0)
	    if (m_basis[0]->knotSpan(1,i2-1) > 0.0)
	    {
	      while (iel < myMNPC.size() && myMNPC[iel].empty()) iel++;

	      int lnod = 0;
	      for (j2 = p2-1; j2 >= 0; j2--)
		for (j1 = p1-1; j1 >= 0; j1--)
		  myMNPC[iel][lnod++] = inod - n1*j2 - j1;

	      iel++;
	    }	
	}
  }

#ifdef SP_DEBUG
  std::cout <<"NEL = "<< nel <<" NNOD = "<< nnod << std::endl;
#endif
  return true;
}


bool ASMs2Dmx::connectPatch (int edge, ASMs2D& neighbor, int nedge, bool revers)
{
  ASMs2Dmx* neighMx = dynamic_cast<ASMs2Dmx*>(&neighbor);
  if (!neighMx) return false;

  if (!this->connectBasis(edge,neighbor,nedge,revers,1,0,0) ||
      !this->connectBasis(edge,neighbor,nedge,revers,2,nb[0],neighMx->nb[0]))
    return false;

  this->addNeighbor(neighMx);
  return true;
}


void ASMs2Dmx::closeEdges (int dir, int, int)
{
  this->ASMs2D::closeEdges(dir,1,1);
  this->ASMs2D::closeEdges(dir,2,nb[0]+1);
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
  size_t lnod0 = geoBasis == 1 ? 0 : m_basis[0]->order_u()*m_basis[0]->order_v();

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
  if (inod <= nb[0])
    cit = m_basis[0]->coefs_begin()
      + (J*m_basis[0]->numCoefs_u()+I) * m_basis[0]->dimension();
  else
    cit = m_basis[1]->coefs_begin()
      + (J*m_basis[1]->numCoefs_u()+I) * m_basis[1]->dimension();

  Vec3 X;
  for (size_t i = 0; i < nsd; i++, cit++)
    X[i] = *cit;

  return X;
}


bool ASMs2Dmx::getSize (int& n1, int& n2, int basis) const
{
  switch (basis)
    {
    case 1:
      n1 = m_basis[0]->numCoefs_u();
      n2 = m_basis[0]->numCoefs_v();
      return true;

    case 2:
      n1 = m_basis[1]->numCoefs_u();
      n2 = m_basis[1]->numCoefs_v();
      return true;
    }

  return this->ASMs2D::getSize(n1,n2);
}


bool ASMs2Dmx::integrate (Integrand& integrand,
			  GlobalIntegral& glInt,
			  const TimeDomain& time)
{
  if (!surf) return true; // silently ignore empty patches

  PROFILE2("ASMs2Dmx::integrate(I)");

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  std::array<Matrix,2> gpar;
  for (int d = 0; d < 2; d++)
    this->getGaussPointParameters(gpar[d],d,nGauss,xg);

  // Evaluate basis function derivatives at all integration points
  std::vector<Go::BasisDerivsSf> spline1, spline2;
  m_basis[0]->computeBasisGrid(gpar[0],gpar[1],spline1);
  m_basis[1]->computeBasisGrid(gpar[0],gpar[1],spline2);

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int n1 = surf->numCoefs_u();
  const int nel1 = n1 - p1 + 1;
  std::vector<size_t> elem_sizes = {(size_t)m_basis[0]->order_u()*m_basis[0]->order_v(),
                                    (size_t)m_basis[1]->order_u()*m_basis[1]->order_v()};

  // === Assembly loop over all elements in the patch ==========================

  bool ok=true;
  for (size_t g=0;g<threadGroups.size() && ok;++g) {
#pragma omp parallel for schedule(static)
    for (size_t t=0;t<threadGroups[g].size();++t) {
      MxFiniteElement fe(elem_sizes);
      Matrix dN1du, dN2du, Xnod, Jac;
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

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Find the parametric direction of the edge normal {-2,-1, 1, 2}
  const int edgeDir = (lIndex+1)/(lIndex%2 ? -2 : 2);

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

  // Evaluate basis function derivatives at all integration points
  std::vector<Go::BasisDerivsSf> spline1, spline2;
  m_basis[0]->computeBasisGrid(gpar[0],gpar[1],spline1);
  m_basis[1]->computeBasisGrid(gpar[0],gpar[1],spline2);

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();

  std::vector<size_t> elem_sizes = {(size_t)m_basis[0]->order_u()*m_basis[0]->order_v(),
                                    (size_t)m_basis[1]->order_u()*m_basis[1]->order_v()};

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;

  MxFiniteElement fe(elem_sizes);
  fe.xi = fe.eta = edgeDir < 0 ? -1.0 : 1.0;
  fe.u = gpar[0](1,1);
  fe.v = gpar[1](1,1);

  Matrix dN1du, dN2du, Xnod, Jac;
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

	if (edgeDir < 0) normal *= -1.0;

	// Cartesian coordinates of current integration point
	X = Xnod * (geoBasis == 1 ? fe.basis(1) : fe.basis(2));
	X.t = time.t;

	// Evaluate the integrand and accumulate element contributions
	fe.detJxW *= 0.5*dS*wg[i];

	ok = integrand.evalBouMx(*A,fe,time,X,normal);
      }

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
  std::vector<Go::BasisPtsSf> spline1, spline2;
  if (regular)
  {
    m_basis[0]->computeBasisGrid(gpar[0],gpar[1],spline1);
    m_basis[1]->computeBasisGrid(gpar[0],gpar[1],spline2);
  }
  else if (gpar[0].size() == gpar[1].size())
  {
    spline1.resize(gpar[0].size());
    spline2.resize(gpar[0].size());
    for (size_t i = 0; i < spline1.size(); i++)
    {
      m_basis[0]->computeBasis(gpar[0][i],gpar[1][i],spline1[i]);
      m_basis[1]->computeBasis(gpar[0][i],gpar[1][i],spline2[i]);
    }
  }
  else
    return false;

  const int p1 = m_basis[0]->order_u();
  const int p2 = m_basis[0]->order_v();
  const int n1 = m_basis[0]->numCoefs_u();
  const int n2 = m_basis[0]->numCoefs_v();

  const int q1 = m_basis[1]->order_u();
  const int q2 = m_basis[1]->order_v();
  const int m1 = m_basis[1]->numCoefs_u();
  const int m2 = m_basis[1]->numCoefs_v();

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
    scatterInd(n1,n2,p1,p2,spline1[i].left_idx,ip);

    utl::gather(ip,nc1,locSol,Xtmp);
    Xtmp.multiply(spline1[i].basisValues,Ytmp);

    if (nc2 > 0)
    {
      ip.clear();
      scatterInd(m1,m2,q1,q2,spline2[i].left_idx,ip);

      utl::gather(ip,nc2,locSol,Xtmp,nc1*nb[0]);
      Xtmp.multiply(spline2[i].basisValues,Ztmp);

      Ytmp.insert(Ytmp.end(),Ztmp.begin(),Ztmp.end());
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
  std::vector<Go::BasisDerivsSf> spline1, spline2;
  if (regular)
  {
    m_basis[0]->computeBasisGrid(gpar[0],gpar[1],spline1);
    m_basis[1]->computeBasisGrid(gpar[0],gpar[1],spline2);
  }
  else if (gpar[0].size() == gpar[1].size())
  {
    spline1.resize(gpar[0].size());
    spline2.resize(gpar[0].size());
    for (size_t i = 0; i < spline1.size(); i++)
    {
      m_basis[0]->computeBasis(gpar[0][i],gpar[1][i],spline1[i]);
      m_basis[1]->computeBasis(gpar[0][i],gpar[1][i],spline2[i]);
    }
  }

  const size_t p1 = m_basis[0]->order_u();
  const size_t p2 = m_basis[0]->order_v();
  const size_t n1 = m_basis[0]->numCoefs_u();
  const size_t n2 = m_basis[0]->numCoefs_v();

  const size_t q1 = m_basis[1]->order_u();
  const size_t q2 = m_basis[1]->order_v();
  const size_t m1 = m_basis[1]->numCoefs_u();
  const size_t m2 = m_basis[1]->numCoefs_v();
  std::vector<size_t> elem_sizes = {p1*p2,q1*q2};

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
    scatterInd(n1,n2,p1,p2,spline1[i].left_idx,ip1);
    scatterInd(m1,m2,q1,q2,spline2[i].left_idx,ip2);
    ip1.insert(ip1.end(), ip2.begin(), ip2.end());
    fe.u = spline1[i].param[0];
    fe.v = spline1[i].param[1];
    fe.iGP = firstIp + i;

    // Fetch associated control point coordinates
    utl::gather(geoBasis == 1 ? ip1 : ip2, nsd, Xnod, Xtmp);

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


void ASMs2Dmx::generateThreadGroups (const Integrand& integrand, bool silence)
{
#ifdef USE_OPENMP
  omp_set_num_threads(1);
#endif
  ASMs2D::generateThreadGroups(integrand, silence);
}
