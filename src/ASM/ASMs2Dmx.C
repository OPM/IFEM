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
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "Utilities.h"
#include "Profiler.h"
#include "Vec3Oper.h"
#include "Vec3.h"


ASMs2Dmx::ASMs2Dmx (unsigned char n_s,
		    unsigned char n_f1, unsigned char n_f2)
  : ASMs2D(n_s), ASMmxBase(n_f1,n_f2)
{
  basis1 = basis2 = 0;
  nf = nf1 + nf2;
}


ASMs2Dmx::ASMs2Dmx (const ASMs2Dmx& patch, char n_f1, char n_f2)
  : ASMs2D(patch), ASMmxBase(patch.nf1,patch.nf2)
{
  basis1 = patch.basis1;
  basis2 = patch.basis2;
  nb1 = patch.nb1;
  nb2 = patch.nb2;
  if (n_f1 >= 0) nf1 = n_f1;
  if (n_f2 >= 0) nf2 = n_f2;
  nf = nf1 + nf2;
}


Go::SplineSurface* ASMs2Dmx::getBasis(int basis) const
{
  if (basis == 2) 
    return basis2;
  else
    return basis1;
}


bool ASMs2Dmx::write (std::ostream& os, int basis) const
{
  if (basis1 && basis == 1)
    os <<"200 1 0 0\n" << *basis1;
  else if (basis2 && basis == 2)
    os <<"200 1 0 0\n" << *basis2;
  else if (surf)
    os <<"200 1 0 0\n" << *surf;
  else
    return false;

  return os.good();
}


void ASMs2Dmx::clear (bool retainGeometry)
{
  // Erase the spline data
  if (retainGeometry)
  {
    if (basis1 && basis1 != surf && !shareFE) delete basis1;
    if (basis2 && basis2 != surf && !shareFE) delete basis2;
  }
  else
  {
    if (basis1 && !shareFE) delete basis1;
    if (basis2 && !shareFE) delete basis2;
    surf = 0;
  }
  basis1 = basis2 = 0;

  // Erase the FE data
  this->ASMs2D::clear(retainGeometry);
}


size_t ASMs2Dmx::getNoNodes (int basis) const
{
  switch (basis)
    {
    case 1: return nb1;
    case 2: return nb2;
    }

  return nb1+nb2;
}


unsigned char ASMs2Dmx::getNoFields (int basis) const
{
  switch (basis)
    {
    case 1: return nf1;
    case 2: return nf2;
    }

  return nf;
}


unsigned char ASMs2Dmx::getNodalDOFs (size_t inod) const
{
  return inod <= nb1 ? nf1 : nf2;
}


unsigned char ASMs2Dmx::getNodalBasis (size_t inod) const
{
  return inod <= nb1 ? 1 : 2;
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


bool ASMs2Dmx::getSolution (Matrix& sField, const Vector& locSol,
			    const IntVec& nodes) const
{
  return this->getSolutionMx(sField,locSol,nodes);
}


bool ASMs2Dmx::generateFEMTopology ()
{
  if (!surf) return false;

  if (!basis1 && !basis2)
  {
    // With mixed methods we need two separate spline spaces
    if (useCpminus1)
    {
      // basis1 should be one degree higher than basis2 and C^p-1 continuous
      int ndim = surf->dimension();
      Go::BsplineBasis b1 = surf->basis(0).extendedBasis(surf->order_u()+1);
      Go::BsplineBasis b2 = surf->basis(1).extendedBasis(surf->order_v()+1);
      // To lower order and regularity this can be used
//       std::vector<double>::const_iterator first =  ++surf->basis(0).begin();
//       std::vector<double>::const_iterator last  =  --surf->basis(0).end();
//       Go::BsplineBasis b1 = Go::BsplineBasis(surf->order_u()-1,first,last);
//       first =  ++surf->basis(1).begin();
//       last  =  --surf->basis(1).end();
//       Go::BsplineBasis b2 = Go::BsplineBasis(surf->order_v()-1,first,last);

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
      basis1 = Go::SurfaceInterpolator::regularInterpolation(b1,b2,
							     ug,vg,XYZ,ndim,
							     false,XYZ);
    }
    else
    {
      // Order-elevate basis1 such that it is of one degree higher than basis2
      // but only C^p-2 continuous
      basis1 = new Go::SplineSurface(*surf);
      basis1->raiseOrder(1,1);
    }
    basis2 = surf;

    // Define which basis that should be used to represent the geometry
    if (geoUsesBasis1) surf = basis1;
  }

  const int n1 = basis1->numCoefs_u();
  const int n2 = basis1->numCoefs_v();
  const int m1 = basis2->numCoefs_u();
  const int m2 = basis2->numCoefs_v();

  if (!nodeInd.empty() && !shareFE)
  {
    if (nodeInd.size() == nb1 + nb2) return true;
    std::cerr <<" *** ASMs2Dmx::generateFEMTopology: Inconsistency between the"
	      <<" number of FE nodes "<< nodeInd.size()
	      <<"\n     and the number of spline coefficients "<< nb1 + nb2
	      <<" in the patch."<< std::endl;
    return false;
  }

  nb1 = n1*n2; // Number of functions in first basis
  nb2 = m1*m2; // Number of functions in second basis

  if (shareFE) return true;

  const int p1 = basis1->order_u();
  const int p2 = basis1->order_v();
  const int q1 = basis2->order_u();
  const int q2 = basis2->order_v();

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

  myMLGE.resize((n1-p1+1)*(n2-p2+1),0);
  myMLGN.resize(nb1 + nb2);
  myMNPC.resize(myMLGE.size());
  myNodeInd.resize(myMLGN.size());

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

  if (geoUsesBasis1)
  {
    // Create nodal connectivities for basis 1
    iel = inod = 0;
    for (i2 = 1; i2 <= n2; i2++)
      for (i1 = 1; i1 <= n1; i1++, inod++)
	if (i1 >= p1 && i2 >= p2)
	{
	  if (basis1->knotSpan(0,i1-1) > 0.0)
	    if (basis1->knotSpan(1,i2-1) > 0.0)
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
	  if (basis2->knotSpan(0,i1-1) > 0.0)
	    if (basis2->knotSpan(1,i2-1) > 0.0)
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
	  if (basis2->knotSpan(0,i1-1) > 0.0)
	    if (basis2->knotSpan(1,i2-1) > 0.0)
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
	  if (basis1->knotSpan(0,i1-1) > 0.0)
	    if (basis1->knotSpan(1,i2-1) > 0.0)
	    {
	      while (iel < myMNPC.size() && myMNPC[iel].empty()) iel++;

	      int lnod = 0;
	      for (j2 = p2-1; j2 >= 0; j2--)
		for (j1 = p1-1; j1 >= 0; j1--) 
		  myMNPC[iel][lnod++] = inod - n1*j2 - j1;

	      iel++;
	    }
  }

#ifdef SP_DEBUG
  std::cout <<"NEL = "<< MLGE.size() <<" NNOD = "<< MLGN.size() << std::endl;
#endif
  return true;
}


bool ASMs2Dmx::connectPatch (int edge, ASMs2D& neighbor, int nedge, bool revers)
{
  ASMs2Dmx* neighMx = dynamic_cast<ASMs2Dmx*>(&neighbor);
  if (!neighMx) return false;

  return this->connectBasis(edge,neighbor,nedge,revers,1,0,0)
    &&   this->connectBasis(edge,neighbor,nedge,revers,2,nb1,neighMx->nb1);
}


void ASMs2Dmx::closeEdges (int dir, int, int)
{
  this->ASMs2D::closeEdges(dir,1,1);
  this->ASMs2D::closeEdges(dir,2,nb1+1);
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
  size_t lnod0 = geoUsesBasis1 ? 0 : basis1->order_u()*basis1->order_v();

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
  if (inod <= nb1)
    cit = basis1->coefs_begin()
      + (J*basis1->numCoefs_u()+I) * basis1->dimension();
  else
    cit = basis2->coefs_begin()
      + (J*basis2->numCoefs_u()+I) * basis2->dimension();

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
      if (!basis1) return false;
      n1 = basis1->numCoefs_u();
      n2 = basis1->numCoefs_v();
      return true;

    case 2:
      if (!basis2) return false;
      n1 = basis2->numCoefs_u();
      n2 = basis2->numCoefs_v();
      return true;
    }

  return this->ASMs2D::getSize(n1,n2);
}


bool ASMs2Dmx::integrate (Integrand& integrand,
			  GlobalIntegral& glInt,
			  const TimeDomain& time,
			  const LintegralVec& locInt)
{
  if (!surf) return true; // silently ignore empty patches
  if (!basis1 || !basis2) return false;

  PROFILE2("ASMs2Dmx::integrate(I)");

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  Matrix gpar[2];
  for (int d = 0; d < 2; d++)
    this->getGaussPointParameters(gpar[d],d,nGauss,xg);

  // Evaluate basis function derivatives at all integration points
  std::vector<Go::BasisDerivsSf> spline1, spline2;
  basis1->computeBasisGrid(gpar[0],gpar[1],spline1);
  basis2->computeBasisGrid(gpar[0],gpar[1],spline2);

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();
  const int nel1 = n1 - p1 + 1;

  MxFiniteElement fe(basis1->order_u()*basis1->order_v(),
		     basis2->order_u()*basis2->order_v());
  Matrix dN1du, dN2du, Xnod, Jac;
  Vec4   X;

  // === Assembly loop over all elements in the patch ==========================

  int iel = 1;
  for (int i2 = p2; i2 <= n2; i2++)
    for (int i1 = p1; i1 <= n1; i1++, iel++)
    {
      fe.iel = MLGE[iel-1];
      if (fe.iel < 1) continue; // zero-area element

      // Get element area in the parameter space
      double dA = this->getParametricArea(iel);
      if (dA < 0.0) return false; // topology error (probably logic error)

      // Set up control point (nodal) coordinates for current element
      if (!this->getElementCoordinates(Xnod,iel)) return false;

      // Initialize element quantities
      IntVec::const_iterator f2start = MNPC[iel-1].begin() + fe.N1.size();
      if (!integrand.initElement(IntVec(MNPC[iel-1].begin(),f2start),
				 IntVec(f2start,MNPC[iel-1].end()),nb1))
	return false;

      // Caution: Unless locInt is empty, we assume it points to an array of
      // LocalIntegral pointers, of length at least the number of elements in
      // the model (as defined by the highest number in the MLGE array).
      // If the array is shorter than this, expect a segmentation fault.
      LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[fe.iel-1];


      // --- Integration loop over all Gauss points in each direction ----------

      int ip = ((i2-p2)*nGauss*nel1 + i1-p1)*nGauss;
      for (int j = 0; j < nGauss; j++, ip += nGauss*(nel1-1))
	for (int i = 0; i < nGauss; i++, ip++)
	{
          // Local element coordinates of current integration point
          fe.xi  = xg[i];
          fe.eta = xg[j];

	  // Parameter values of current integration point
	  fe.u = gpar[0](i+1,i1-p1+1);
	  fe.v = gpar[1](j+1,i2-p2+1);

	  // Fetch basis function derivatives at current integration point
	  extractBasis(spline1[ip],fe.N1,dN1du);
	  extractBasis(spline2[ip],fe.N2,dN2du);

	  // Compute Jacobian inverse of the coordinate mapping and
	  // basis function derivatives w.r.t. Cartesian coordinates
	  if (geoUsesBasis1)
	  {
	    fe.detJxW = utl::Jacobian(Jac,fe.dN1dX,Xnod,dN1du);
	    fe.dN2dX.multiply(dN2du,Jac); // dN2dX = dN2du * J^-1
	  }
	  else
	  {
	    fe.detJxW = utl::Jacobian(Jac,fe.dN2dX,Xnod,dN2du);
	    fe.dN1dX.multiply(dN1du,Jac); // dN1dX = dN1du * J^-1
	  }
	  if (fe.detJxW == 0.0) continue; // skip singular points

	  // Cartesian coordinates of current integration point
	  X = Xnod * (geoUsesBasis1 ? fe.N1 : fe.N2);
	  X.t = time.t;

	  // Evaluate the integrand and accumulate element contributions
	  fe.detJxW *= 0.25*dA*wg[i]*wg[j];
	  if (!integrand.evalIntMx(elmInt,fe,time,X))
	    return false;
	}
      
      // Assembly of global system integral
      if (!glInt.assemble(elmInt,fe.iel))
	return false;
    }

  return true;
}


bool ASMs2Dmx::integrate (Integrand& integrand, int lIndex,
			  GlobalIntegral& glInt,
			  const TimeDomain& time,
			  const LintegralVec& locInt)
{
  if (!surf) return true; // silently ignore empty patches
  if (!basis1 || !basis2) return false;

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
  Matrix gpar[2];
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
  basis1->computeBasisGrid(gpar[0],gpar[1],spline1);
  basis2->computeBasisGrid(gpar[0],gpar[1],spline2);

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();

  MxFiniteElement fe(basis1->order_u()*basis1->order_v(),
		     basis2->order_u()*basis2->order_v());
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
      double dS = this->getParametricLength(iel,t1);
      if (dS < 0.0) return false; // topology error (probably logic error)

      // Set up control point coordinates for current element
      if (!this->getElementCoordinates(Xnod,iel)) return false;

      // Initialize element quantities
      IntVec::const_iterator f2start = MNPC[iel-1].begin() + fe.N1.size();
      if (!integrand.initElementBou(IntVec(MNPC[iel-1].begin(),f2start),
				    IntVec(f2start,MNPC[iel-1].end()),nb1))
	return false;

      // Caution: Unless locInt is empty, we assume it points to an array of
      // LocalIntegral pointers, of length at least the number of elements in
      // the model (as defined by the highest number in the MLGE array).
      // If the array is shorter than this, expect a segmentation fault.
      LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[fe.iel-1];


      // --- Integration loop over all Gauss points along the edge -------------

      int ip = (t1 == 1 ? i2-p2 : i1-p1)*nGauss;
      for (int i = 0; i < nGauss; i++, ip++)
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
	extractBasis(spline1[ip],fe.N1,dN1du);
	extractBasis(spline2[ip],fe.N2,dN2du);

	// Compute Jacobian inverse of the coordinate mapping and
	// basis function derivatives w.r.t. Cartesian coordinates
	if (geoUsesBasis1)
	{
	  fe.detJxW = utl::Jacobian(Jac,normal,fe.dN1dX,Xnod,dN1du,t1,t2);
	  fe.dN2dX.multiply(dN2du,Jac); // dN2dX = dN2du * J^-1
	}
	else
	{
	  fe.detJxW = utl::Jacobian(Jac,normal,fe.dN2dX,Xnod,dN2du,t1,t2);
	  fe.dN1dX.multiply(dN1du,Jac); // dN1dX = dN1du * J^-1
	}
	if (fe.detJxW == 0.0) continue; // skip singular points

	if (edgeDir < 0) normal *= -1.0;

	// Cartesian coordinates of current integration point
	X = Xnod * (geoUsesBasis1 ? fe.N1 : fe.N2);
	X.t = time.t;

	// Evaluate the integrand and accumulate element contributions
	fe.detJxW *= 0.5*dS*wg[i];
	if (!integrand.evalBouMx(elmInt,fe,time,X,normal))
	  return false;
      }

      // Assembly of global system integral
      if (!glInt.assemble(elmInt,fe.iel))
	return false;
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

  // Check if this point matches any of the basis1 control points (nodes)
  Vec3 Xnod;
  size_t inod = 1;
  RealArray::const_iterator cit = basis1->coefs_begin();
  for (int i = 0; cit != basis1->coefs_end(); cit++, i++)
  {
    if (i < nsd) Xnod[i] = *cit;
    if (i+1 == basis1->dimension())
      if (X.equal(Xnod,0.001))
        return inod;
      else
      {
        inod++;
        i = -1;
      }
  }

  return 0;
}


bool ASMs2Dmx::evalSolution (Matrix& sField, const Vector& locSol,
			     const RealArray* gpar, bool regular) const
{
  if (!basis1 || !basis2) return false;

  // Evaluate the basis functions at all points
  std::vector<Go::BasisPtsSf> spline1, spline2;
  if (regular)
  {
    basis1->computeBasisGrid(gpar[0],gpar[1],spline1);
    basis2->computeBasisGrid(gpar[0],gpar[1],spline2);
  }
  else if (gpar[0].size() == gpar[1].size())
  {
    spline1.resize(gpar[0].size());
    spline2.resize(gpar[0].size());
    for (size_t i = 0; i < spline1.size(); i++)
    {
      basis1->computeBasis(gpar[0][i],gpar[1][i],spline1[i]);
      basis2->computeBasis(gpar[0][i],gpar[1][i],spline2[i]);
    }
  }
  else
    return false;

  const int p1 = basis1->order_u();
  const int p2 = basis1->order_v();
  const int n1 = basis1->numCoefs_u();
  const int n2 = basis1->numCoefs_v();

  const int q1 = basis2->order_u();
  const int q2 = basis2->order_v();
  const int m1 = basis2->numCoefs_u();
  const int m2 = basis2->numCoefs_v();

  size_t nc1 = nf1;
  size_t nc2 = 0;
  if (nc1*nb1 < locSol.size())
    nc2 = (locSol.size() - nc1*nb1)/nb2;
  else
    nc1 = locSol.size()/nb1;

  if (nc1*nb1 + nc2*nb2 != locSol.size())
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

      utl::gather(ip,nc2,locSol,Xtmp,nc1*nb1);
      Xtmp.multiply(spline2[i].basisValues,Ztmp);

      Ytmp.insert(Ytmp.end(),Ztmp.begin(),Ztmp.end());
    }
    sField.fillColumn(1+i,Ytmp);
  }

  return true;
}


bool ASMs2Dmx::evalSolution (Matrix& sField, const Integrand& integrand,
			     const RealArray* gpar, bool regular) const
{
  sField.resize(0,0);

  if (!basis1 || !basis2) return false;

  // Evaluate the basis functions and their derivatives at all points
  std::vector<Go::BasisDerivsSf> spline1, spline2;
  if (regular)
  {
    basis1->computeBasisGrid(gpar[0],gpar[1],spline1);
    basis2->computeBasisGrid(gpar[0],gpar[1],spline2);
  }
  else if (gpar[0].size() == gpar[1].size())
  {
    std::vector<Go::BasisDerivsSf> tmpSpline(1);
    spline1.resize(gpar[0].size());
    spline2.resize(gpar[0].size());
    for (size_t i = 0; i < spline1.size(); i++)
    {
      basis1->computeBasisGrid(RealArray(1,gpar[0][i]),
			       RealArray(1,gpar[1][i]),
			       tmpSpline);
      spline1[i] = tmpSpline.front();
      basis2->computeBasisGrid(RealArray(1,gpar[0][i]),
			       RealArray(1,gpar[1][i]),
			       tmpSpline);
      spline2[i] = tmpSpline.front();
    }
    // TODO: Request a GoTools method replacing the above (see ASMs2D)
  }

  const int p1 = basis1->order_u();
  const int p2 = basis1->order_v();
  const int n1 = basis1->numCoefs_u();
  const int n2 = basis1->numCoefs_v();

  const int q1 = basis2->order_u();
  const int q2 = basis2->order_v();
  const int m1 = basis2->numCoefs_u();
  const int m2 = basis2->numCoefs_v();

  // Fetch nodal (control point) coordinates
  Matrix Xnod, Xtmp;
  this->getNodalCoordinates(Xnod);

  Vector N1(p1*p2), N2(q1*q2), solPt;
  Matrix dN1du, dN1dX, dN2du, dN2dX, Jac;
  Vec3   X;

  // Evaluate the secondary solution field at each point
  size_t nPoints = spline1.size();
  for (size_t i = 0; i < nPoints; i++)
  {
    // Fetch indices of the non-zero basis functions at this point
    IntVec ip1, ip2;
    scatterInd(n1,n2,p1,p2,spline1[i].left_idx,ip1);
    scatterInd(m1,m2,q1,q2,spline2[i].left_idx,ip2);

    // Fetch associated control point coordinates
    utl::gather(geoUsesBasis1 ? ip1 : ip2, nsd, Xnod, Xtmp);

    // Fetch basis function derivatives at current integration point
    extractBasis(spline1[i],N1,dN1du);
    extractBasis(spline2[i],N2,dN2du);

    // Compute Jacobian inverse of the coordinate mapping and
    // basis function derivatives w.r.t. Cartesian coordinates
    if (geoUsesBasis1)
      if (utl::Jacobian(Jac,dN1dX,Xtmp,dN1du) == 0.0) // Jac = (Xtmp * dN1du)^-1
	continue; // skip singular points
      else
	dN2dX.multiply(dN2du,Jac); // dN2dX = dN2du * J^-1
    else
      if (utl::Jacobian(Jac,dN2dX,Xtmp,dN2du) == 0.0) // Jac = (Xtmp * dN2du)^-1
	continue; // skip singular points
      else
	dN1dX.multiply(dN1du,Jac); // dN1dX = dN1du * J^-1

    // Cartesian coordinates of current integration point
    X = Xtmp * (geoUsesBasis1 ? N1 : N2);

    // Now evaluate the solution field
    if (!integrand.evalSol(solPt,N1,N2,dN1dX,dN2dX,X,ip1,ip2))
      return false;
    else if (sField.empty())
      sField.resize(solPt.size(),nPoints,true);

    sField.fillColumn(1+i,solPt);
  }

  return true;
}
