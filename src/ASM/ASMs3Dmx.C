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

#include "ASMs3Dmx.h"
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


ASMs3Dmx::ASMs3Dmx (const char* fName, bool checkRHS,
		    char n_f1, unsigned char n_f2)
  : ASMs3D(fName,checkRHS), ASMmxBase(n_f1 < 0 ? -n_f1 : n_f1,n_f2, n_f1 < 0)
{
  basis1 = basis2 = 0;
  nf = nf1 + nf2;
}


ASMs3Dmx::ASMs3Dmx (std::istream& is, bool checkRHS,
		    char n_f1, unsigned char n_f2)
  : ASMs3D(is,checkRHS), ASMmxBase(n_f1 < 0 ? -n_f1 : n_f1, n_f2, n_f1 < 0)
{
  basis1 = basis2 = 0;
  nf = nf1 + nf2;
}


void ASMs3Dmx::clear ()
{
  // Erase the spline data
  if (basis1) delete basis1;
  if (basis2) delete basis2;
  svol = basis1 = basis2 = 0;

  // Erase the FE data
  ASMs3D::clear();
}


unsigned char ASMs3Dmx::getNoFields (int basis) const
{
  switch (basis)
    {
    case 1: return nf1;
    case 2: return nf2;
    }

  return nf;
}


unsigned char ASMs3Dmx::getNodalDOFs (size_t inod) const
{
  return inod <= nb1 ? nf1 : nf2;
}


unsigned char ASMs3Dmx::getNodalBasis (size_t inod) const
{
  return inod <= nb1 ? 1 : 2;
}


void ASMs3Dmx::initMADOF (const int* sysMadof)
{
  this->initMx(MLGN,sysMadof);
}


void ASMs3Dmx::extractNodeVec (const Vector& globRes, Vector& nodeVec,
			       unsigned char) const
{
  this->extractNodeVecMx(globRes,nodeVec);
}


bool ASMs3Dmx::getSolution (Matrix& sField, const Vector& locSol,
			    const IntVec& nodes) const
{
  return this->getSolutionMx(sField,locSol,nodes);
}


bool ASMs3Dmx::generateFEMTopology ()
{
  if (!svol) return false;

  if (!basis1 && !basis2)
  {
    // With mixed methods we need two separates spline spaces
    basis1 = new Go::SplineVolume(*svol);
    basis2 = svol;

    // Order-elevate basis1 such that it is of one degree higher than basis2
    basis1->raiseOrder(1,1,1);

    // Define which basis that should be used to represent the geometry
    if (geoUsesBasis1) svol = basis1;
  }

  const int n1 = basis1->numCoefs(0);
  const int n2 = basis1->numCoefs(1);
  const int n3 = basis1->numCoefs(2);
  const int m1 = basis2->numCoefs(0);
  const int m2 = basis2->numCoefs(1);
  const int m3 = basis2->numCoefs(2);
  if (!nodeInd.empty())
  {
    if (nodeInd.size() == nb1 + nb2) return true;
    std::cerr <<" *** ASMs3Dmx::generateFEMTopology: Inconsistency between the"
	      <<" number of FE nodes "<< nodeInd.size()
	      <<"\n     and the number of spline coefficients "<< nb1 + nb2
	      <<" in the patch."<< std::endl;
    return false;
  }

  nb1 = n1*n2*n3; // Number of functions in first basis
  nb2 = m1*m2*m3; // Number of functions in second basis

  const int p1 = basis1->order(0);
  const int p2 = basis1->order(1);
  const int p3 = basis1->order(2);
  const int q1 = basis2->order(0);
  const int q2 = basis2->order(1);
  const int q3 = basis2->order(2);

#ifdef SP_DEBUG
  std::cout <<"numCoefs: "<< n1 <<" "<< n2 <<" "<< n3
	    <<", "<< m1 <<" "<< m2 <<" "<< m3;
  std::cout <<"\norder: "<< p1 <<" "<< p2 <<" "<< p3
	    <<", "<< q1 <<" "<< q2 <<" "<< q3;
  for (int d = 0; d < 3; d++)
  {
    std::cout <<"\nd"<< char('u'+d) <<':';
    for (int i = 0; i < basis1->numCoefs(d); i++)
      std::cout <<' '<< basis1->knotSpan(d,i);
    for (int j = 0; j < basis2->numCoefs(d); j++)
      std::cout <<' '<< basis2->knotSpan(d,j);
  }
  std::cout << std::endl;
#endif
  // Consistency checks, just to be fool-proof
  if (m1 <  2 || m2 <  2 || m3 <  2) return false;
  if (q1 <  1 || q2 <  1 || q3 <  1) return false;
  if (p1 > n1 || p2 > n2 || p3 > n3) return false;
  if (q1 > m1 || q2 > m2 || q3 > m3) return false;

  MLGE.resize((n1-p1+1)*(n2-p2+1)*(n3-p3+1),0);
  MLGN.resize(nb1 + nb2);
  MNPC.resize(MLGE.size());
  nodeInd.resize(MLGN.size());

  int i1, i2, i3, j1, j2, j3;
  size_t iel, inod = 0;
  for (i3 = 0; i3 < n3; i3++)
    for (i2 = 0; i2 < n2; i2++)
      for (i1 = 0; i1 < n1; i1++)
      {
	nodeInd[inod].I = i1;
	nodeInd[inod].J = i2;
	nodeInd[inod].K = i3;
	MLGN[inod++]    = ++gNod;
      }

  for (i3 = 0; i3 < m3; i3++)
    for (i2 = 0; i2 < m2; i2++)
      for (i1 = 0; i1 < m1; i1++)
      {
	nodeInd[inod].I = i1;
	nodeInd[inod].J = i2;
	nodeInd[inod].K = i3;
	MLGN[inod++]    = ++gNod;
      }

  if (geoUsesBasis1)
  {
    // Create nodal connectivities for basis 1
    iel = inod = 0;
    for (i3 = 1; i3 <= n3; i3++)
      for (i2 = 1; i2 <= n2; i2++)
	for (i1 = 1; i1 <= n1; i1++, inod++)
	  if (i1 >= p1 && i2 >= p2 && i3 >= p3)
	  {
	    if (basis1->knotSpan(0,i1-1) > 0.0)
	      if (basis1->knotSpan(1,i2-1) > 0.0)
		if (basis1->knotSpan(2,i3-1) > 0.0)
		{
		  MLGE[iel] = ++gEl; // global element number over all patches
		  MNPC[iel].resize(p1*p2*p3+q1*q2*q3,0);

		  int lnod = 0;
		  for (j3 = p3-1; j3 >= 0; j3--)
		    for (j2 = p2-1; j2 >= 0; j2--)
		      for (j1 = p1-1; j1 >= 0; j1--)
			MNPC[iel][lnod++] = inod - n1*n2*j3 - n1*j2 - j1;
		}

	    iel++;
	  }

    // Create nodal connectivities for basis 2
    iel = 0;
    for (i3 = 1; i3 <= m3; i3++)
      for (i2 = 1; i2 <= m2; i2++)
	for (i1 = 1; i1 <= m1; i1++, inod++)
	  if (i1 >= q1 && i2 >= q2 && i3 >= q3)
	    if (basis2->knotSpan(0,i1-1) > 0.0)
	      if (basis2->knotSpan(1,i2-1) > 0.0)
		if (basis2->knotSpan(2,i3-1) > 0.0)
		{
		  while (iel < MNPC.size() && MNPC[iel].empty()) iel++;

		  int lnod = p1*p2*p3;
		  for (j3 = q3-1; j3 >= 0; j3--)
		    for (j2 = q2-1; j2 >= 0; j2--)
		      for (j1 = q1-1; j1 >= 0; j1--)
			MNPC[iel][lnod++] = inod - m1*m2*j3 - m1*j2 - j1;

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
	    if (basis2->knotSpan(0,i1-1) > 0.0)
	      if (basis2->knotSpan(1,i2-1) > 0.0)
		if (basis2->knotSpan(2,i3-1) > 0.0)
		{
		  MLGE[iel] = ++gEl; // global element number over all patches
		  MNPC[iel].resize(p1*p2*p3+q1*q2*q3,0);

		  int lnod = p1*p2*p3;
		  for (j3 = q3-1; j3 >= 0; j3--)
		    for (j2 = q2-1; j2 >= 0; j2--)
		      for (j1 = q1-1; j1 >= 0; j1--)
			MNPC[iel][lnod++] = inod - m1*m2*j3 - m1*j2 - j1;
		}

	    iel++;
	  }

    // Create nodal connectivities for basis 1
    iel = inod = 0;
    for (i3 = 1; i3 <= n3; i3++)
      for (i2 = 1; i2 <= n2; i2++)
	for (i1 = 1; i1 <= n1; i1++, inod++)
	  if (i1 >= p1 && i2 >= p2 && i3 >= p3)
	    if (basis1->knotSpan(0,i1-1) > 0.0)
	      if (basis1->knotSpan(1,i2-1) > 0.0)
		if (basis1->knotSpan(2,i3-1) > 0.0)
		{
		  while (iel < MNPC.size() && MNPC[iel].empty()) iel++;

		  int lnod = 0;
		  for (j3 = p3-1; j3 >= 0; j3--)
		    for (j2 = p2-1; j2 >= 0; j2--)
		      for (j1 = p1-1; j1 >= 0; j1--)
			MNPC[iel][lnod++] = inod - n1*n2*j3 - n1*j2 - j1;

		  iel++;
		}
  }

#ifdef SP_DEBUG
  std::cout <<"NEL = "<< MLGE.size() <<" NNOD = "<< MLGN.size() << std::endl;
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

  return this->connectBasis(face,neighbor,nface,norient,1,0,0)
    &&   this->connectBasis(face,neighbor,nface,norient,2,nb1,neighMx->nb1);
}


void ASMs3Dmx::closeFaces (int dir, int, int)
{
  this->ASMs3D::closeFaces(dir,1,1);
  this->ASMs3D::closeFaces(dir,2,nb1+1);
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
  size_t lnod0 = basis1->order(0)*basis1->order(1)*basis1->order(2);
  if (geoUsesBasis1) lnod0 = 0;

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
  if (inod <= nb1)
    cit = basis1->coefs_begin()
      + ((K*basis1->numCoefs(1)+J)*basis1->numCoefs(0)+I) * basis1->dimension();
  else
    cit = basis2->coefs_begin()
      + ((K*basis2->numCoefs(1)+J)*basis2->numCoefs(0)+I) * basis2->dimension();

  return Vec3(*cit,*(cit+1),*(cit+2));
}


bool ASMs3Dmx::getSize (int& n1, int& n2, int& n3, int basis) const
{
  switch (basis)
    {
    case 1:
      if (!basis1) return false;
      n1 = basis1->numCoefs(0);
      n2 = basis1->numCoefs(1);
      n3 = basis1->numCoefs(2);
      return true;

    case 2:
      if (!basis2) return false;
      n1 = basis2->numCoefs(0);
      n2 = basis2->numCoefs(1);
      n3 = basis2->numCoefs(2);
      return true;
    }

  return this->ASMs3D::getSize(n1,n2,n3);
}


bool ASMs3Dmx::integrate (Integrand& integrand,
			  GlobalIntegral& glInt,
			  const TimeDomain& time,
			  const LintegralVec& locInt)
{
  if (!svol) return true; // silently ignore empty patches
  if (!basis1 || !basis2) return false;

  PROFILE2("ASMs3Dmx::integrate(I)");

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  int dir;
  Matrix gpar[3];
  for (dir = 0; dir < 3; dir++)
  {
    int pm1 = svol->order(dir) - 1;
    RealArray::const_iterator uit = svol->basis(dir).begin() + pm1;
    double ucurr, uprev = *(uit++);
    int nCol = svol->numCoefs(dir) - pm1;
    gpar[dir].resize(nGauss,nCol);
    for (int j = 1; j <= nCol; uit++, j++)
    {
      ucurr = *uit;
      for (int i = 1; i <= nGauss; i++)
	gpar[dir](i,j) = 0.5*((ucurr-uprev)*xg[i-1] + ucurr+uprev);
      uprev = ucurr;
    }
  }

  // Evaluate basis function derivatives at all integration points
  std::vector<Go::BasisDerivs> spline1, spline2;
  basis1->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline1);
  basis2->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline2);

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);

  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;

  MxFiniteElement fe(basis1->order(0)*basis1->order(1)*basis1->order(2),
		     basis2->order(0)*basis2->order(1)*basis2->order(2));
  Matrix dN1du, dN2du, Xnod, Jac;
  Vec4   X;


  // === Assembly loop over all elements in the patch ==========================

  int iel = 1;
  for (int i3 = p3; i3 <= n3; i3++)
    for (int i2 = p2; i2 <= n2; i2++)
      for (int i1 = p1; i1 <= n1; i1++, iel++)
      {
	if (MLGE[iel-1] < 1) continue; // zero-volume element

	// Get element volume in the parameter space
	double dV = this->getParametricVolume(iel);
	if (dV < 0.0) return false; // topology error (probably logic error)

	// Set up control point (nodal) coordinates for current element
	if (!this->getElementCoordinates(Xnod,iel)) return false;

	// Initialize element quantities
	IntVec::iterator f2start = MNPC[iel-1].begin() + fe.N1.size();
	if (!integrand.initElement(IntVec(MNPC[iel-1].begin(),f2start),
				   IntVec(f2start,MNPC[iel-1].end()),nb1))
	  return false;

	// Caution: Unless locInt is empty, we assume it points to an array of
	// LocalIntegral pointers, of length at least the number of elements in
	// the model (as defined by the highest number in the MLGE array).
	// If the array is shorter than this, expect a segmentation fault.
	LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[MLGE[iel-1]-1];


	// --- Integration loop over all Gauss points in each direction --------

	int ip = (((i3-p3)*nGauss*nel2 + i2-p2)*nGauss*nel1 + i1-p1)*nGauss;
	for (int k = 0; k < nGauss; k++, ip += nGauss*(nel2-1)*nGauss*nel1)
	  for (int j = 0; j < nGauss; j++, ip += nGauss*(nel1-1))
	    for (int i = 0; i < nGauss; i++, ip++)
	    {
	      // Parameter values of current integration point
	      fe.u = gpar[0](i+1,i1-p1+1);
	      fe.v = gpar[1](j+1,i2-p2+1);
	      fe.w = gpar[2](k+1,i3-p3+1);

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
	      fe.detJxW *= 0.125*dV*wg[i]*wg[j]*wg[k];
	      if (!integrand.evalIntMx(elmInt,fe,time,X))
		return false;
	    }

	// Assembly of global system integral
	if (!glInt.assemble(elmInt,MLGE[iel-1]))
	  return false;
      }

  return true;
}


bool ASMs3Dmx::integrate (Integrand& integrand, int lIndex,
			  GlobalIntegral& glInt,
			  const TimeDomain& time,
			  const LintegralVec& locInt)
{
  if (!svol) return true; // silently ignore empty patches
  if (!basis1 || !basis2) return false;

  PROFILE2("ASMs3Dmx::integrate(B)");

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Find the parametric direction of the face normal {-3,-2,-1, 1, 2, 3}
  const int faceDir = (lIndex+1)/(lIndex%2 ? -2 : 2);

  const int t1 = 1 + abs(faceDir)%3; // first tangent direction
  const int t2 = 1 + t1%3;           // second tangent direction

  // Compute parameter values of the Gauss points over the whole patch face
  Matrix gpar[3];
  for (int d = 0; d < 3; d++)
    if (-1-d == faceDir)
    {
      gpar[d].resize(1,1);
      gpar[d](1,1) = svol->startparam(d);
    }
    else if (1+d == faceDir)
    {
      gpar[d].resize(1,1);
      gpar[d](1,1) = svol->endparam(d);
    }
    else
    {
      int pm1 = svol->order(d) - 1;
      RealArray::const_iterator uit = svol->basis(d).begin() + pm1;
      double ucurr, uprev = *(uit++);
      int nCol = svol->numCoefs(d) - pm1;
      gpar[d].resize(nGauss,nCol);
      for (int j = 1; j <= nCol; uit++, j++)
      {
	ucurr = *uit;
	for (int i = 1; i <= nGauss; i++)
	  gpar[d](i,j) = 0.5*((ucurr-uprev)*xg[i-1] + ucurr+uprev);
	uprev = ucurr;
      }
    }

  // Evaluate basis function derivatives at all integration points
  std::vector<Go::BasisDerivs> spline1, spline2;
  basis1->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline1);
  basis2->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline2);

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;

  MxFiniteElement fe(basis1->order(0)*basis1->order(1)*basis1->order(2),
		     basis2->order(0)*basis2->order(1)*basis2->order(2));
  fe.u = gpar[0](1,1);
  fe.v = gpar[1](1,1);
  fe.w = gpar[2](1,1);

  Matrix dN1du, dN2du, Xnod, Jac;
  Vec4   X;
  Vec3   normal;


  // === Assembly loop over all elements on the patch face =====================

  int iel = 1;
  for (int i3 = p3; i3 <= n3; i3++)
    for (int i2 = p2; i2 <= n2; i2++)
      for (int i1 = p1; i1 <= n1; i1++, iel++)
      {
	if (MLGE[iel-1] < 1) continue; // zero-volume element

	// Skip elements that are not on current boundary face
	bool skipMe = false;
	switch (faceDir)
	  {
	  case -1: if (i1 > p1) skipMe = true; break;
	  case  1: if (i1 < n1) skipMe = true; break;
	  case -2: if (i2 > p2) skipMe = true; break;
	  case  2: if (i2 < n2) skipMe = true; break;
	  case -3: if (i3 > p3) skipMe = true; break;
	  case  3: if (i3 < n3) skipMe = true; break;
	  }
	if (skipMe) continue;

	// Get element face area in the parameter space
	double dA = this->getParametricArea(iel,abs(faceDir));
	if (dA < 0.0) return false; // topology error (probably logic error)

	// Set up control point coordinates for current element
	if (!this->getElementCoordinates(Xnod,iel)) return false;

	// Initialize element quantities
	IntVec::iterator f2start = MNPC[iel-1].begin() + fe.N1.size();
	if (!integrand.initElementBou(IntVec(MNPC[iel-1].begin(),f2start),
				      IntVec(f2start,MNPC[iel-1].end()),nb1))
	  return false;

	// Define some loop control variables depending on which face we are on
	int nf1, j1, j2;
	switch (abs(faceDir))
	  {
	  case 1: nf1 = nel2; j2 = i3-p3; j1 = i2-p2; break;
	  case 2: nf1 = nel1; j2 = i3-p3; j1 = i1-p1; break;
	  case 3: nf1 = nel1; j2 = i2-p2; j1 = i1-p1; break;
	  default: nf1 = j1 = j2 = 0;
	  }

	// Caution: Unless locInt is empty, we assume it points to an array of
	// LocalIntegral pointers, of length at least the number of elements in
	// the model (as defined by the highest number in the MLGE array).
	// If the array is shorter than this, expect a segmentation fault.
	LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[MLGE[iel-1]-1];


	// --- Integration loop over all Gauss points in each direction --------

	int k1, k2, k3;
	int ip = (j2*nGauss*nf1 + j1)*nGauss;
	for (int j = 0; j < nGauss; j++, ip += nGauss*(nf1-1))
	  for (int i = 0; i < nGauss; i++, ip++)
	  {
	    // Parameter values of current integration point
	    switch (abs(faceDir)) {
	    case 1: k2 = i+1; k3 = j+1; k1 = 0; break;
	    case 2: k1 = i+1; k3 = j+1; k2 = 0; break;
	    case 3: k1 = i+1; k2 = j+1; k3 = 0; break;
	    default: k1 = k2 = k3 = 0;
	    }
	    if (gpar[0].size() > 1) fe.u = gpar[0](k1,i1-p1+1);
	    if (gpar[1].size() > 1) fe.v = gpar[1](k2,i2-p2+1);
	    if (gpar[2].size() > 1) fe.w = gpar[2](k3,i3-p3+1);

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

	    if (faceDir < 0) normal *= -1.0;

	    // Cartesian coordinates of current integration point
	    X = Xnod * (geoUsesBasis1 ? fe.N1 : fe.N2);
	    X.t = time.t;

	    // Evaluate the integrand and accumulate element contributions
	    fe.detJxW *= 0.25*dA*wg[i]*wg[j];
	    if (!integrand.evalBouMx(elmInt,fe,time,X,normal))
	      return false;
	  }

	// Assembly of global system integral
	if (!glInt.assemble(elmInt,MLGE[iel-1]))
	  return false;
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

  // Check if this point matches any of the basis1 control points (nodes)
  Vec3 Xnod;
  size_t inod = 1;
  RealArray::const_iterator cit = basis1->coefs_begin();
  for (i = 0; cit != basis1->coefs_end(); cit++, i++)
  {
    if (i < 3) Xnod[i] = *cit;
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


bool ASMs3Dmx::evalSolution (Matrix& sField, const Vector& locSol,
			     const RealArray* gpar, bool regular) const
{
  if (!basis1 || !basis2) return false;

  // Evaluate the basis functions at all points
  std::vector<Go::BasisPts> spline1, spline2;
  if (regular)
  {
    basis1->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline1);
    basis2->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline2);
  }
  else if (gpar[0].size() == gpar[1].size() && gpar[0].size() == gpar[2].size())
  {
    std::vector<Go::BasisPts> tmpSpline(1);
    spline1.resize(gpar[0].size());
    spline2.resize(gpar[0].size());
    for (size_t i = 0; i < spline1.size(); i++)
    {
      basis1->computeBasisGrid(RealArray(1,gpar[0][i]),
			       RealArray(1,gpar[1][i]),
			       RealArray(1,gpar[2][i]),
			       tmpSpline);
      spline1[i] = tmpSpline.front();
      basis2->computeBasisGrid(RealArray(1,gpar[0][i]),
			       RealArray(1,gpar[1][i]),
			       RealArray(1,gpar[2][i]),
			       tmpSpline);
      spline2[i] = tmpSpline.front();
    }
    // TODO: Request a GoTools method replacing the above (see ASMs3D)
  }

  const int p1 = basis1->order(0);
  const int p2 = basis1->order(1);
  const int p3 = basis1->order(2);
  const int n1 = basis1->numCoefs(0);
  const int n2 = basis1->numCoefs(1);
  const int n3 = basis1->numCoefs(2);

  const int q1 = basis2->order(0);
  const int q2 = basis2->order(1);
  const int q3 = basis2->order(2);
  const int m1 = basis2->numCoefs(0);
  const int m2 = basis2->numCoefs(1);
  const int m3 = basis2->numCoefs(2);

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
    scatterInd(n1,n2,n3,p1,p2,p3,spline1[i].left_idx,ip);

    utl::gather(ip,nc1,locSol,Xtmp);
    Xtmp.multiply(spline1[i].basisValues,Ytmp);

    if (nc2 > 0)
    {
      ip.clear();
      scatterInd(m1,m2,m3,q1,q2,q3,spline2[i].left_idx,ip);

      utl::gather(ip,nc2,locSol,Xtmp,nc1*nb1);
      Xtmp.multiply(spline2[i].basisValues,Ztmp);

      Ytmp.insert(Ytmp.end(),Ztmp.begin(),Ztmp.end());
    }
    sField.fillColumn(1+i,Ytmp);
  }

  return true;
}


bool ASMs3Dmx::evalSolution (Matrix& sField, const Integrand& integrand,
			     const RealArray* gpar, bool regular) const
{
  sField.resize(0,0);

  if (!basis1 || !basis2) return false;

  // Evaluate the basis functions and their derivatives at all points
  std::vector<Go::BasisDerivs> spline1, spline2;
  if (regular)
  {
    basis1->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline1);
    basis2->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline2);
  }
  else if (gpar[0].size() == gpar[1].size() && gpar[0].size() == gpar[2].size())
  {
    std::vector<Go::BasisDerivs> tmpSpline(1);
    spline1.resize(gpar[0].size());
    spline2.resize(gpar[0].size());
    for (size_t i = 0; i < spline1.size(); i++)
    {
      basis1->computeBasisGrid(RealArray(1,gpar[0][i]),
			       RealArray(1,gpar[1][i]),
			       RealArray(1,gpar[2][i]),
			       tmpSpline);
      spline1[i] = tmpSpline.front();
      basis2->computeBasisGrid(RealArray(1,gpar[0][i]),
			       RealArray(1,gpar[1][i]),
			       RealArray(1,gpar[2][i]),
			       tmpSpline);
      spline2[i] = tmpSpline.front();
    }
    // TODO: Request a GoTools method replacing the above (see ASMs3D)
  }

  const int p1 = basis1->order(0);
  const int p2 = basis1->order(1);
  const int p3 = basis1->order(2);
  const int n1 = basis1->numCoefs(0);
  const int n2 = basis1->numCoefs(1);
  const int n3 = basis1->numCoefs(2);

  const int q1 = basis2->order(0);
  const int q2 = basis2->order(1);
  const int q3 = basis2->order(2);
  const int m1 = basis2->numCoefs(0);
  const int m2 = basis2->numCoefs(1);
  const int m3 = basis2->numCoefs(2);

  // Fetch nodal (control point) coordinates
  Matrix Xnod, Xtmp;
  this->getNodalCoordinates(Xnod);

  Vector N1(p1*p2*p3), N2(q1*q2*q3), solPt;
  Matrix dN1du, dN1dX, dN2du, dN2dX, Jac;
  Vec3   X;

  // Evaluate the secondary solution field at each point
  size_t nPoints = spline1.size();
  for (size_t i = 0; i < nPoints; i++)
  {
    // Fetch indices of the non-zero basis functions at this point
    IntVec ip1, ip2;
    scatterInd(n1,n2,n3,p1,p2,p3,spline1[i].left_idx,ip1);
    scatterInd(m1,m2,m3,q1,q2,q3,spline2[i].left_idx,ip2);

    // Fetch associated control point coordinates
    utl::gather(geoUsesBasis1 ? ip1 : ip2, 3, Xnod, Xtmp);

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
