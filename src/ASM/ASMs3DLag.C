// $Id$
//==============================================================================
//!
//! \file ASMs3DLag.C
//!
//! \date Feb 10 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for assembly of structured 3D Lagrange FE models.
//!
//==============================================================================

#include "GoTools/trivariate/SplineVolume.h"

#include "ASMs3DLag.h"
#include "Lagrange.h"
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "ElementBlock.h"
#include "Utilities.h"
#include "Vec3Oper.h"


ASMs3DLag::ASMs3DLag (const char* fName, bool checkRHS, unsigned char n_f)
  : ASMs3D(fName,checkRHS,n_f)
{
  nx = ny = nz = 0;
}


ASMs3DLag::ASMs3DLag (std::istream& is, bool checkRHS, unsigned char n_f)
  : ASMs3D(is,checkRHS,n_f)
{
  nx = ny = nz = 0;
}


void ASMs3DLag::clear (bool retainGeometry)
{
  coord.clear();
  nx = ny = nz = 0;

  this->ASMs3D::clear(retainGeometry);
}


bool ASMs3DLag::generateFEMTopology ()
{
  if (!svol) return false;

  // Order of basis in the three parametric directions (order = degree + 1)
  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  // Evaluate the parametric values
  RealArray gpar[3];
  if (!this->getGridParameters(gpar[0],0,p1-1)) return false;
  if (!this->getGridParameters(gpar[1],1,p2-1)) return false;
  if (!this->getGridParameters(gpar[2],2,p3-1)) return false;

  // Number of nodes in each direction
  nx = gpar[0].size();
  ny = gpar[1].size();
  nz = gpar[2].size();

  if (!coord.empty())
    return coord.size() == nx*ny*nz;

  // Number of elements in each direction
  const int nelx = (nx-1)/(p1-1);
  const int nely = (ny-1)/(p2-1);
  const int nelz = (nz-1)/(p3-1);

  // Evaluate the nodal coordinates in the physical space
  RealArray XYZ(svol->dimension()*nx*ny*nz);
  svol->gridEvaluator(gpar[0],gpar[1],gpar[2],XYZ);

  size_t i1, j1;
  MLGN.resize(nx*ny*nz);
  coord.resize(nx*ny*nz);
  for (i1 = j1 = 0; i1 < coord.size(); i1++)
  {
    MLGN[i1] = ++gNod;
    coord[i1][0] = XYZ[j1++];
    coord[i1][1] = XYZ[j1++];
    coord[i1][2] = XYZ[j1++];
  }

  // Number of elements in patch
  const int nel = nelx*nely*nelz;
  // Number of nodes per element
  const int nen = p1*p2*p3;
  // Number of nodes in a xy-surface of an element
  const int ct  = p1*p2;

  // Connectivity array: local --> global node relation
  MLGE.resize(nel);
  MNPC.resize(nel);

  int i, j, k, a, b, c, iel = 0;
  for (k = 0; k < nelz; k++)
    for (j = 0; j < nely; j++)
      for (i = 0; i < nelx; i++, iel++)
      {
	MLGE[iel] = ++gEl;
	MNPC[iel].resize(nen);
	// First node in current element
	int corner = (p3-1)*(nx*ny)*k + (p2-1)*nx*j + (p1-1)*i;

	for (c = 0; c < p3; c++)
	{
	  int cornod = ct*c;
	  MNPC[iel][cornod] = corner + c*nx*ny;
	  for (b = 1; b < p2; b++)
	  {
	    int facenod = cornod + b*p1;
	    MNPC[iel][facenod] = MNPC[iel][cornod] + b*nx;
	    for (a = 1; a < p1; a++)
	    {
	      MNPC[iel][facenod+a] = MNPC[iel][facenod] + a;
	      MNPC[iel][cornod+a] = MNPC[iel][cornod] + a;
	    }
	  }
	}
      }

  return true;
}


bool ASMs3DLag::getElementCoordinates (Matrix& X, int iel) const
{
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs3DLag::getElementCoordinates: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return false;
  }

  // Number of nodes per element
  const size_t nen = svol->order(0)*svol->order(1)*svol->order(2);

  X.resize(3,nen);
  for (size_t i = 0; i < nen; i++)
    X.fillColumn(i+1,coord[MNPC[iel-1][i]].ptr());

  return true;
}


void ASMs3DLag::getNodalCoordinates (Matrix& X) const
{
  X.resize(3,coord.size());

  for (size_t inod = 0; inod < coord.size(); inod++)
    X.fillColumn(inod+1,coord[inod].ptr());
}


Vec3 ASMs3DLag::getCoord (size_t inod) const
{
  if (inod < 1 || inod > coord.size()) return Vec3();

  return coord[inod-1];
}


bool ASMs3DLag::getSize (int& n1, int& n2, int& n3, int) const
{
  n1 = nx;
  n2 = ny;
  n3 = nz;

  return true;
}


bool ASMs3DLag::integrate (Integrand& integrand,
			   GlobalIntegral& glInt,
			   const TimeDomain& time,
			   const LintegralVec& locInt)
{
  if (!svol) return true; // silently ignore empty patches

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Points for selective reduced integration
  const double* xr = 0;
  int nRed = integrand.getIntegrandType() - 10;
  if (nRed < 1)
    nRed = nRed < 0 ? nGauss : 0;
  else if (!(xr = GaussQuadrature::getCoord(nRed)))
    return false;

  // Get parametric coordinates of the elements
  RealArray upar, vpar, wpar;
  this->getGridParameters(upar,0,1);
  this->getGridParameters(vpar,1,1);
  this->getGridParameters(wpar,2,1);

  // Number of elements in each direction
  const int nelx = upar.size() - 1;
  const int nely = vpar.size() - 1;
  const int nelz = wpar.size() - 1;

  // Order of basis in the three parametric directions (order = degree + 1)
  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  FiniteElement fe(p1*p2*p3);
  Matrix dNdu, Xnod, Jac;
  Vec4   X;


  // === Assembly loop over all elements in the patch ==========================

  int iel = 1;
  for (int i3 = 0; i3 < nelz; i3++)
    for (int i2 = 0; i2 < nely; i2++)
      for (int i1 = 0; i1 < nelx; i1++, iel++)
      {
	fe.iel = MLGE[iel-1];

	// Set up nodal point coordinates for current element
	if (!this->getElementCoordinates(Xnod,iel)) return false;

	if (integrand.getIntegrandType() == 4)
	{
	  // Compute the element "center" (average of element node coordinates)
	  X = 0.0;
	  for (size_t i = 1; i <= 3; i++)
	    for (size_t j = 1; j <= Xnod.cols(); j++)
	      X[i-1] += Xnod(i,j);

	  X *= 1.0/(double)Xnod.cols();
	}

	// Initialize element quantities
	if (!integrand.initElement(MNPC[iel-1],X,nRed*nRed*nRed))
	  return false;

	// Caution: Unless locInt is empty, we assume it points to an array of
	// LocalIntegral pointers, of length at least the number of elements in
	// the model (as defined by the highest number in the MLGE array).
	// If the array is shorter than this, expect a segmentation fault.
	LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[fe.iel-1];


	// --- Selective reduced integration loop ------------------------------

	if (integrand.getIntegrandType() > 10)
	  for (int k = 0; k < nRed; k++)
	    for (int j = 0; j < nRed; j++)
	      for (int i = 0; i < nRed; i++)
	      {
		// Local element coordinates of current integration point
		fe.xi   = xr[i];
		fe.eta  = xr[j];
		fe.zeta = xr[k];

		// Parameter value of current integration point
		fe.u = 0.5*(upar[i1]*(1.0-xr[i]) + upar[i1+1]*(1.0+xr[i]));
		fe.v = 0.5*(vpar[i2]*(1.0-xr[j]) + vpar[i2+1]*(1.0+xr[j]));
		fe.w = 0.5*(wpar[i3]*(1.0-xr[k]) + wpar[i3+1]*(1.0+xr[k]));

		// Compute basis function derivatives at current point
		// using tensor product of one-dimensional Lagrange polynomials
		if (!Lagrange::computeBasis(fe.N,dNdu,
					    p1,xr[i],p2,xr[j],p3,xr[k]))
		  return false;

		// Compute Jacobian inverse and derivatives
		fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

		// Compute the reduced integration terms of the integrand
		if (!integrand.reducedInt(fe))
		  return false;
	      }


	// --- Integration loop over all Gauss points in each direction --------

	for (int k = 0; k < nGauss; k++)
	  for (int j = 0; j < nGauss; j++)
	    for (int i = 0; i < nGauss; i++)
	    {
	      // Local element coordinates of current integration point
	      fe.xi   = xg[i];
	      fe.eta  = xg[j];
	      fe.zeta = xg[k];

	      // Parameter value of current integration point
	      fe.u = 0.5*(upar[i1]*(1.0-xg[i]) + upar[i1+1]*(1.0+xg[i]));
	      fe.v = 0.5*(vpar[i2]*(1.0-xg[j]) + vpar[i2+1]*(1.0+xg[j]));
	      fe.w = 0.5*(wpar[i3]*(1.0-xg[k]) + wpar[i3+1]*(1.0+xg[k]));

	      // Compute basis function derivatives at current integration point
	      // using tensor product of one-dimensional Lagrange polynomials
	      if (!Lagrange::computeBasis(fe.N,dNdu,p1,xg[i],p2,xg[j],p3,xg[k]))
		return false;

	      // Compute Jacobian inverse of coordinate mapping and derivatives
	      fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);
	      if (fe.detJxW == 0.0) continue; // skip singular points

	      // Cartesian coordinates of current integration point
	      X = Xnod * fe.N;
	      X.t = time.t;

	      // Evaluate the integrand and accumulate element contributions
	      fe.detJxW *= wg[i]*wg[j]*wg[k];
	      if (!integrand.evalInt(elmInt,fe,time,X))
		return false;
	    }

	// Finalize the element quantities
 	if (!integrand.finalizeElement(elmInt,time))
 	  return false;

	// Assembly of global system integral
	if (!glInt.assemble(elmInt,fe.iel))
	  return false;
      }

  return true;
}


bool ASMs3DLag::integrate (Integrand& integrand, int lIndex,
			   GlobalIntegral& glInt,
			   const TimeDomain& time,
			   const LintegralVec& locInt)
{
  if (!svol) return true; // silently ignore empty patches

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Find the parametric direction of the face normal {-3,-2,-1, 1, 2, 3}
  const int faceDir = (lIndex+1)/(lIndex%2 ? -2 : 2);

  const int t0 = abs(faceDir); // unsigned normal direction of the face
  const int t1 = 1 + t0%3; // first tangent direction of the face
  const int t2 = 1 + t1%3; // second tangent direction of the face

  // Order of basis in the three parametric directions (order = degree + 1)
  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  // Number of elements in each direction
  const int nelx = (nx-1)/(p1-1);
  const int nely = (ny-1)/(p2-1);
  const int nelz = (nz-1)/(p3-1);

  // Get parametric coordinates of the elements
  RealArray upar, vpar, wpar;
  if (t0 == 1)
    upar.resize(1,faceDir < 0 ? svol->startparam(0) : svol->endparam(0));
  else if (t0 == 2)
    vpar.resize(1,faceDir < 0 ? svol->startparam(1) : svol->endparam(1));
  else if (t0 == 3)
    wpar.resize(1,faceDir < 0 ? svol->startparam(2) : svol->endparam(2));

  if (upar.empty()) this->getGridParameters(upar,0,1);
  if (vpar.empty()) this->getGridParameters(vpar,1,1);
  if (wpar.empty()) this->getGridParameters(wpar,2,1);

  FiniteElement fe(p1*p2*p3);
  fe.u = upar.front();
  fe.v = vpar.front();
  fe.w = wpar.front();

  Matrix dNdu, Xnod, Jac;
  Vec4   X;
  Vec3   normal;
  double xi[3];


  // === Assembly loop over all elements on the patch face =====================

  int iel = 1;
  for (int i3 = 0; i3 < nelz; i3++)
    for (int i2 = 0; i2 < nely; i2++)
      for (int i1 = 0; i1 < nelx; i1++, iel++)
      {
	fe.iel = MLGE[iel-1];

	// Skip elements that are not on current boundary face
	bool skipMe = false;
	switch (faceDir)
	  {
	  case -1: if (i1 > 0)      skipMe = true; break;
	  case  1: if (i1 < nelx-1) skipMe = true; break;
	  case -2: if (i2 > 0)      skipMe = true; break;
	  case  2: if (i2 < nely-1) skipMe = true; break;
	  case -3: if (i3 > 0)      skipMe = true; break;
	  case  3: if (i3 < nelz-1) skipMe = true; break;
	  }
	if (skipMe) continue;

	// Set up nodal point coordinates for current element
	if (!this->getElementCoordinates(Xnod,iel)) return false;

	// Initialize element quantities
	if (!integrand.initElementBou(MNPC[iel-1])) return false;

	// Caution: Unless locInt is empty, we assume it points to an array of
	// LocalIntegral pointers, of length at least the number of elements in
	// the model (as defined by the highest number in the MLGE array).
	// If the array is shorter than this, expect a segmentation fault.
	LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[fe.iel-1];


	// --- Integration loop over all Gauss points in each direction --------

	int k1, k2, k3;
	for (int j = 0; j < nGauss; j++)
	  for (int i = 0; i < nGauss; i++)
	  {
	    // Local element coordinates of current integration point
	    xi[t0-1] = faceDir < 0 ? -1.0 : 1.0;
	    xi[t1-1] = xg[i];
	    xi[t2-1] = xg[j];
	    fe.xi   = xi[0];
	    fe.eta  = xi[1];
	    fe.zeta = xi[2];

	    // Local element coordinates and parameter values
	    // of current integration point
	    switch (abs(faceDir)) {
	    case 1: k2 = i; k3 = j; k1 = -1; break;
	    case 2: k1 = i; k3 = j; k2 = -1; break;
	    case 3: k1 = i; k2 = j; k3 = -1; break;
	    default: k1 = k2 = k3 = -1;
	    }
	    if (upar.size() > 1)
	      fe.u = 0.5*(upar[i1]*(1.0-xg[k1]) + upar[i1+1]*(1.0+xg[k1]));
	    if (vpar.size() > 1)
	      fe.v = 0.5*(vpar[i2]*(1.0-xg[k2]) + vpar[i2+1]*(1.0+xg[k2]));
	    if (wpar.size() > 1)
	      fe.w = 0.5*(wpar[i3]*(1.0-xg[k3]) + wpar[i3+1]*(1.0+xg[k3]));

	    // Compute the basis functions and their derivatives, using
	    // tensor product of one-dimensional Lagrange polynomials
	    if (!Lagrange::computeBasis(fe.N,dNdu,p1,xi[0],p2,xi[1],p3,xi[2]))
	      return false;

	    // Compute basis function derivatives and the face normal
	    fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,Xnod,dNdu,t1,t2);
	    if (fe.detJxW == 0.0) continue; // skip singular points

	    if (faceDir < 0) normal *= -1.0;

	    // Cartesian coordinates of current integration point
	    X = Xnod * fe.N;
	    X.t = time.t;

	    // Evaluate the integrand and accumulate element contributions
	    fe.detJxW *= wg[i]*wg[j];
	    if (!integrand.evalBou(elmInt,fe,time,X,normal))
	      return false;
	  }

	// Assembly of global system integral
	if (!glInt.assemble(elmInt,fe.iel))
	  return false;
      }

  return true;
}


bool ASMs3DLag::integrateEdge (Integrand& integrand, int lEdge,
			       GlobalIntegral& glInt,
			       const TimeDomain& time)
{
  if (!svol) return true; // silently ignore empty patches

  // Parametric direction of the edge {0, 1, 2}
  const int lDir = (lEdge-1)/4;

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Order of basis in the three parametric directions (order = degree + 1)
  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  // Number of elements in each direction
  const int nelx = (nx-1)/(p1-1);
  const int nely = (ny-1)/(p2-1);
  const int nelz = (nz-1)/(p3-1);

  FiniteElement fe(p1*p2*p3);
  Matrix dNdu, Xnod, Jac;
  Vec4   X;
  Vec3   tangent;
  double xi[3];

  switch (lEdge)
    {
    case  1: xi[1] = -1.0; xi[2] = -1.0; break;
    case  2: xi[1] =  1.0; xi[2] = -1.0; break;
    case  3: xi[1] = -1.0; xi[2] =  1.0; break;
    case  4: xi[1] =  1.0; xi[2] =  1.0; break;
    case  5: xi[0] = -1.0; xi[2] = -1.0; break;
    case  6: xi[0] =  1.0; xi[2] = -1.0; break;
    case  7: xi[0] = -1.0; xi[2] =  1.0; break;
    case  8: xi[0] =  1.0; xi[2] =  1.0; break;
    case  9: xi[0] = -1.0; xi[1] = -1.0; break;
    case 10: xi[0] =  1.0; xi[1] = -1.0; break;
    case 11: xi[0] = -1.0; xi[1] =  1.0; break;
    case 12: xi[0] =  1.0; xi[1] =  1.0; break;
    }


  // === Assembly loop over all elements on the patch edge =====================

  int iel = 1;
  for (int i3 = 0; i3 < nelz; i3++)
    for (int i2 = 0; i2 < nely; i2++)
      for (int i1 = 0; i1 < nelx; i1++, iel++)
      {
	fe.iel = MLGE[iel-1];

	// Skip elements that are not on current boundary edge
	bool skipMe = false;
	switch (lEdge)
	  {
	  case  1: if (i2 > 0      || i3 > 0)      skipMe = true; break;
	  case  2: if (i2 < nely-1 || i3 > 0)      skipMe = true; break;
	  case  3: if (i2 > 0      || i3 < nelz-1) skipMe = true; break;
	  case  4: if (i2 < nely-1 || i3 < nelz-1) skipMe = true; break;
	  case  5: if (i1 > 0      || i3 > 0)      skipMe = true; break;
	  case  6: if (i1 < nelx-1 || i3 > 0)      skipMe = true; break;
	  case  7: if (i1 > 0      || i3 < nelz-1) skipMe = true; break;
	  case  8: if (i1 < nelx-1 || i3 < nelz-1) skipMe = true; break;
	  case  9: if (i1 > 0      || i2 > 0)      skipMe = true; break;
	  case 10: if (i1 < nelx-1 || i2 > 0)      skipMe = true; break;
	  case 11: if (i1 > 0      || i2 < nely-1) skipMe = true; break;
	  case 12: if (i1 < nelx-1 || i2 < nely-1) skipMe = true; break;
	  }
	if (skipMe) continue;

	// Set up nodal point coordinates for current element
	if (!this->getElementCoordinates(Xnod,iel)) return false;

	// Initialize element quantities
	if (!integrand.initElementBou(MNPC[iel-1])) return false;


	// --- Integration loop over all Gauss points along the edge -----------

	LocalIntegral* elmInt = 0;
	for (int i = 0; i < nGauss; i++)
	{
	  // Gauss point coordinates on the edge
	  xi[lDir] = xg[i];

	  // Compute the basis functions and their derivatives, using
	  // tensor product of one-dimensional Lagrange polynomials
	  if (!Lagrange::computeBasis(fe.N,dNdu,p1,xi[0],p2,xi[1],p3,xi[2]))
	    return false;

	  // Compute basis function derivatives and the edge tangent
	  fe.detJxW = utl::Jacobian(Jac,tangent,fe.dNdX,Xnod,dNdu,1+lDir);
	  if (fe.detJxW == 0.0) continue; // skip singular points

	  // Cartesian coordinates of current integration point
	  X = Xnod * fe.N;
	  X.t = time.t;

	  // Evaluate the integrand and accumulate element contributions
	  if (!integrand.evalBou(elmInt,fe,time,X,tangent))
	    return false;
	}

	// Assembly of global system integral
	if (!glInt.assemble(elmInt,fe.iel))
	  return false;
      }

  return true;
}


int ASMs3DLag::evalPoint (const double* xi, double* param, Vec3& X) const
{
  if (!svol) return -3;

  // Evaluate the parametric values of the point and nodes
  RealArray u[3];
  for (int d = 0; d < 3; d++)
  {
    param[d] = (1.0-xi[d])*svol->startparam(d) + xi[d]*svol->endparam(d);
    if (!this->getGridParameters(u[d],d,svol->order(d)-1)) return -3;
  }

  // Search for the closest node
  size_t i = utl::find_closest(u[0],param[0]);
  size_t j = utl::find_closest(u[1],param[1]);
  size_t k = utl::find_closest(u[2],param[2]);
  size_t n = u[0].size()*(u[1].size()*k + j) + i;
  X = coord[n];

  return 1+n;
}


bool ASMs3DLag::tesselate (ElementBlock& grid, const int* npe) const
{
  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);
  if (p1 != npe[0] || p2 != npe[1] || p2 != npe[2])
  {
    int* newnpe = const_cast<int*>(npe);
    std::cout <<"\nLagrange elements: The number of visualization points are "
	      << p1 <<" "<< p2 <<" "<< p3 <<" by default\n"<< std::endl;
    newnpe[0] = p1;
    newnpe[1] = p2;
    newnpe[2] = p3;
  }

  return this->ASMs3D::tesselate(grid,npe);
}


bool ASMs3DLag::evalSolution (Matrix& sField, const Vector& locSol,
			      const int*) const
{
  size_t nPoints = coord.size();
  size_t nComp = locSol.size() / nPoints;
  if (nComp*nPoints != locSol.size())
    return false;

  size_t i, n, ip = 0;
  sField.resize(nComp,nPoints);
  for (n = 1; n <= nPoints; n++)
    for (i = 1; i <= nComp; i++)
      sField(i,n) = locSol(++ip);

  return true;
}


bool ASMs3DLag::evalSolution (Matrix&, const Vector&,
			      const RealArray*, bool) const
{
  std::cerr <<" *** ASMs3DLag::evalSolution(Matrix&,const Vector&,"
	    <<"const RealArray*,bool): Not implemented."<< std::endl;
  return false;
}


bool ASMs3DLag::evalSolution (Matrix& sField, const Integrand& integrand,
			      const int*, bool) const
{
  sField.resize(0,0);

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);
  double incx = 2.0/double(p1-1);
  double incy = 2.0/double(p2-1);
  double incz = 2.0/double(p3-1);

  size_t nPoints = coord.size();
  IntVec check(nPoints,0);

  Vector N(p1*p2*p3), solPt;
  std::vector<Vector> globSolPt(nPoints);
  Matrix dNdu, dNdX, Xnod, Jac;

  // Evaluate the secondary solution field at each point
  const int nel = this->getNoElms();
  for (int iel = 1; iel <= nel; iel++)
  {
    const IntVec& mnpc = MNPC[iel-1];
    this->getElementCoordinates(Xnod,iel);

    int i, j, k, loc = 0;
    for (k = 0; k < p3; k++)
      for (j = 0; j < p2; j++)
	for (i = 0; i < p1; i++, loc++)
	{
	  double xi   = -1.0 + i*incx;
	  double eta  = -1.0 + j*incy;
	  double zeta = -1.0 + k*incz;
	  if (!Lagrange::computeBasis(N,dNdu,p1,xi,p2,eta,p3,zeta))
	    return false;

	  // Compute the Jacobian inverse
	  if (utl::Jacobian(Jac,dNdX,Xnod,dNdu) == 0.0) // Jac = (Xnod*dNdu)^-1
	    continue; // skip singular points

	  // Now evaluate the solution field
	  if (!integrand.evalSol(solPt,N,dNdX,Xnod*N,mnpc))
	    return false;
	  else if (sField.empty())
	    sField.resize(solPt.size(),nPoints,true);

	  if (++check[mnpc[loc]] == 1)
	    globSolPt[mnpc[loc]] = solPt;
	  else
	    globSolPt[mnpc[loc]] += solPt;
	}
  }

  for (size_t i = 0; i < nPoints; i++)
    sField.fillColumn(1+i,globSolPt[i]/=check[i]);

  return true;
}


bool ASMs3DLag::evalSolution (Matrix&, const Integrand&,
			      const RealArray*, bool) const
{
  std::cerr <<" *** ASMs3DLag::evalSolution(Matrix&,const Integrand&,"
	    <<"const RealArray*,bool): Not implemented."<< std::endl;
  return false;
}
