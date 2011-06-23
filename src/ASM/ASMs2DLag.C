// $Id$
//==============================================================================
//!
//! \file ASMs2DLag.C
//!
//! \date Mar 20 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for assembly of structured 2D Lagrange FE models.
//!
//==============================================================================

#include "GoTools/geometry/SplineSurface.h"

#include "ASMs2DLag.h"
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


ASMs2DLag::ASMs2DLag (const char* fileName,
		      unsigned char n_s, unsigned char n_f)
  : ASMs2D(fileName,n_s,n_f)
{
  nx = ny = 0;
}


ASMs2DLag::ASMs2DLag (std::istream& is, unsigned char n_s, unsigned char n_f)
  : ASMs2D(is,n_s,n_f)
{
  nx = ny = 0;
}


void ASMs2DLag::clear ()
{
  coord.clear();
  nx = ny = 0;
  ASMs2D::clear();
}


bool ASMs2DLag::generateFEMTopology ()
{
  if (!surf) return false;

  // Order of the basis in the two parametric directions (order = degree + 1)
  const int p1 = surf->order_u();
  const int p2 = surf->order_v();

  // Evaluate the parametric values
  RealArray gpar[2];
  if (!this->getGridParameters(gpar[0],0,p1-1)) return false;
  if (!this->getGridParameters(gpar[1],1,p2-1)) return false;

  // Number of nodes in each direction
  nx = gpar[0].size();
  ny = gpar[1].size();

  if (!coord.empty())
    return coord.size() == nx*ny;

  // Number of elements in each direction
  const int nelx = (nx-1)/(p1-1);
  const int nely = (ny-1)/(p2-1);

  // Evaluate the nodal coordinates in the physical space
  RealArray XYZ(surf->dimension()*nx*ny);
  surf->gridEvaluator(XYZ,gpar[0],gpar[1]);

  size_t i1, j1;
  MLGN.resize(nx*ny);
  coord.resize(nx*ny);
  for (i1 = j1 = 0; i1 < coord.size(); i1++)
  {
    MLGN[i1] = ++gNod;
    for (size_t d = 0; d < nsd; d++)
      coord[i1][d] = XYZ[j1+d];
    j1 += surf->dimension();
  }

  // Number of elements in patch
  const int nel = nelx*nely;
  // Number of nodes per element
  const int nen = p1*p2;

  // Connectivity array: local --> global node relation
  MLGE.resize(nel);
  MNPC.resize(nel);

  int i, j, a, b, iel = 0;
  for (j = 0; j < nely; j++)
    for (i = 0; i < nelx; i++, iel++)
    {
      MLGE[iel] = ++gEl;
      MNPC[iel].resize(nen);
      // First node in current element
      int corner = (p2-1)*nx*j + (p1-1)*i;

      for (b = 0; b < p2; b++)
      {
	int facenod = b*p1;
	MNPC[iel][facenod] = corner + b*nx;
	for (a = 1; a < p1; a++)
	  MNPC[iel][facenod+a] = MNPC[iel][facenod] + a;
      }
    }

  return true;
}


bool ASMs2DLag::getElementCoordinates (Matrix& X, int iel) const
{
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs2DLag::getElementCoordinates: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return false;
  }

  // Number of nodes per element
  const size_t nen = surf->order_u()*surf->order_v();

  X.resize(nsd,nen);
  for (size_t i = 0; i < nen; i++)
    X.fillColumn(i+1,coord[MNPC[iel-1][i]].ptr());

  return true;
}


void ASMs2DLag::getNodalCoordinates (Matrix& X) const
{
  X.resize(nsd,coord.size());

  for (size_t inod = 0; inod < coord.size(); inod++)
    X.fillColumn(inod+1,coord[inod].ptr());
}


Vec3 ASMs2DLag::getCoord (size_t inod) const
{
  if (inod < 1 || inod > coord.size()) return Vec3();

  return coord[inod-1];
}


bool ASMs2DLag::getSize (int& n1, int& n2, int) const
{
  n1 = nx;
  n2 = ny;

  return true;
}


bool ASMs2DLag::integrate (Integrand& integrand,
			   GlobalIntegral& glInt,
			   const TimeDomain& time,
			   const LintegralVec& locInt)
{
  if (!surf) return true; // silently ignore empty patches

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Get parametric coordinates of the elements
  RealArray upar, vpar;
  this->getGridParameters(upar,0,1);
  this->getGridParameters(vpar,1,1);

  // Number of elements in each direction
  const int nelx = upar.size() - 1;
  const int nely = vpar.size() - 1;

  // Order of basis in the two parametric directions (order = degree + 1)
  const int p1 = surf->order_u();
  const int p2 = surf->order_v();

  FiniteElement fe(p1*p2);
  Matrix dNdu, Xnod, Jac;
  Vec4   X;


  // === Assembly loop over all elements in the patch ==========================

  int iel = 1;
  for (int i2 = 0; i2 < nely; i2++)
    for (int i1 = 0; i1 < nelx; i1++, iel++)
    {
      // Set up control point coordinates for current element
      if (!this->getElementCoordinates(Xnod,iel)) return false;

      if (integrand.getIntegrandType() == 4)
      {
	// Compute the element "center" (average of element node coordinates)
	X = 0.0;
	for (size_t i = 1; i <= nsd; i++)
	  for (size_t j = 1; j <= Xnod.cols(); j++)
	    X[i-1] += Xnod(i,j);

	X *= 1.0/(double)Xnod.cols();
      }

      // Initialize element quantities
      if (!integrand.initElement(MNPC[iel-1],X,nGauss*nGauss)) return false;

      // Caution: Unless locInt is empty, we assume it points to an array of
      // LocalIntegral pointers, of length at least the number of elements in
      // the model (as defined by the highest number in the MLGE array).
      // If the array is shorter than this, expect a segmentation fault.
      LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[MLGE[iel-1]-1];


      // --- Integration loop over all Gauss points in each direction ----------

      for (int j = 0; j < nGauss; j++)
	for (int i = 0; i < nGauss; i++)
	{
	  // Parameter value of current integration point
	  fe.u = 0.5*(upar[i1]*(1.0-xg[i]) + upar[i1+1]*(1.0+xg[i]));
	  fe.v = 0.5*(vpar[i2]*(1.0-xg[j]) + vpar[i2+1]*(1.0+xg[j]));

	  // Local coordinates of current integration point
	  fe.xi  = xg[i];
	  fe.eta = xg[j];

	  // Compute basis function derivatives at current integration point
	  // using tensor product of one-dimensional Lagrange polynomials
	  if (!Lagrange::computeBasis(fe.N,dNdu,p1,xg[i],p2,xg[j]))
	    return false;

	  // Compute Jacobian inverse of coordinate mapping and derivatives
	  fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);
	  if (fe.detJxW == 0.0) continue; // skip singular points

	  // Cartesian coordinates of current integration point
	  X = Xnod * fe.N;
	  X.t = time.t;

	  // Evaluate the integrand and accumulate element contributions
	  fe.detJxW *= wg[i]*wg[j];
	  if (!integrand.evalInt(elmInt,fe,time,X))
	    return false;
	}

      // Finalize the element quantities
      if (!integrand.finalizeElement(elmInt,time))
	return false;

      // Assembly of global system integral
      if (!glInt.assemble(elmInt,MLGE[iel-1]))
	return false;
    }

  return true;
}


bool ASMs2DLag::integrate (Integrand& integrand, int lIndex,
			   GlobalIntegral& glInt,
			   const TimeDomain& time,
			   const LintegralVec& locInt)
{
  if (!surf) return true; // silently ignore empty patches

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Find the parametric direction of the edge normal {-2,-1, 1, 2}
  const int edgeDir = (lIndex+1)/(lIndex%2 ? -2 : 2);

  const int t1 = abs(edgeDir); // tangent direction normal to the patch edge
  const int t2 = 3-t1;         // tangent direction along the patch edge

  // Order of basis in the three parametric directions (order = degree + 1)
  const int p1 = surf->order_u();
  const int p2 = surf->order_v();

  // Number of elements in each direction
  const int nelx = (nx-1)/(p1-1);
  const int nely = (ny-1)/(p2-1);

  // Get parametric coordinates of the elements
  FiniteElement fe(p1*p2);
  RealArray upar, vpar;
  if (t1 == 1)
  {
    fe.u  = edgeDir < 0 ? surf->startparam_u() : surf->endparam_u();
    fe.xi = edgeDir < 0 ? -1.0 : 1.0; 
    this->getGridParameters(vpar,1,1);
  }
  else if (t1 == 2)
  {
    this->getGridParameters(upar,0,1);
    fe.v   = edgeDir < 0 ? surf->startparam_v() : surf->endparam_v();
    fe.eta = edgeDir < 0 ? -1.0 : 1.0;
  }

  Matrix dNdu, Xnod, Jac;
  Vec4   X;
  Vec3   normal;
  double xi[2];


  // === Assembly loop over all elements on the patch edge =====================

  int iel = 1;
  for (int i2 = 0; i2 < nely; i2++)
    for (int i1 = 0; i1 < nelx; i1++, iel++)
    {
      // Skip elements that are not on current boundary edge
      bool skipMe = false;
      switch (edgeDir)
	{
	case -1: if (i1 > 0)      skipMe = true; break;
	case  1: if (i1 < nelx-1) skipMe = true; break;
	case -2: if (i2 > 0)      skipMe = true; break;
	case  2: if (i2 < nely-1) skipMe = true; break;
	}
      if (skipMe) continue;

      // Set up control point coordinates for current element
      if (!this->getElementCoordinates(Xnod,iel)) return false;

      // Initialize element quantities
      if (!integrand.initElementBou(MNPC[iel-1])) return false;

      // Caution: Unless locInt is empty, we assume it points to an array of
      // LocalIntegral pointers, of length at least the number of elements in
      // the model (as defined by the highest number in the MLGE array).
      // If the array is shorter than this, expect a segmentation fault.
      LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[MLGE[iel-1]-1];


      // --- Integration loop over all Gauss points along the edge -------------

      for (int i = 0; i < nGauss; i++)
      {
	// Gauss point coordinates along the edge
	xi[t1-1] = edgeDir < 0 ? -1.0 : 1.0;
	xi[t2-1] = xg[i];

	// Parameter values and local coordinates of current integration point
	if (upar.size() > 1) {
	  fe.u = 0.5*(upar[i1]*(1.0-xg[i]) + upar[i1+1]*(1.0+xg[i]));
	  fe.xi = xg[i];
	}
	if (vpar.size() > 1) {
	  fe.v = 0.5*(vpar[i2]*(1.0-xg[i]) + vpar[i2+1]*(1.0+xg[i]));
	  fe.eta = xg[i];
	}

	// Compute the basis functions and their derivatives, using
	// tensor product of one-dimensional Lagrange polynomials
	if (!Lagrange::computeBasis(fe.N,dNdu,p1,xi[0],p2,xi[1]))
	  return false;

	// Compute basis function derivatives and the edge normal
	fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,Xnod,dNdu,t1,t2);
	if (fe.detJxW == 0.0) continue; // skip singular points

	if (edgeDir < 0) normal *= -1.0;

	// Cartesian coordinates of current integration point
	X = Xnod * fe.N;
	X.t = time.t;

	// Evaluate the integrand and accumulate element contributions
	fe.detJxW *= wg[i];
	if (!integrand.evalBou(elmInt,fe,time,X,normal))
	  return false;
      }

      // Assembly of global system integral
      if (!glInt.assemble(elmInt,MLGE[iel-1]))
	return false;
    }

  return true;
}


int ASMs2DLag::evalPoint (const double* xi, double* param, Vec3& X) const
{
  if (!surf) return -2;

  param[0] = (1.0-xi[0])*surf->startparam_u() + xi[0]*surf->endparam_u();
  param[1] = (1.0-xi[1])*surf->startparam_v() + xi[1]*surf->endparam_v();

  // Evaluate the parametric values of the nodes
  RealArray u, v;
  if (!this->getGridParameters(u,0,surf->order_u()-1)) return -2;
  if (!this->getGridParameters(v,1,surf->order_v()-1)) return -2;

  // Search for the closest node
  size_t i = utl::find_closest(u,param[0]);
  size_t j = utl::find_closest(v,param[1]);
  size_t n = u.size()*j + i;
  X = coord[n];

  return 1+n;
}


bool ASMs2DLag::tesselate (ElementBlock& grid, const int* npe) const
{
  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  if (p1 != npe[0] || p2 != npe[1])
  {
    int* newnpe = const_cast<int*>(npe);
    std::cout <<"\nLagrange elements: The number of visualization points are "
	      << p1 <<" "<< p2 <<" by default\n"<< std::endl;
    newnpe[0] = p1;
    newnpe[1] = p2;
  }

  return this->ASMs2D::tesselate(grid,npe);
}


bool ASMs2DLag::evalSolution (Matrix& sField, const Vector& locSol,
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


bool ASMs2DLag::evalSolution (Matrix&, const Vector&,
			      const RealArray*, bool) const
{
  std::cerr <<" *** ASMs2DLag::evalSolution(Matrix&,const Vector&,"
	    <<"const RealArray*,bool): Not implemented."<< std::endl;
  return false;
}


bool ASMs2DLag::evalSolution (Matrix& sField, const Integrand& integrand,
			      const int*, bool) const
{
  sField.resize(0,0);

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  double incx = 2.0/double(p1-1);
  double incy = 2.0/double(p2-1);

  size_t nPoints = coord.size();
  IntVec check(nPoints,0);

  Vector N(p1*p2), solPt;
  std::vector<Vector> globSolPt(nPoints);
  Matrix dNdu, dNdX, Xnod, Jac;

  // Evaluate the secondary solution field at each point
  const int nel = this->getNoElms();
  for (int iel = 1; iel <= nel; iel++)
  {
    const IntVec& mnpc = MNPC[iel-1];
    this->getElementCoordinates(Xnod,iel);

    int i, j, loc = 0;
    for (j = 0; j < p2; j++)
      for (i = 0; i < p1; i++, loc++)
      {
	double xi  = -1.0 + i*incx;
	double eta = -1.0 + j*incy;
	if (!Lagrange::computeBasis(N,dNdu,p1,xi,p2,eta))
	  return false;

	// Compute the Jacobian inverse
	if (utl::Jacobian(Jac,dNdX,Xnod,dNdu) == 0.0) // Jac = (Xnod * dNdu)^-1
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


bool ASMs2DLag::evalSolution (Matrix&, const Integrand&,
			      const RealArray*, bool) const
{
  std::cerr <<" *** ASMs2DLag::evalSolution(Matrix&,const Integrand&,"
	    <<"const RealArray*,bool): Not implemented."<< std::endl;
  return false;
}
