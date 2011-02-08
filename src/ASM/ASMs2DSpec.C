// $Id: ASMs2DSpec.C,v 1.6 2010-12-29 18:36:06 kmo Exp $
//==============================================================================
//!
//! \file ASMs2DSpec.C
//!
//! \date Mar 22 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for assembly of structured 2D Spectral FE models.
//!
//==============================================================================

#include "GoTools/geometry/SplineSurface.h"

#include "ASMs2DSpec.h"
#include "TimeDomain.h"
#include "GlobalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "Vec3Oper.h"
#include "Legendre.h"


bool ASMs2DSpec::getGridParameters (RealArray& prm, int dir,
				    int nSegPerSpan) const
{
  if (!surf) return false;

  // Evaluate the Gauss-Lobatto-Legendre points in this direction
  Vector dummy, xGLL;
  int p = dir == 0 ? surf->order_u() : surf->order_v();
  if (!Legendre::GLL(dummy,xGLL,p)) return false;

  if (xGLL.size() != nSegPerSpan+1)
  {
    nSegPerSpan = xGLL.size() - 1;
    std::cout <<"Spectral elements: Number of nodes per knot-span in "
	      << char('u'+dir) <<"-directon reset to "<< nSegPerSpan
	      <<" (GLL points)"<< std::endl;
  }

  std::vector<double>::const_iterator uit = surf->basis(dir).begin();
  double ucurr, uprev = *(uit++);
  while (uit != surf->basis(dir).end())
  {
    ucurr = *(uit++);
    if (ucurr > uprev)
      for (int i = 1; i <= nSegPerSpan; i++)
	prm.push_back(0.5*(ucurr-uprev)*(1.0+xGLL(i)) + uprev);
    uprev = ucurr;
  }

  prm.push_back(surf->basis(dir).endparam());
  return true;
}


/*!
  \brief Establishes matrices with basis functions and 1st derivatives.
*/

static void evalBasis (int i, int j, int p1, int p2,
		       const Matrix& der1, const Matrix& der2,
		       Vector& N, Matrix& dNdu)
{
  int a, b, n = 1;
  for (b = 1; b <= p2; b++)
    for (a = 1; a <= p1; a++, n++)
    {
       N  (n)   = (a == i && b == j) ? 1.0 : 0.0;
      dNdu(n,1) = (b == j) ? der1(i,a) : 0.0;
      dNdu(n,2) = (a == i) ? der2(j,b) : 0.0;
    }
}


bool ASMs2DSpec::integrate (Integrand& integrand,
			    GlobalIntegral& glInt,
			    const TimeDomain& time,
			    const LintegralVec& locInt)
{
  if (!surf) return true; // silently ignore empty patches

  // Order of basis in the two parametric directions (order = degree + 1)
  const int p1 = surf->order_u();
  const int p2 = surf->order_v();

  // Evaluate integration points (= nodal points) and weights
  Vector wg1,xg1,wg2,xg2;
  if (!Legendre::GLL(wg1,xg1,p1)) return false;
  if (!Legendre::GLL(wg2,xg2,p2)) return false;

  Matrix D1, D2;
  if (!Legendre::basisDerivatives(p1,D1)) return false;
  if (!Legendre::basisDerivatives(p2,D2)) return false;

  Vector  N (p1*p2);
  Matrix dNdu(p1*p2,2);
  Matrix dNdX, Xnod, Jac;
  Vec4   X;


  // === Assembly loop over all elements in the patch ==========================

  for (int iel = 1; iel <= this->getNoElms(); iel++)
    {
      // Set up control point coordinates for current element
      if (!this->getElementCoordinates(Xnod,iel)) return false;

      // Initialize element quantities
      if (!integrand.initElement(MNPC[iel-1])) return false;

      // Caution: Unless locInt is empty, we assume it points to an array of
      // LocalIntegral pointers, of length at least the number of elements in
      // the model (as defined by the highest number in the MLGE array).
      // If the array is shorter than this, expect a segmentation fault.
      LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[MLGE[iel-1]-1];


      // --- Integration loop over all Gauss points in each direction ----------

      int count = 1;
      for (int j = 1; j <= p2; j++)
	for (int i = 1; i <= p1; i++, count++)
	{
	  // Weight of current integration point
	  double weight = wg1(i)*wg2(j);

	  // Evaluate the basis functions and gradients using
	  // tensor product of one-dimensional Lagrange polynomials
	  evalBasis(i,j,p1,p2,D1,D2,N,dNdu);

	  // Compute Jacobian inverse of coordinate mapping and derivatives
	  double dA = utl::Jacobian(Jac,dNdX,Xnod,dNdu);
	  if (dA == 0.0) continue; // skip singular points

	  // Cartesian coordinates of current integration point
	  X.x = Xnod(1,count);
	  X.y = Xnod(2,count);
	  X.t = time.t;

	  // Evaluate the integrand and accumulate element contributions
	  if (!integrand.evalInt(elmInt,time,dA*weight,N,dNdX,X))
	    return false;
	}

      // Assembly of global system integral
      if (!glInt.assemble(elmInt,MLGE[iel-1]))
	return false;
    }

  return true;
}


bool ASMs2DSpec::integrate (Integrand& integrand, int lIndex,
			    GlobalIntegral& glInt,
			    const TimeDomain& time,
			    const LintegralVec& locInt)
{
  if (!surf) return true; // silently ignore empty patches

  // Find the parametric direction of the edge normal {-2,-1, 1, 2}
  int edgeDir = lIndex = (lIndex+1)/(lIndex%2 ? -2 : 2);

  const int t1 = abs(edgeDir);   // Tangent direction normal to the patch edge
  const int t2 = 3-abs(edgeDir); // Tangent direction along the patch edge

  // Order of basis in the two parametric directions (order = degree + 1)
  int p[2];
  p[0] = surf->order_u();
  p[1] = surf->order_v();

  // Number of elements in each direction
  int n1, n2;
  this->getSize(n1,n2);
  const int nelx = (n1-1)/(p[0]-1);
  const int nely = (n2-1)/(p[1]-1);

  // Evaluate integration points and weights

  Vector wg[2],xg1,xg2;
  if (!Legendre::GLL(wg[0],xg1,p[0])) return false;
  if (!Legendre::GLL(wg[1],xg2,p[1])) return false;

  Matrix D1, D2;
  if (!Legendre::basisDerivatives(p[0],D1)) return false;
  if (!Legendre::basisDerivatives(p[1],D2)) return false;

  int nen = p[0]*p[1];

  Vector N(nen);
  Matrix dNdu(nen,2);
  Matrix dNdX, Xnod, Jac;
  Vec4   X;
  Vec3   normal;
  int    xi[2];


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

      for (int i = 0; i < p[t2-1]; i++)
      {
	// Weight of current integration point
	double weight = wg[t2-1][i];

	// "Coordinates" along the edge
	xi[t1-1] = edgeDir < 0 ? 1 : p[t1-1];
	xi[t2-1] = i+1;

	// Evaluate the basis functions and gradients using
	// tensor product of one-dimensional Lagrange polynomials
	evalBasis(xi[0],xi[1],p[0],p[1],D1,D2,N,dNdu);

	// Compute basis function derivatives and the edge normal
	double dS = utl::Jacobian(Jac,normal,dNdX,Xnod,dNdu,t1,t2);
	if (dS == 0.0) continue; // skip singular points

	if (edgeDir < 0) normal *= -1.0;

	// Cartesian coordinates of current integration point
	X = Xnod * N;
	X.t = time.t;

	// Evaluate the integrand and accumulate element contributions
	if (!integrand.evalBou(elmInt,time,dS*weight,N,dNdX,X,normal))
	  return false;
      }

      // Assembly of global system integral
      if (!glInt.assemble(elmInt,MLGE[iel-1]))
	return false;
    }

  return true;
}


bool ASMs2DSpec::evalSolution (Matrix& sField, const Integrand& integrand,
			       const int*) const
{
  sField.resize(0,0);

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();

  Vector wg1,xg1,wg2,xg2;
  if (!Legendre::GLL(wg1,xg1,p1)) return false;
  if (!Legendre::GLL(wg2,xg2,p2)) return false;

  Matrix D1, D2;
  if (!Legendre::basisDerivatives(p1,D1)) return false;
  if (!Legendre::basisDerivatives(p2,D2)) return false;

  size_t nPoints = this->getNoNodes();
  IntVec check(nPoints,0);

  Vector N(p1*p2), solPt;
  std::vector<Vector> globSolPt(nPoints);
  Matrix dNdu(p1*p2,2), dNdX, Xnod, Jac;

  // Evaluate the secondary solution field at each point
  for (int iel = 1; iel <= this->getNoElms(); iel++)
  {
    const IntVec& mnpc = MNPC[iel-1];
    this->getElementCoordinates(Xnod,iel);

    int i, j, loc = 0;
    for (j = 0; j < p2; j++)
      for (i = 0; i < p1; i++, loc++)
      {
	evalBasis(i+1,j+1,p1,p2,D1,D2,N,dNdu);

	// Compute the Jacobian inverse
	if (utl::Jacobian(Jac,dNdX,Xnod,dNdu) == 0.0) // Jac = (Xnod * dNdu)^-1
	  continue; // skip singular points

	// Now evaluate the solution field
	if (!integrand.evalSol(solPt,N,dNdX,Xnod.getColumn(loc+1),mnpc))
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
