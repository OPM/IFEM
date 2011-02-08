// $Id: ASMs1DLag.C,v 1.5 2010-12-29 18:39:37 kmo Exp $
//==============================================================================
//!
//! \file ASMs1DLag.C
//!
//! \date Apr 20 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for assembly of structured 1D Lagrange FE models.
//!
//==============================================================================

#include "GoTools/geometry/SplineCurve.h"

#include "ASMs1DLag.h"
#include "Lagrange.h"
#include "TimeDomain.h"
#include "GlobalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "ElementBlock.h"
#include "Utilities.h"
#include "Vec3Oper.h"


ASMs1DLag::ASMs1DLag (const char* fileName,
		      unsigned char n_s, unsigned char n_f)
  : ASMs1D(fileName,n_s,n_f)
{
  nx = 0;
}


ASMs1DLag::ASMs1DLag (std::istream& is, unsigned char n_s, unsigned char n_f)
  : ASMs1D(is,n_s,n_f)
{
  nx = 0;
}


void ASMs1DLag::clear ()
{
  coord.clear();
  nx = 0;
  ASMs1D::clear();
}


bool ASMs1DLag::generateFEMTopology ()
{
  if (!curv) return false;

  // Order of the basis
  const int p1 = curv->order();

  // Evaluate the parametric values
  RealArray gpar;
  if (!this->getGridParameters(gpar,p1-1)) return false;

  // Number of nodes
  nx = gpar.size();

  if (!coord.empty())
    return coord.size() == nx;

  // Number of elements
  const int nel = (nx-1)/(p1-1);

  MLGN.resize(nx);
  coord.resize(nx);

  // Evaluate the nodal coordinates
  for (int i = 0; i < nx; i++)
  {
    Go::Point pt;
    curv->point(pt,gpar[i]);
    for (int k = 0; k < pt.size(); k++)
      coord[i][k] = pt[k];
    MLGN[i] = ++gNod;
  }

  // Connectivety array: local --> global node relation
  MLGE.resize(nel);
  MNPC.resize(nel);

  int a, iel;
  for (iel = 0; iel < nel; iel++)
  {
    MLGE[iel] = ++gEl;
    // Element array
    MNPC[iel].resize(p1);
    // First node in current element
    int first = (p1-1)*iel;

    for (a = 0; a < p1; a++)
      MNPC[iel][a] = first + a;
  }

  return true;
}


bool ASMs1DLag::getElementCoordinates (Matrix& X, int iel) const
{
  if (iel < 1 || iel > MNPC.size())
  {
    std::cerr <<" *** ASMs1DLag::getElementCoordinates: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return false;
  }

  const IntVec& ien = MNPC[iel-1];
  X.resize(nsd,ien.size());

  for (size_t i = 0; i < ien.size(); i++)
    X.fillColumn(i+1,coord[ien[i]].ptr());

  return true;
}


void ASMs1DLag::getNodalCoordinates (Matrix& X) const
{
  X.resize(nsd,coord.size());

  for (size_t inod = 0; inod < coord.size(); inod++)
    X.fillColumn(inod+1,coord[inod].ptr());
}


Vec3 ASMs1DLag::getCoord (size_t inod) const
{
  if (inod < 1 || inod > coord.size()) return Vec3();

  return coord[inod-1];
}


bool ASMs1DLag::integrate (Integrand& integrand,
			   GlobalIntegral& glInt,
			   const TimeDomain& time,
			   const LintegralVec& locInt)
{
  if (!curv) return true; // silently ignore empty patches

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Order of basis (order = degree + 1)
  const int p1 = curv->order();

  Vector N(p1);
  Matrix dNdu, dNdX, Xnod, Jac;
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

      for (int i = 0; i < nGauss; i++)
	{
	  // Compute basis function derivatives at current integration point
	  if (!Lagrange::computeBasis(N,dNdu,p1,xg[i]))
	    return false;

	  // Compute Jacobian inverse of coordinate mapping and derivatives
	  double detJ = utl::Jacobian(Jac,dNdX,Xnod,dNdu);

	  // Cartesian coordinates of current integration point
	  X = Xnod * N;
	  X.t = time.t;

	  // Evaluate the integrand and accumulate element contributions
	  if (!integrand.evalInt(elmInt,time,detJ*wg[i],N,dNdX,X))
	    return false;
	}

      // Assembly of global system integral
      if (!glInt.assemble(elmInt,MLGE[iel-1]))
	return false;
    }

  return true;
}


bool ASMs1DLag::integrate (Integrand& integrand, int lIndex,
			   GlobalIntegral& glInt,
			   const TimeDomain& time,
			   const LintegralVec& locInt)
{
  if (!curv) return true; // silently ignore empty patches

  // Number of elements
  const int p1  = curv->order();
  const int nel = (nx-1)/(p1-1);

  // Integration of boundary point

  int iel;
  double x;
  switch (lIndex)
    {
    case 1:
      iel = 1;
      x = -1.0;
      break;

    case 2:
      iel = nel;
      x = 1.0;
    }

  // Set up control point coordinates for current element
  Matrix Xnod;
  if (!this->getElementCoordinates(Xnod,iel)) return false;

  // Initialize element quantities
  if (!integrand.initElementBou(MNPC[iel-1])) return false;

  // Caution: Unless locInt is empty, we assume it points to an array of
  // LocalIntegral pointers, of length at least the number of elements in
  // the model (as defined by the highest number in the MLGE array).
  // If the array is shorter than this, expect a segmentation fault.
  LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[MLGE[iel-1]-1];

  // Evaluate basis functions and corresponding derivatives
  Vector N;
  Matrix dNdu;
  if (!Lagrange::computeBasis(N,dNdu,p1,x)) return false;

  // Cartesian coordinates of current integration point
  Vec4 X(Xnod*N,time.t);

  // Compute basis function derivatives
  Matrix Jac, dNdX;
  utl::Jacobian(Jac,dNdX,Xnod,dNdu);

  // Set up the normal vector
  Vec3 normal;
  if (lIndex == 1)
    normal.x = -copysign(1.0,Jac(1,1));
  else
    normal.x = copysign(1.0,Jac(1,1));

  // Evaluate the integrand and accumulate element contributions
  if (!integrand.evalBou(elmInt,time,1.0,N,dNdX,X,normal))
    return false;

  // Assembly of global system integral
  return glInt.assemble(elmInt,MLGE[iel-1]);
}


bool ASMs1DLag::tesselate (ElementBlock& grid, const int*) const
{
  int i, l;

  grid.resize(nx);
  for (i = 0; i < grid.getNoNodes(); i++)
    grid.setCoor(i,coord[i].x,coord[i].y,coord[i].z);

  // Establish the block grid topology
  int n[2], ip = 0;
  n[0] = 0;
  n[1] = n[0] + 1;

  for (i = 1; i < nx; i++)
    for (l = 0; l < 2; l++)
      grid.setNode(ip++,n[l]++);

  return true;
}


bool ASMs1DLag::evalSolution (Matrix& sField, const Vector& locSol,
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


bool ASMs1DLag::evalSolution (Matrix& sField, const Integrand& integrand,
			      const int*) const
{
  sField.resize(0,0);

  const int p1 = curv->order();
  double incx = 2.0/double(p1-1);

  size_t nPoints = coord.size();
  IntVec check(nPoints,0);

  Vector N(p1), solPt;
  std::vector<Vector> globSolPt(nPoints);
  Matrix dNdu, dNdX, Xnod, Jac;

  // Evaluate the secondary solution field at each point
  for (int iel = 1; iel <= this->getNoElms(); iel++)
  {
    const IntVec& mnpc = MNPC[iel-1];
    this->getElementCoordinates(Xnod,iel);

    for (int loc = 0; loc < p1; loc++)
    {
      double xi = -1.0 + loc*incx;
      if (!Lagrange::computeBasis(N,dNdu,p1,xi))
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
