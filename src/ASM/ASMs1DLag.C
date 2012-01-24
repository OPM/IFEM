// $Id$
//==============================================================================
//!
//! \file ASMs1DLag.C
//!
//! \date Apr 20 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for assembly of 1D Lagrange FE models.
//!
//==============================================================================

#include "GoTools/geometry/SplineCurve.h"

#include "ASMs1DLag.h"
#include "Lagrange.h"
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "ElementBlock.h"
#include "Utilities.h"
#include "Vec3Oper.h"


ASMs1DLag::ASMs1DLag (unsigned char n_s, unsigned char n_f)
  : ASMs1D(n_s,n_f), coord(myCoord)
{
  nx = 0;
}


ASMs1DLag::ASMs1DLag (const ASMs1DLag& patch, unsigned char n_f)
  : ASMs1D(patch,n_f), coord(patch.myCoord)
{
  nx = coord.size();
}


void ASMs1DLag::clear (bool retainGeometry)
{
  myCoord.clear();
  nx = 0;

  this->ASMs1D::clear(retainGeometry);
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

  myMLGN.resize(nx);
  myCoord.resize(nx);

  // Evaluate the nodal coordinates
  Go::Point pt;
  for (size_t i = 0; i < nx; i++)
  {
    curv->point(pt,gpar[i]);
    for (int k = 0; k < pt.size(); k++)
      myCoord[i][k] = pt[k];
    myMLGN[i] = ++gNod;
  }

  // Connectivity array: local --> global node relation
  myMLGE.resize(nel);
  myMNPC.resize(nel);

  int a, iel;
  for (iel = 0; iel < nel; iel++)
  {
    myMLGE[iel] = ++gEl;
    // Element array
    myMNPC[iel].resize(p1);
    // First node in current element
    int first = (p1-1)*iel;

    for (a = 0; a < p1; a++)
      myMNPC[iel][a] = first + a;
  }

  return true;
}


Vec3 ASMs1DLag::getCoord (size_t inod) const
{
  if (inod < 1 || inod > coord.size()) return Vec3();

  return coord[inod-1];
}


bool ASMs1DLag::getElementCoordinates (Matrix& X, int iel) const
{
  if (iel < 1 || (size_t)iel > MNPC.size())
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


bool ASMs1DLag::updateCoords (const Vector& displ)
{
  if (shareFE) return true;

  if (displ.size() != nsd*coord.size())
  {
    std::cerr <<" *** ASMs1DLag::updateCoords: Invalid dimension "
	      << displ.size() <<" on displ, should be "
	      << nsd*coord.size() << std::endl;
    return false;
  }

  const double* u = displ.ptr();
  for (size_t inod = 0; inod < myCoord.size(); inod++, u += nsd)
    myCoord[inod] += RealArray(u,u+nsd);

  return true;
}


bool ASMs1DLag::integrate (Integrand& integrand,
			   GlobalIntegral& glInt,
			   const TimeDomain& time)
{
  if (!curv) return true; // silently ignore empty patches

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
  RealArray gpar;
  this->getGridParameters(gpar,1);

  // Order of basis (order = degree + 1)
  const int p1 = curv->order();

  FiniteElement fe(p1);
  Matrix dNdu, Xnod, Jac;
  Vec4   X;


  // === Assembly loop over all elements in the patch ==========================

  const int nel = this->getNoElms();
  for (int iel = 1; iel <= nel; iel++)
  {
    // Set up nodal point coordinates for current element
    if (!this->getElementCoordinates(Xnod,iel)) return false;

    // Initialize element quantities
    fe.iel = MLGE[iel-1];
    LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel);
    if (!integrand.initElement(MNPC[iel-1],X,nRed,*A)) return false;

    if (integrand.getIntegrandType() > 10)

      // --- Selective reduced integration loop --------------------------------

      for (int i = 0; i < nRed; i++)
      {
	// Local element coordinates of current integration point
	fe.xi  = xr[i];

	// Parameter value of current integration point
	fe.u = 0.5*(gpar[iel-1]*(1.0-xr[i]) + gpar[iel]*(1.0+xr[i]));

	// Compute basis function derivatives at current integration point
	if (!Lagrange::computeBasis(fe.N,dNdu,p1,xr[i]))
	  return false;

	// Compute Jacobian inverse and derivatives
	fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

	// Cartesian coordinates of current integration point
	X = Xnod * fe.N;
	X.t = time.t;

	// Compute the reduced integration terms of the integrand
	if (!integrand.reducedInt(*A,fe,X))
	  return false;
      }


    // --- Integration loop over all Gauss points in each direction ------------

    int jp = (iel-1)*nGauss;
    fe.iGP = firstIp + jp; // Global integration point counter

    for (int i = 0; i < nGauss; i++, fe.iGP++)
    {
      // Local element coordinate of current integration point
      fe.xi = xg[i];

      // Parameter value of current integration point
      fe.u = 0.5*(gpar[iel-1]*(1.0-xg[i]) + gpar[iel]*(1.0+xg[i]));

      // Compute basis function derivatives at current integration point
      if (!Lagrange::computeBasis(fe.N,dNdu,p1,xg[i]))
	return false;

      // Compute Jacobian inverse of coordinate mapping and derivatives
      fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu)*wg[i];

      // Cartesian coordinates of current integration point
      X = Xnod * fe.N;
      X.t = time.t;

      // Evaluate the integrand and accumulate element contributions
      if (!integrand.evalInt(*A,fe,time,X))
	return false;
    }

    // Finalize the element quantities
    if (!integrand.finalizeElement(*A,time,firstIp+jp))
      return false;

    // Assembly of global system integral
    if (!glInt.assemble(A->ref(),fe.iel))
      return false;

    A->destruct();
  }

  return true;
}


bool ASMs1DLag::integrate (Integrand& integrand, int lIndex,
			   GlobalIntegral& glInt,
			   const TimeDomain& time)
{
  if (!curv) return true; // silently ignore empty patches

  // Integration of boundary point

  FiniteElement fe(curv->order());
  int iel;
  switch (lIndex)
    {
    case 1:
      fe.xi = -1.0;
      fe.u = curv->startparam();
      iel = 1;
      break;

    case 2:
      fe.xi = 1.0;
      fe.u = curv->endparam();
      iel = this->getNoElms();

    default:
      return false;
    }

  // Set up nodal point coordinates for current element
  Matrix Xnod;
  if (!this->getElementCoordinates(Xnod,iel)) return false;

  // Initialize element quantities
  fe.iGP = firstBp[lIndex];
  fe.iel = MLGE[iel-1];
  LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel,true);
  if (!integrand.initElementBou(MNPC[iel-1],*A)) return false;

  // Evaluate basis functions and corresponding derivatives
  Matrix dNdu, Jac;
  if (!Lagrange::computeBasis(fe.N,dNdu,curv->order(),fe.xi)) return false;

  // Cartesian coordinates of current integration point
  Vec4 X(Xnod*fe.N,time.t);

  // Compute basis function derivatives
  utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

  // Set up the normal vector
  Vec3 normal;
  if (lIndex == 1)
    normal.x = -copysign(1.0,Jac(1,1));
  else
    normal.x = copysign(1.0,Jac(1,1));

  // Evaluate the integrand and accumulate element contributions
  if (!integrand.evalBou(*A,fe,time,X,normal))
    return false;

  // Assembly of global system integral
  bool result = glInt.assemble(A->ref(),fe.iel);

  A->destruct();
  return result;
}


int ASMs1DLag::evalPoint (const double* xi, double* param, Vec3& X) const
{
  if (!curv) return -1;

  param[0] = (1.0-xi[0])*curv->startparam() + xi[0]*curv->endparam();

  // Evaluate the parametric values of the nodes
  RealArray u;
  if (!this->getGridParameters(u,curv->order()-1)) return -1;

  // Search for the closest node
  size_t n = utl::find_closest(u,param[0]);
  X = coord[n];

  return 1+n;
}


bool ASMs1DLag::tesselate (ElementBlock& grid, const int* npe) const
{
  const int p1 = curv->order();
  if (p1 != npe[0])
  {
    std::cout <<"\nLagrange elements: The number of visualization points is "
	      << p1 <<" by default\n"<< std::endl;
    const_cast<int*>(npe)[0] = p1;
  }

  size_t i, l;

  grid.resize(nx);
  for (i = 0; i < grid.getNoNodes(); i++)
    grid.setCoor(i,coord[i].x,coord[i].y,coord[i].z);

  // Establish the block grid topology
  int n[2], ie = 1, ip = 0;
  n[0] = 0;
  n[1] = n[0] + 1;

  for (i = 1; i < nx; i++)
  {
    for (l = 0; l < 2; l++)
      grid.setNode(ip++,n[l]++);
    grid.setElmId(i,ie);
    if (i%(p1-1) == 0) ie++;
  }

  return true;
}


bool ASMs1DLag::evalSolution (Matrix& sField, const Vector& locSol,
			      const int*) const
{
  return this->evalSolution(sField,locSol,(const RealArray*)0,true);
}


bool ASMs1DLag::evalSolution (Matrix& sField, const Vector& locSol,
			      const RealArray*, bool) const
{
  size_t nPoints = coord.size();
  size_t nComp = locSol.size() / nPoints;
  if (nComp*nPoints != locSol.size())
    return false;

  sField.resize(nComp,nPoints);
  const double* u = locSol.ptr();
  for (size_t n = 1; n <= nPoints; n++, u += nComp)
    sField.fillColumn(n,u);

  return true;
}


bool ASMs1DLag::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			      const int*, bool) const
{
  return this->evalSolution(sField,integrand,(const RealArray*)0,true);
}


bool ASMs1DLag::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			      const RealArray*, bool) const
{
  sField.resize(0,0);
  if (!curv) return false;

  const int p1 = curv->order();
  double incx = 2.0/double(p1-1);

  size_t nPoints = coord.size();
  IntVec check(nPoints,0);

  Vector N(p1), solPt;
  std::vector<Vector> globSolPt(nPoints);
  Matrix dNdu, dNdX, Xnod, Jac;

  // Evaluate the secondary solution field at each point
  const int nel = this->getNoElms();
  for (int iel = 1; iel <= nel; iel++)
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
