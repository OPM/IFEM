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


ASMs1DLag::ASMs1DLag (unsigned char n_s, unsigned char n_f)
  : ASMs1D(n_s,n_f), coord(myCoord)
{
  nx = 0;
  p1 = 0;
}


ASMs1DLag::ASMs1DLag (const ASMs1DLag& patch, unsigned char n_f)
  : ASMs1D(patch,n_f), coord(patch.myCoord)
{
  nx = patch.nx;
  p1 = patch.p1;
}


ASMs1DLag::ASMs1DLag (const ASMs1DLag& patch)
  : ASMs1D(patch), coord(myCoord), myCoord(patch.coord)
{
  nx = patch.nx;
  p1 = patch.p1;
}


void ASMs1DLag::clear (bool retainGeometry)
{
  myCoord.clear();
  nx = 0;
  p1 = 0;

  this->ASMs1D::clear(retainGeometry);
}


bool ASMs1DLag::generateOrientedFEModel (const Vec3& Zaxis)
{
  if (!curv) return false;
  if (!proj) proj = curv;

  // Order of the basis
  p1 = curv->order();

  // Evaluate the parametric values
  RealArray gpar;
  if (!this->getGridParameters(gpar,p1-1)) return false;

  // Number of nodes in the patch
  nnod = nx = gpar.size();

  if (!coord.empty())
    return coord.size() == nnod;

  // Number of elements in the patch
  nel = (nx-1)/(p1-1);

  myMLGN.resize(nnod);
  myCoord.resize(nnod);
  if (nsd == 3 && nf == 6)
  {
    // This is a 3D beam problem, allocate the nodal/element rotation tensors.
    // The nodal rotations are updated during the simulation according to the
    // deformation state, whereas the element tensors are kept constant.
    myCS.resize(nel,Tensor(3));
    myT.resize(nnod,Tensor(3,true)); // Initialize nodal rotations to unity
  }

  // Evaluate the nodal coordinates
  Go::Point pt;
  for (size_t i = 0; i < nnod; i++)
  {
    curv->point(pt,gpar[i]);
    for (int k = 0; k < pt.size(); k++)
      myCoord[i][k] = pt[k];
    myMLGN[i] = ++gNod;
  }

  // Connectivity array: local --> global node relation
  myMLGE.resize(nel);
  myMNPC.resize(nel);

  for (size_t iel = 0; iel < nel; iel++)
  {
    myMLGE[iel] = ++gEl;
    // Element array
    myMNPC[iel].resize(p1);
    // First node in current element
    myMNPC[iel].front() = (p1-1)*iel;

    for (int a = 1; a < p1; a++)
      myMNPC[iel][a] = myMNPC[iel][a-1] + 1;
  }

  return myCS.empty() ? true : this->initLocalElementAxes(Zaxis);
}


Vec3 ASMs1DLag::getCoord (size_t inod) const
{
  if (inod < 1 || inod > coord.size()) return Vec3();

  return coord[inod-1];
}


bool ASMs1DLag::getElementCoordinates (Matrix& X, int iel, bool) const
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


/*!
  \brief Extracts the element end points from the element coordinates.
*/

static double getEndPoints (const Matrix& Xnod, Vec3Vec& XC)
{
  XC.resize(2);
  for (size_t i = 0; i < Xnod.rows() && i < 3; i++)
  {
    XC[0][i] = Xnod(i+1,1);
    XC[1][i] = Xnod(i+1,Xnod.cols());
  }

  return ASM1D::getElementSize(XC);
}


bool ASMs1DLag::integrate (Integrand& integrand,
			   GlobalIntegral& glInt,
			   const TimeDomain& time)
{
  if (this->empty()) return true; // silently ignore empty patches

  // Get Gaussian quadrature points and weights
  const int     ng = this->getNoGaussPt(p1);
  const double* xg = GaussQuadrature::getCoord(ng);
  const double* wg = GaussQuadrature::getWeight(ng);
  if (!xg || !wg) return false;

  // Get the reduced integration quadrature points, if needed
  const double* xr = nullptr;
  const double* wr = nullptr;
  int nRed = integrand.getReducedIntegration(ng);
  if (nRed > 0)
  {
    xr = GaussQuadrature::getCoord(nRed);
    wr = GaussQuadrature::getWeight(nRed);
    if (!xr || !wr) return false;
  }
  else if (nRed < 0)
    nRed = ng; // The integrand needs to know nGauss

  // Get parametric coordinates of the elements
  RealArray gpar;
  this->getGridParameters(gpar,1);

  FiniteElement fe(p1);
  Matrix dNdu, Jac;
  Vec4   X;


  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t iel = 0; iel < nel && ok; iel++)
  {
    fe.iel = MLGE[iel];
    LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel);
    if (!A) continue; // no integrand contributions for this element

    // Set up nodal point coordinates for current element
    ok = this->getElementCoordinates(fe.Xn,1+iel);

    if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
      fe.h = getEndPoints(fe.Xn,fe.XC);

    if (integrand.getIntegrandType() & Integrand::NODAL_ROTATIONS)
    {
      this->getElementNodalRotations(fe.Tn,iel);
      if (!elmCS.empty()) fe.Te = elmCS[iel];
    }

    // Initialize element quantities
    ok &= integrand.initElement(MNPC[iel],fe,X,nRed,*A);

    if (xr)
    {
      // --- Selective reduced integration loop --------------------------------

      for (int i = 0; i < nRed && ok; i++)
      {
	// Local element coordinates of current integration point
	fe.xi  = xr[i];

	// Parameter value of current integration point
        if (!gpar.empty())
          fe.u = 0.5*(gpar[iel]*(1.0-xr[i]) + gpar[iel+1]*(1.0+xr[i]));

        if (integrand.getIntegrandType() & Integrand::NO_DERIVATIVES)
          ok = Lagrange::computeBasis(fe.N,p1,xr[i]);
        else
        {
          // Compute basis function derivatives at current integration point
          ok = Lagrange::computeBasis(fe.N,dNdu,p1,xr[i]);
          // Compute Jacobian inverse and derivatives
          fe.detJxW = utl::Jacobian(Jac,fe.dNdX,fe.Xn,dNdu)*wr[i];
	}

	// Cartesian coordinates of current integration point
	X = fe.Xn * fe.N;
	X.t = time.t;

	// Compute the reduced integration terms of the integrand
	ok &= integrand.reducedInt(*A,fe,X);
      }
    }


    // --- Integration loop over all Gauss points in each direction ------------

    int jp = iel*ng;
    fe.iGP = firstIp + jp; // Global integration point counter

    for (int i = 0; i < ng && ok; i++, fe.iGP++)
    {
      // Local element coordinate of current integration point
      fe.xi = xg[i];

      // Parameter value of current integration point
      if (!gpar.empty())
        fe.u = 0.5*(gpar[iel]*(1.0-xg[i]) + gpar[iel+1]*(1.0+xg[i]));

      if (integrand.getIntegrandType() & Integrand::NO_DERIVATIVES)
        ok = Lagrange::computeBasis(fe.N,p1,xg[i]);
      else
      {
        // Compute basis function derivatives at current integration point
        ok = Lagrange::computeBasis(fe.N,dNdu,p1,xg[i]);
        // Compute Jacobian inverse of coordinate mapping and derivatives
        fe.detJxW = utl::Jacobian(Jac,fe.dNdX,fe.Xn,dNdu)*wg[i];
      }

      // Cartesian coordinates of current integration point
      X = fe.Xn * fe.N;
      X.t = time.t;

      // Evaluate the integrand and accumulate element contributions
      ok &= integrand.evalInt(*A,fe,time,X);
    }

    // Finalize the element quantities
    if (ok && !integrand.finalizeElement(*A,fe,time,firstIp+jp))
      ok = false;

    // Assembly of global system integral
    if (ok && !glInt.assemble(A->ref(),fe.iel))
      ok = false;

    A->destruct();
  }

  return ok;
}


bool ASMs1DLag::integrate (Integrand& integrand, int lIndex,
			   GlobalIntegral& glInt,
			   const TimeDomain& time)
{
  if (this->empty()) return true; // silently ignore empty patches

  // Extract the Neumann order flag (1 or higher) for the integrand
  integrand.setNeumannOrder(1 + lIndex/10);

  // Integration of boundary point

  FiniteElement fe(p1);
  size_t iel = 0;
  switch (lIndex%10)
    {
    case 1:
      fe.xi = -1.0;
      fe.u = curv ? curv->startparam() : 0.0;
      break;

    case 2:
      fe.xi = 1.0;
      fe.u = curv ? curv->endparam() : 1.0;
      iel = nel-1;
      break;

    default:
      return false;
    }

  fe.iel = MLGE[iel];
  LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel,true);
  if (!A) return true; // no integrand contributions for this element

  // Set up nodal point coordinates for current element
  bool ok = this->getElementCoordinates(fe.Xn,1+iel);

  if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
    fe.h = getEndPoints(fe.Xn,fe.XC);

  if (integrand.getIntegrandType() & Integrand::NODAL_ROTATIONS)
  {
    this->getElementNodalRotations(fe.Tn,iel);
    if (!elmCS.empty()) fe.Te = elmCS[iel];
  }

  // Initialize element quantities
  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex%10);
  fe.iGP = iit == firstBp.end() ? 0 : iit->second;
  ok &= integrand.initElementBou(MNPC[iel],*A);

  // Evaluate basis functions and corresponding derivatives
  Vec3 normal;
  if (integrand.getIntegrandType() & Integrand::NO_DERIVATIVES)
    ok &= Lagrange::computeBasis(fe.N,p1,fe.xi);
  else
  {
    // Compute basis function derivatives
    Matrix dNdu, Jac;
    ok &= Lagrange::computeBasis(fe.N,dNdu,p1,fe.xi);
    utl::Jacobian(Jac,fe.dNdX,fe.Xn,dNdu);

    // Set up the normal vector
    if (lIndex%10 == 1)
      normal.x = -copysign(1.0,Jac(1,1));
    else
      normal.x = copysign(1.0,Jac(1,1));
  }

  // Cartesian coordinates of current integration point
  Vec4 X(fe.Xn*fe.N,time.t);

  // Evaluate the integrand and accumulate element contributions
  if (ok && !integrand.evalBou(*A,fe,time,X,normal))
    ok = false;

  // Assembly of global system integral
  if (ok && !glInt.assemble(A->ref(),fe.iel))
    ok = false;

  A->destruct();
  return ok;
}


int ASMs1DLag::evalPoint (const double* xi, double* param, Vec3& X) const
{
  if (curv)
    param[0] = (1.0-xi[0])*curv->startparam() + xi[0]*curv->endparam();
  else
    param[0] = xi[0];

  // Evaluate the parametric values of the nodes
  RealArray u;
  if (!this->getGridParameters(u,p1-1)) return -1;

  // Search for the closest node
  size_t n = utl::find_closest(u,param[0]);
  X = coord[n];

  return 1+n;
}


bool ASMs1DLag::tesselate (ElementBlock& grid, const int* npe) const
{
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
                              const int*, int) const
{
  return this->evalSolution(sField,locSol,nullptr,false,0,0);
}


bool ASMs1DLag::evalSolution (Matrix& sField, const Vector& locSol,
                              const RealArray*, bool, int, int) const
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
			      const int*, char) const
{
  return this->evalSolution(sField,integrand,(const RealArray*)nullptr);
}


bool ASMs1DLag::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			      const RealArray*, bool) const
{
  sField.resize(0,0);

  double incx = 2.0/double(p1-1);

  size_t nPoints = coord.size();
  IntVec check(nPoints,0);

  FiniteElement fe(p1);
  Vector        solPt;
  Vectors       globSolPt(nPoints);
  Matrix        dNdu, Jac;

  // Evaluate the secondary solution field at each point
  for (size_t iel = 0; iel < nel; iel++)
  {
    const IntVec& mnpc = MNPC[iel];
    this->getElementCoordinates(fe.Xn,1+iel);

    for (int loc = 0; loc < p1; loc++)
    {
      fe.xi = -1.0 + loc*incx;
      if (integrand.getIntegrandType() & Integrand::NO_DERIVATIVES)
      {
        if (!Lagrange::computeBasis(fe.N,p1,fe.xi))
          return false;
      }
      else
      {
        if (!Lagrange::computeBasis(fe.N,dNdu,p1,fe.xi))
          return false;

        // Compute the Jacobian inverse
        fe.detJxW = utl::Jacobian(Jac,fe.dNdX,fe.Xn,dNdu);
      }

      // Now evaluate the solution field
      if (!integrand.evalSol(solPt,fe,fe.Xn*fe.N,mnpc))
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
    sField.fillColumn(1+i,globSolPt[i] /= check[i]);

  return true;
}


bool ASMs1DLag::write (std::ostream& os, int) const
{
  return this->writeLagBasis(os,"line");
}
