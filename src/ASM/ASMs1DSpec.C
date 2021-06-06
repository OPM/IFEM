// $Id$
//==============================================================================
//!
//! \file ASMs1DSpec.C
//!
//! \date Mar 22 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for assembly of structured 1D Spectral FE models.
//!
//==============================================================================

#include "GoTools/geometry/SplineCurve.h"

#include "ASMs1DSpec.h"
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "Vec3Oper.h"
#include "Legendre.h"
#include "Lagrange.h"


bool ASMs1DSpec::getGridParameters (RealArray& prm, int nSegPerSpan) const
{
  if (!curv) return false;

  // Evaluate the Gauss-Lobatto-Legendre points
  Vector dummy, xGLL;
  if (!Legendre::GLL(dummy,xGLL,curv->order())) return false;

  if (xGLL.size() != (size_t)(nSegPerSpan+1))
  {
    nSegPerSpan = xGLL.size() - 1;
    std::cout <<"Spectral elements: Number of nodes per knot-span reset to "
	      << nSegPerSpan <<" (GLL points)"<< std::endl;
  }

  RealArray::const_iterator uit = curv->basis().begin();
  double uprev = *(uit++);
  while (uit != curv->basis().end())
  {
    double ucurr = *(uit++);
    if (ucurr > uprev)
      for (int i = 1; i <= nSegPerSpan; i++)
	prm.push_back(0.5*(ucurr-uprev)*(1.0+xGLL(i)) + uprev);
    uprev = ucurr;
  }

  prm.push_back(curv->basis().endparam());
  return true;
}


bool ASMs1DSpec::integrate (Integrand& integrand,
			    GlobalIntegral& glInt,
			    const TimeDomain& time)
{
  if (this->empty()) return true; // silently ignore empty patches

  // Evaluate integration points and weights

  Vector wg1, xg1, points1;
  Matrix D1, dNdu, Jac;
  bool ok = Legendre::GLL(wg1,points1,p1);
  int nGP = nGauss;

  if (nGauss < 1) // using the nodal points themselves as integration points
    ok &= Legendre::basisDerivatives(nGP=p1,D1);
  else // using Gauss-Legendre scheme with nGauss integration points
    ok &= Legendre::GL(wg1,xg1,nGP=nGauss);

  FiniteElement fe(p1);
  Vec4 X(Vec3(),time.t);


  // === Assembly loop over all elements in the patch ==========================

  for (size_t iel = 1; iel <= nel && ok; iel++)
  {
    fe.iel = MLGE[iel-1];
    LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel);
    if (!A) continue; // no integrand contributions for this element

    // Set up control point coordinates for current element
    ok = this->getElementCoordinates(fe.Xn,iel);

    // Initialize element quantities
    ok &= integrand.initElement(MNPC[iel-1],*A);

    // --- Integration loop over integration points ----------------------------

    for (int i = 0; i < nGP && ok; i++)
    {
      // Compute basis function derivatives at current integration point
      if (nGauss < 1)
      {
	fe.N.fill(0.0);
	fe.N(i+1) = 1.0;
	dNdu.fillColumn(1,D1.getRow(i+1));
      }
      else
	ok = Lagrange::computeBasis(fe.N,&dNdu,points1,xg1[i]);

      // Compute Jacobian inverse of coordinate mapping and derivatives
      fe.detJxW = utl::Jacobian(Jac,fe.dNdX,fe.Xn,dNdu);

      // Cartesian coordinates of current integration point
      X.assign(fe.Xn*fe.N);

      // Evaluate the integrand and accumulate element contributions
      fe.detJxW *= wg1[i];
      ok &= integrand.evalInt(*A,fe,time,X);
    }

    // Assembly of global system integral
    if (ok && !glInt.assemble(A->ref(),fe.iel))
      ok = false;

    A->destruct();
  }

  return ok;
}


bool ASMs1DSpec::integrate (Integrand& integrand, int lIndex,
			    GlobalIntegral& glInt,
			    const TimeDomain& time)
{
  return false; // not implemented
}


bool ASMs1DSpec::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			       const RealArray*, bool) const
{
  sField.resize(0,0);

  Matrix D1;
  if (!Legendre::basisDerivatives(p1,D1))
    return false;

  size_t nPoints = this->getNoNodes();
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

    for (int i = 0; i < p1; i++)
    {
      fe.N.fill(0.0);
      fe.N(i+1) = 1.0;
      dNdu.fillColumn(1,D1.getRow(i+1));

      // Compute the Jacobian inverse
      fe.detJxW = utl::Jacobian(Jac,fe.dNdX,fe.Xn,dNdu);

      // Now evaluate the solution field
      if (!integrand.evalSol(solPt,fe,fe.Xn.getColumn(i+1),mnpc))
	return false;
      else if (sField.empty())
	sField.resize(solPt.size(),nPoints,true);

      if (++check[mnpc[i]] == 1)
	globSolPt[mnpc[i]] = solPt;
      else
	globSolPt[mnpc[i]] += solPt;
    }
  }

  for (size_t i = 0; i < nPoints; i++)
    sField.fillColumn(1+i,globSolPt[i] /= check[i]);

  return true;
}
