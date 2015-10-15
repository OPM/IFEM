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
  double ucurr, uprev = *(uit++);
  while (uit != curv->basis().end())
  {
    ucurr = *(uit++);
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
  if (!curv) return true; // silently ignore empty patches

  // Order of basis (order = degree + 1)
  const int p1 = curv->order();
  const int n1 = nGauss < 1 ? p1 : nGauss;

  // Evaluate integration points and weights

  Vector wg1, xg1, points1;
  if (!Legendre::GLL(wg1,points1,p1))
    return false;

  Matrix D1;
  if (nGauss < 1)
  {
    // We are using the nodal points themselves as integration points
    if (!Legendre::basisDerivatives(n1,D1))
      return false;
  }
  else
    // Using Gauss-Legendre scheme with nGauss points
    if (!Legendre::GL(wg1,xg1,n1))
      return false;

  FiniteElement fe(p1);
  Matrix dNdu, Xnod, Jac;
  Vec4   X;


  // === Assembly loop over all elements in the patch ==========================

  const int nel = this->getNoElms();
  for (int iel = 1; iel <= nel; iel++)
  {
    // Set up control point coordinates for current element
    if (!this->getElementCoordinates(Xnod,iel)) return false;

    // Initialize element quantities
    fe.iel = MLGE[iel-1];
    LocalIntegral* A = integrand.getLocalIntegral(p1,fe.iel);
    if (!integrand.initElement(MNPC[iel-1],*A)) return false;

    // --- Integration loop over integration points ----------------------------

    for (int i = 0; i < n1; i++)
    {
      // Compute basis function derivatives at current integration point
      if (nGauss < 1)
      {
	fe.N.fill(0.0);
	fe.N(i+1) = 1.0;
	dNdu.fillColumn(1,D1.getRow(i+1));
      }
      else
	if (!Lagrange::computeBasis(fe.N,&dNdu,points1,xg1[i]))
	  return false;

      // Compute Jacobian inverse of coordinate mapping and derivatives
      fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

      // Cartesian coordinates of current integration point
      X = Xnod*fe.N;
      X.t = time.t;

      // Evaluate the integrand and accumulate element contributions
      fe.detJxW *= wg1[i];
      if (!integrand.evalInt(*A,fe,time,X))
	return false;
    }

    // Assembly of global system integral
    if (!glInt.assemble(A->ref(),fe.iel))
      return false;

    A->destruct();
  }

  return true;
}


bool ASMs1DSpec::integrate (Integrand& integrand, int lIndex,
			    GlobalIntegral& glInt,
			    const TimeDomain& time)
{
  if (!curv) return true; // silently ignore empty patches

  return false; // not implemented
}


bool ASMs1DSpec::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			       const RealArray*, bool) const
{
  sField.resize(0,0);
  if (!curv) return false;

  const int p1 = curv->order();

  Matrix D1;
  if (!Legendre::basisDerivatives(p1,D1))
    return false;

  size_t nPoints = this->getNoNodes();
  IntVec check(nPoints,0);

  FiniteElement fe(p1);
  Vector        solPt;
  Vectors       globSolPt(nPoints);
  Matrix        dNdu(p1,1), Xnod, Jac;

  // Evaluate the secondary solution field at each point
  const int nel = this->getNoElms();
  for (int iel = 1; iel <= nel; iel++)
  {
    const IntVec& mnpc = MNPC[iel-1];
    this->getElementCoordinates(Xnod,iel);

    for (int i = 0; i < p1; i++)
    {
      fe.N.fill(0.0);
      fe.N(i+1) = 1.0;
      dNdu.fillColumn(1,D1.getRow(i+1));

      // Compute the Jacobian inverse
      fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

      // Now evaluate the solution field
      if (!integrand.evalSol(solPt,fe,Xnod.getColumn(i+1),mnpc))
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
