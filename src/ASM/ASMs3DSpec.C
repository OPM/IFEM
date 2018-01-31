// $Id$
//==============================================================================
//!
//! \file ASMs3DSpec.C
//!
//! \date Mar 22 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for assembly of structured 3D Spectral FE models.
//!
//==============================================================================

#include "GoTools/trivariate/SplineVolume.h"

#include "ASMs3DSpec.h"
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "Vec3Oper.h"
#include "Legendre.h"
#include <array>


bool ASMs3DSpec::getGridParameters (RealArray& prm, int dir,
				    int nSegPerSpan) const
{
  if (!svol) return false;

  // Evaluate the Gauss-Lobatto-Legendre points in this direction
  Vector dummy, xGLL;
  if (!Legendre::GLL(dummy,xGLL,svol->order(dir))) return false;

  if (xGLL.size() != (size_t)(nSegPerSpan+1))
  {
    nSegPerSpan = xGLL.size() - 1;
    std::cout <<"Spectral elements: Number of nodes per knot-span in "
	      << char('u'+dir) <<"-directon reset to "<< nSegPerSpan
	      <<" (GLL points)"<< std::endl;
  }

  RealArray::const_iterator uit = svol->basis(dir).begin();
  double ucurr, uprev = *(uit++);
  while (uit != svol->basis(dir).end())
  {
    ucurr = *(uit++);
    if (ucurr > uprev)
      for (int i = 1; i <= nSegPerSpan; i++)
	prm.push_back(0.5*(ucurr-uprev)*(1.0+xGLL(i)) + uprev);
    uprev = ucurr;
  }

  prm.push_back(svol->basis(dir).endparam());
  return true;
}


/*!
  \brief Establishes matrices with basis functions and 1st derivatives.
*/

static void evalBasis (int i, int j, int k, int p1, int p2, int p3,
		       const Matrix& der1, const Matrix& der2,
		       const Matrix& der3, Vector& N, Matrix& dNdu)
{
  int a, b, c, n = 1;
  for (c = 1; c <= p3; c++)
    for (b = 1; b <= p2; b++)
      for (a = 1; a <= p1; a++, n++)
      {
	 N  (n)   = (a == i && b == j && c == k) ? 1.0 : 0.0;
	dNdu(n,1) = (b == j && c == k) ? der1(i,a) : 0.0;
	dNdu(n,2) = (a == i && c == k) ? der2(j,b) : 0.0;
	dNdu(n,3) = (a == i && b == j) ? der3(k,c) : 0.0;
      }
}


bool ASMs3DSpec::integrate (Integrand& integrand,
			    GlobalIntegral& glInt,
			    const TimeDomain& time)
{
  if (this->empty()) return true; // silently ignore empty patches

  // Evaluate integration points (=nodal points) and weights

  Vector wg1,xg1,wg2,xg2,wg3,xg3;
  if (!Legendre::GLL(wg1,xg1,p1)) return false;
  if (!Legendre::GLL(wg2,xg2,p2)) return false;
  if (!Legendre::GLL(wg3,xg3,p3)) return false;

  Matrix D1, D2, D3;
  if (!Legendre::basisDerivatives(p1,D1)) return false;
  if (!Legendre::basisDerivatives(p2,D2)) return false;
  if (!Legendre::basisDerivatives(p3,D3)) return false;


  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t g = 0; g < threadGroupsVol.size() && ok; g++)
  {
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < threadGroupsVol[g].size(); t++)
    {
      FiniteElement fe(p1*p2*p3);
      Matrix   dNdu(p1*p2*p3,3), Xnod, Jac;
      Vec4     X;
      for (size_t l = 0; l < threadGroupsVol[g][t].size(); l++)
      {
        int iel = threadGroupsVol[g][t][l]+1;

        // Set up nodal point coordinates for current element
        if (!this->getElementCoordinates(Xnod,iel))
        {
          ok = false;
          break;
        }

	// Initialize element quantities
        fe.iel = MLGE[iel-1];
        LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel);
        if (!integrand.initElement(MNPC[iel-1],*A))
        {
          A->destruct();
          ok = false;
          break;
        }


        // --- Integration loop over all Gauss points in each direction --------

        int count = 1;
        for (int k = 1; k <= p3; k++)
          for (int j = 1; j <= p2; j++)
            for (int i = 1; i <= p1; i++, count++)
            {
              // Evaluate the basis functions and gradients using tensor product
              // of the one-dimensional Lagrange polynomials
              evalBasis(i,j,k,p1,p2,p3,D1,D2,D3,fe.N,dNdu);

              // Compute Jacobian inverse of coordinate mapping and derivatives
              fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);
              if (fe.detJxW == 0.0) continue; // skip singular points

              // Cartesian coordinates of current integration point
              X.x = Xnod(1,count);
              X.y = Xnod(2,count);
              X.z = Xnod(3,count);
              X.t = time.t;

              // Evaluate the integrand and accumulate element contributions
              fe.detJxW *= wg1(i)*wg2(j)*wg3(k);
              if (!integrand.evalInt(*A,fe,time,X))
                ok = false;
            }

	// Assembly of global system integral
	if (ok && !glInt.assemble(A->ref(),fe.iel))
          ok = false;

        A->destruct();
      }
    }
  }

  return ok;
}


bool ASMs3DSpec::integrate (Integrand& integrand, int lIndex,
			    GlobalIntegral& glInt,
			    const TimeDomain& time)
{
  if (this->empty()) return true; // silently ignore empty patches

  std::map<char,ThreadGroups>::const_iterator tit;
  if ((tit = threadGroupsFace.find(lIndex)) == threadGroupsFace.end())
  {
    std::cerr <<" *** ASMs3DSpec::integrate: No thread groups for face "<<lIndex
	      << std::endl;
    return false;
  }
  const ThreadGroups& threadGrp = tit->second;

  // Find the parametric direction of the face normal {-3,-2,-1, 1, 2, 3}
  const int faceDir = (lIndex+1)/(lIndex%2 ? -2 : 2);

  const int t0 = abs(faceDir); // unsigned normal direction of the face
  const int t1 = 1 + t0%3; // first tangent direction of the face
  const int t2 = 1 + t1%3; // second tangent direction of the face

  // Evaluate integration points (=nodal points) and weights

  std::array<Vector,3> xg, wg;
  std::array<Matrix,3> D;
  std::array<int,3> p({{p1,p2,p3}});
  for (int d = 0; d < 3; d++)
  {
    if (!Legendre::GLL(wg[d],xg[d],p[d])) return false;
    if (!Legendre::basisDerivatives(p[d],D[d])) return false;
  }
  int nen = p1*p2*p3;


  // === Assembly loop over all elements on the patch face =====================

  bool ok = true;
  for (size_t g = 0; g < threadGrp.size() && ok; g++)
  {
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < threadGrp[g].size(); t++)
    {
      FiniteElement fe(nen);
      Matrix dNdu(nen,3), Xnod, Jac;
      Vec4   X;
      Vec3   normal;
      int    xi[3];
      for (size_t l = 0; l < threadGrp[g][t].size(); l++)
      {
        int iel = threadGrp[g][t][l];

	// Set up nodal point coordinates for current element
        if (!this->getElementCoordinates(Xnod,++iel))
        {
          ok = false;
          break;
        }

	// Initialize element quantities
        fe.iel = MLGE[iel-1];
        LocalIntegral* A = integrand.getLocalIntegral(nen,fe.iel,true);
        if (!integrand.initElementBou(MNPC[iel-1],*A))
        {
          A->destruct();
          ok = false;
          break;
        }


	// --- Integration loop over all Gauss points in each direction --------

	for (int j = 0; j < p[t2-1]; j++)
	  for (int i = 0; i < p[t1-1]; i++)
	  {
	    // "Coordinates" on the face
	    xi[t0-1] = faceDir < 0 ? 1 : p[t0-1];
	    xi[t1-1] = i+1;
	    xi[t2-1] = j+1;

	    // Compute the basis functions and their derivatives, using
	    // tensor product of one-dimensional Lagrange polynomials
	    evalBasis(xi[0],xi[1],xi[2],p1,p2,p3,D[0],D[1],D[2],fe.N,dNdu);

	    // Compute basis function derivatives and the face normal
	    fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,Xnod,dNdu,t1,t2);
	    if (fe.detJxW == 0.0) continue; // skip singular points

	    if (faceDir < 0) normal *= -1.0;

	    // Cartesian coordinates of current integration point
	    X = Xnod * fe.N;
	    X.t = time.t;

	    // Evaluate the integrand and accumulate element contributions
	    fe.detJxW *= wg[t1-1][i]*wg[t2-1][j];
            if (!integrand.evalBou(*A,fe,time,X,normal))
              ok = false;
	  }

        // Finalize the element quantities
        if (ok && !integrand.finalizeElementBou(*A,fe,time))
          ok = false;

	// Assembly of global system integral
        if (ok && !glInt.assemble(A->ref(),fe.iel))
          ok = false;

        A->destruct();
      }
    }
  }

  return ok;
}


bool ASMs3DSpec::integrateEdge (Integrand& integrand, int lEdge,
				GlobalIntegral& glInt,
				const TimeDomain& time)
{
  if (this->empty()) return true; // silently ignore empty patches

  // Parametric direction of the edge {0, 1, 2}
  const int lDir = (lEdge-1)/4;
  // Order of basis in the edge directions (order = degree + 1)
  const int pe = lDir == 0 ? p1 : (lDir == 1 ? p2 : p3);

  // Number of elements in each direction
  int n1, n2, n3;
  this->getSize(n1,n2,n3);
  const int nelx = (n1-1)/(p1-1);
  const int nely = (n2-1)/(p2-1);
  const int nelz = (n3-1)/(p3-1);

  // Evaluate integration points (=nodal points) and weights

  std::array<Vector,3> wg, xg;
  if (!Legendre::GLL(wg[0],xg[0],p1)) return false;
  if (!Legendre::GLL(wg[1],xg[1],p2)) return false;
  if (!Legendre::GLL(wg[2],xg[2],p3)) return false;

  Matrix D1, D2, D3;
  if (!Legendre::basisDerivatives(p1,D1)) return false;
  if (!Legendre::basisDerivatives(p2,D2)) return false;
  if (!Legendre::basisDerivatives(p3,D3)) return false;

  const int nen = p1*p2*p3;

  FiniteElement fe(nen);
  Matrix dNdu(nen,3), Xnod, Jac;
  Vec4   X;
  Vec3   tangent;
  int    xi[3];

  switch (lEdge)
    {
    case  1: xi[1] =  1; xi[2] =  1; break;
    case  2: xi[1] = p2; xi[2] =  1; break;
    case  3: xi[1] =  1; xi[2] = p3; break;
    case  4: xi[1] = p2; xi[2] = p3; break;
    case  5: xi[0] =  1; xi[2] =  1; break;
    case  6: xi[0] = p1; xi[2] =  1; break;
    case  7: xi[0] =  1; xi[2] = p3; break;
    case  8: xi[0] = p1; xi[2] = p3; break;
    case  9: xi[0] =  1; xi[1] =  1; break;
    case 10: xi[0] = p1; xi[1] =  1; break;
    case 11: xi[0] =  1; xi[1] = p2; break;
    case 12: xi[0] = p1; xi[1] = p2; break;
    }


  // === Assembly loop over all elements on the patch edge =====================

  int iel = 1;
  for (int i3 = 0; i3 < nelz; i3++)
    for (int i2 = 0; i2 < nely; i2++)
      for (int i1 = 0; i1 < nelx; i1++, iel++)
      {
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
	fe.iel = MLGE[iel-1];
        LocalIntegral* A = integrand.getLocalIntegral(nen,fe.iel,true);
        if (!integrand.initElementBou(MNPC[iel-1],*A)) return false;


	// --- Integration loop over all Gauss points along the edge -----------

	for (int i = 0; i < pe; i++)
	{
	  // "Coordinate" on the edge
	  xi[lDir] = i+1;

	  // Compute the basis functions and their derivatives, using
	  // tensor product of one-dimensional Lagrange polynomials
	  evalBasis(xi[0],xi[1],xi[2],p1,p2,p3,D1,D2,D3,fe.N,dNdu);

	  // Compute basis function derivatives and the edge tangent
	  fe.detJxW = utl::Jacobian(Jac,tangent,fe.dNdX,Xnod,dNdu,1+lDir);
	  if (fe.detJxW == 0.0) continue; // skip singular points

	  // Cartesian coordinates of current integration point
	  X = Xnod * fe.N;
	  X.t = time.t;

	  // Evaluate the integrand and accumulate element contributions
	  fe.detJxW *= wg[lDir][i];
          if (!integrand.evalBou(*A,fe,time,X,tangent))
	    return false;
	}

        // Finalize the element quantities
        if (!integrand.finalizeElementBou(*A,fe,time))
          return false;

	// Assembly of global system integral
	if (!glInt.assemble(A->ref(),fe.iel))
	  return false;
      }

  return true;
}


bool ASMs3DSpec::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			       const RealArray*, bool) const
{
  sField.resize(0,0);

  Vector wg1,xg1,wg2,xg2,wg3,xg3;
  if (!Legendre::GLL(wg1,xg1,p1)) return false;
  if (!Legendre::GLL(wg2,xg2,p2)) return false;
  if (!Legendre::GLL(wg3,xg3,p3)) return false;

  Matrix D1, D2, D3;
  if (!Legendre::basisDerivatives(p1,D1)) return false;
  if (!Legendre::basisDerivatives(p2,D2)) return false;
  if (!Legendre::basisDerivatives(p3,D3)) return false;

  size_t nPoints = this->getNoNodes();
  IntVec check(nPoints,0);

  FiniteElement fe(p1*p2*p3);
  Vector        solPt;
  Vectors       globSolPt(nPoints);
  Matrix        dNdu(p1*p2*p3,3), Xnod, Jac;

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
	  evalBasis(i+1,j+1,k+1,p1,p2,p3,D1,D2,D3,fe.N,dNdu);

	  // Compute the Jacobian inverse
	  fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

	  // Now evaluate the solution field
	  if (!integrand.evalSol(solPt,fe,Xnod.getColumn(loc+1),mnpc))
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
