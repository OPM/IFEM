// $Id: ASMs3DmxLag.C,v 1.2 2010-12-30 15:02:02 kmo Exp $
//==============================================================================
//!
//! \file ASMs3DmxLag.C
//!
//! \date Dec 28 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of structured 3D Lagrange mixed FE models.
//!
//==============================================================================

#include "GoTools/trivariate/SplineVolume.h"

#include "ASMs3DmxLag.h"
#include "Lagrange.h"
#include "TimeDomain.h"
#include "GlobalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "Utilities.h"
#include "Vec3Oper.h"


ASMs3DmxLag::ASMs3DmxLag (const char* fileName, bool checkRHS,
			  unsigned char n_f1, unsigned char n_f2)
  : ASMs3DLag(fileName,checkRHS), ASMmxBase(n_f1,n_f2)
{
  nx2 = ny2 = nz2 = 0;
  nf = nf1 + nf2;
}


ASMs3DmxLag::ASMs3DmxLag (std::istream& is, bool checkRHS,
			  unsigned char n_f1, unsigned char n_f2)
  : ASMs3DLag(is,checkRHS), ASMmxBase(n_f1,n_f2)
{
  nx2 = ny2 = nz2 = 0;
  nf = nf1 + nf2;
}


ASMs3DmxLag::ASMs3DmxLag (bool checkRHS,
			  unsigned char n_f1, unsigned char n_f2)
  : ASMs3DLag(checkRHS), ASMmxBase(n_f1,n_f2)
{
  nx2 = ny2 = nz2 = 0;
  nf = nf1 + nf2;
}


void ASMs3DmxLag::clear ()
{
  nx2 = ny2 = nz2 = 0;
  ASMs3DLag::clear();
}


unsigned char ASMs3DmxLag::getNoFields (int basis) const
{
  switch (basis)
    {
    case 1: return nf1;
    case 2: return nf2;
    }

  return nf;
}


unsigned char ASMs3DmxLag::getNodalDOFs (size_t inod) const
{
  return inod <= nb1 ? nf1 : nf2;
}


void ASMs3DmxLag::initMADOF (const int* sysMadof)
{
  this->init(MLGN,sysMadof);
}


void ASMs3DmxLag::extractNodeVec (const Vector& globRes, Vector& nodeVec,
				  unsigned char) const
{
  this->extrNodeVec(globRes,nodeVec);
}


bool ASMs3DmxLag::generateFEMTopology ()
{
  // Generate/check FE data for the geometry/basis1
  bool haveFEdata = !MLGN.empty();
  bool basis1IsOK = this->ASMs3DLag::generateFEMTopology();
  if (haveFEdata || !basis1IsOK) return basis1IsOK;

  // Order of 2nd basis in the two parametric directions (order = degree + 1)
  const int p1 = svol->order(0)-1;
  const int p2 = svol->order(1)-1;
  const int p3 = svol->order(2)-1;
  if (p1 < 2 || p2 < 2 || p3 < 2)
  {
    std::cerr <<" *** ASMs3DmxLag::generateFEMTopology: Too low order "<< p1
	      <<","<< p2 <<","<< p3 <<" for the second basis."<< std::endl;
    return false;
  }

  // Evaluate the parametric values
  RealArray gpar[3];
  if (!this->getGridParameters(gpar[0],0,p1-1)) return false;
  if (!this->getGridParameters(gpar[1],1,p2-1)) return false;
  if (!this->getGridParameters(gpar[2],1,p3-1)) return false;

  // Number of nodes in each direction
  nx2 = gpar[0].size();
  ny2 = gpar[1].size();
  nz2 = gpar[2].size();

  // Add nodes for second basis (coordinates are not needed)
  nb1 = MLGN.size();
  nb2 = nx2*ny2*nz2;
  MLGN.reserve(nb1+nb2);
  int i, j, k, a, b, c, iel;
  for (k = 0; k < nz2; k++)
    for (j = 0; j < ny2; j++)
      for (i = 0; i < nx2; i++)
	MLGN.push_back(++gNod);

  // Number of elements in each direction
  const int nelx = (nx2-1)/(p1-1);
  const int nely = (ny2-1)/(p2-1);
  const int nelz = (nz2-1)/(p3-1);

  // Add connectivety for second basis: local --> global node relation
  for (k = iel = 0; k < nelz; k++)
    for (j = 0; j < nely; j++)
      for (i = 0; i < nelx; i++, iel++)
      {
	size_t nen1 = MNPC[iel].size();
	MNPC[iel].resize(nen1+p1*p2*p3);

	// First node in current element
	int corner = nb1 + (p3-1)*nx2*ny2*k + (p2-1)*nx2*j + (p1-1)*i;

	for (c = 0; c < p3; c++)
	{
	  int cornod = nen1 + p1*p2*c;
	  MNPC[iel][cornod] = corner + nx2*ny2*c;
	  for (b = 1; b < p2; b++)
	  {
	    int facenod = cornod + p1*b;
	    MNPC[iel][facenod] = MNPC[iel][cornod] + nx2*b;
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


bool ASMs3DmxLag::getSize (int& n1, int& n2, int& n3, int basis) const
{
  if (basis != 2)
    return this->ASMs3DLag::getSize(n1,n2,n3,1);

  n1 = nx2;
  n2 = ny2;
  n3 = nz2;

  return true;
}


bool ASMs3DmxLag::integrate (Integrand& integrand,
			     GlobalIntegral& glInt,
			     const TimeDomain& time,
			     const LintegralVec& locInt)
{
  if (!svol) return true; // silently ignore empty patches

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Order of basis in the two parametric directions (order = degree + 1)
  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);
  // The second basis is assumed one order lower in each direction
  const int q1 = p1 - 1;
  const int q2 = p2 - 1;
  const int q3 = p3 - 1;

  Vector N1(p1*p2*p3), N2(q1*q2*q3);
  Matrix dN1du, dN2du, dN1dX, dN2dX, Xnod, Jac;
  Vec4   X;


  // === Assembly loop over all elements in the patch ==========================

  for (int iel = 1; iel <= this->getNoElms(); iel++)
    {
        // Set up control point coordinates for current element
        if (!this->getElementCoordinates(Xnod,iel)) return false;

	// Initialize element quantities
	IntVec::iterator f2start = MNPC[iel-1].begin() + N1.size();
	if (!integrand.initElement(IntVec(MNPC[iel-1].begin(),f2start),
				   IntVec(f2start,MNPC[iel-1].end()),nb1))
	  return false;

	// Caution: Unless locInt is empty, we assume it points to an array of
	// LocalIntegral pointers, of length at least the number of elements in
	// the model (as defined by the highest number in the MLGE array).
	// If the array is shorter than this, expect a segmentation fault.
	LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[MLGE[iel-1]-1];


	// --- Integration loop over all Gauss points in each direction --------

	for (int k = 0; k < nGauss; k++)
	  for (int j = 0; j < nGauss; j++)
	    for (int i = 0; i < nGauss; i++)
	    {
	      // Weight of current integration point
	      double weight = wg[i]*wg[j]*wg[k];

	      // Compute basis function derivatives at current integration point
	      // using tensor product of one-dimensional Lagrange polynomials
	      if (!Lagrange::computeBasis(N1,dN1du,p1,xg[i],p2,xg[j],p3,xg[k]))
		return false;
	      if (!Lagrange::computeBasis(N2,dN2du,q1,xg[i],q2,xg[j],q3,xg[k]))
		return false;

	      // Compute Jacobian inverse of coordinate mapping and derivatives
	      double dA = utl::Jacobian(Jac,dN1dX,Xnod,dN1du);
	      if (dA == 0.0) continue; // skip singular points

	      dN2dX.multiply(dN2du,Jac); // dN2dX = dN2du * J^-1

	      // Cartesian coordinates of current integration point
	      X = Xnod * N1;
	      X.t = time.t;

	      // Evaluate the integrand and accumulate element contributions
	      if (!integrand.evalInt(elmInt,time,dA*weight,N1,N2,dN1dX,dN2dX,X))
		return false;
	    }

	// Assembly of global system integral
	if (!glInt.assemble(elmInt,MLGE[iel-1]))
	  return false;
    }

  return true;
}


bool ASMs3DmxLag::integrate (Integrand& integrand, int lIndex,
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
  // The second basis is assumed one order lower in each direction
  const int q1 = p1 - 1;
  const int q2 = p2 - 1;
  const int q3 = p3 - 1;

  // Number of elements in each direction
  const int nelx = (nx2-1)/(q1-1);
  const int nely = (ny2-1)/(q2-1);
  const int nelz = (nz2-1)/(q3-1);

  Vector N1(p1*p2*p3), N2(q1*q2*q3);
  Matrix dN1du, dN2du, dN1dX, dN2dX, Xnod, Jac;
  Vec4   X;
  Vec3   normal;
  double xi[3];


  // === Assembly loop over all elements on the patch face =====================

  int iel = 1;
  for (int i3 = 0; i3 < nelz; i3++)
    for (int i2 = 0; i2 < nely; i2++)
      for (int i1 = 0; i1 < nelx; i1++, iel++)
      {
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

	// Set up control point coordinates for current element
	if (!this->getElementCoordinates(Xnod,iel)) return false;

	// Initialize element quantities
	IntVec::iterator f2start = MNPC[iel-1].begin() + N1.size();
	if (!integrand.initElementBou(IntVec(MNPC[iel-1].begin(),f2start),
				      IntVec(f2start,MNPC[iel-1].end()),nb1))
	  return false;

	// Caution: Unless locInt is empty, we assume it points to an array of
	// LocalIntegral pointers, of length at least the number of elements in
	// the model (as defined by the highest number in the MLGE array).
	// If the array is shorter than this, expect a segmentation fault.
	LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[MLGE[iel-1]-1];


	// --- Integration loop over all Gauss points in each direction --------

	for (int j = 0; j < nGauss; j++)
	  for (int i = 0; i < nGauss; i++)
	  {
	    // Weight of current integration point
	    double weight = wg[i]*wg[j];

	    // Gauss point coordinates on the face
	    xi[t0-1] = faceDir < 0 ? -1.0 : 1.0;
	    xi[t1-1] = xg[i];
	    xi[t2-1] = xg[j];

	    // Compute the basis functions and their derivatives, using
	    // tensor product of one-dimensional Lagrange polynomials
	    if (!Lagrange::computeBasis(N1,dN1du,p1,xi[0],p2,xi[1],p3,xi[2]))
	      return false;
	    if (!Lagrange::computeBasis(N2,dN2du,q1,xi[0],q2,xi[1],q3,xi[2]))
	      return false;

	    // Compute basis function derivatives and the edge normal
	    double dS = utl::Jacobian(Jac,normal,dN1dX,Xnod,dN1du,t1,t2);
	    if (dS == 0.0) continue; // skip singular points

	    dN2dX.multiply(dN2du,Jac); // dN2dX = dN2du * J^-1

	    if (faceDir < 0) normal *= -1.0;

	    // Cartesian coordinates of current integration point
	    X = Xnod * N1;
	    X.t = time.t;

	    // Evaluate the integrand and accumulate element contributions
	    if (!integrand.evalBou(elmInt,time,dS*weight,
				   N1,N2,dN1dX,dN2dX,X,normal))
	      return false;
	  }

	// Assembly of global system integral
	if (!glInt.assemble(elmInt,MLGE[iel-1]))
	  return false;
      }

  return true;
}


bool ASMs3DmxLag::evalSolution (Matrix& sField, const Vector& locSol,
				const int*) const
{
  size_t nc1 = nf1;
  size_t nc2 = 0;
  if (nc1*nb1 < locSol.size())
    nc2 = (locSol.size() - nc1*nb1)/nb2;
  else
    nc1 = locSol.size()/nb1;

  if (nc1*nb1 + nc2*nb2 != locSol.size())
    return false;

  // TODO: Add evaluation second field at the nodes of the first field
  size_t nPoints = nb1;
  size_t nComp = nc1;
  size_t i, n, ip = 0;
  sField.resize(nComp,nPoints);
  for (n = 1; n <= nPoints; n++)
    for (i = 1; i <= nComp; i++)
      sField(i,n) = locSol(++ip);

  return true;
}


bool ASMs3DmxLag::evalSolution (Matrix& sField, const Integrand& integrand,
				const int*) const
{
  sField.resize(0,0);

  if (!svol) return false;

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);
  const int q1 = p1 - 1;
  const int q2 = p2 - 1;
  const int q3 = p3 - 1;
  double incx = 2.0/double(q1);
  double incy = 2.0/double(q2);
  double incz = 2.0/double(q3);

  size_t nPoints = nb1;
  IntVec check(nPoints,0);

  Vector N1(p1*p2*p3), N2(q1*q2*q3), solPt;
  std::vector<Vector> globSolPt(nPoints);
  Matrix dN1du, dN2du, dN1dX, dN2dX, Xnod, Jac;

  // Evaluate the secondary solution field at each point
  for (int iel = 1; iel <= this->getNoElms(); iel++)
  {
    IntVec::const_iterator f2start = MNPC[iel-1].begin() + p1*p2*p3;
    IntVec mnpc1(MNPC[iel-1].begin(),f2start);
    IntVec mnpc2(f2start,MNPC[iel-1].end());

    this->getElementCoordinates(Xnod,iel);

    int i, j, k, loc = 0;
    for (k = 0; k < p3; k++)
      for (j = 0; j < p2; j++)
	for (i = 0; i < p1; i++, loc++)
	{
	  double xi   = -1.0 + i*incx;
	  double eta  = -1.0 + j*incy;
	  double zeta = -1.0 + k*incz;
	  if (!Lagrange::computeBasis(N1,dN1du,p1,xi,p2,eta,p3,zeta))
	    return false;
	  if (!Lagrange::computeBasis(N2,dN2du,q1,xi,q2,eta,q3,zeta))
	    return false;

	  // Compute the Jacobian inverse
	  if (utl::Jacobian(Jac,dN1dX,Xnod,dN1du) == 0.0) // Jac=(Xnod*dN1du)^-1
	    continue; // skip singular points
	  else
	    dN2dX.multiply(dN2du,Jac); // dN2dX = dN2du * J^-1

	  // Now evaluate the solution field
	  if (!integrand.evalSol(solPt,N1,N2,dN1dX,dN2dX,Xnod*N1,mnpc1,mnpc2))
	    return false;
	  else if (sField.empty())
	    sField.resize(solPt.size(),nPoints,true);

	  if (++check[mnpc1[loc]] == 1)
	    globSolPt[mnpc1[loc]] = solPt;
	  else
	    globSolPt[mnpc1[loc]] += solPt;
	}
  }

  for (size_t i = 0; i < nPoints; i++)
    sField.fillColumn(1+i,globSolPt[i]/=check[i]);

  return true;
}
