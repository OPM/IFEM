// $Id$
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
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "Utilities.h"
#include "Vec3Oper.h"


ASMs3DmxLag::ASMs3DmxLag (unsigned char n_f1, unsigned char n_f2)
  : ASMs3DLag(n_f1+n_f2), ASMmxBase(n_f1,n_f2)
{
  nx2 = ny2 = nz2 = 0;
}


ASMs3DmxLag::ASMs3DmxLag (const ASMs3DmxLag& patch, char n_f1, char n_f2)
  : ASMs3DLag(patch), ASMmxBase(patch.nf1,patch.nf2)
{
  nx2 = patch.nx2;
  ny2 = patch.ny2;
  nz2 = patch.nz2;
  nb1 = patch.nb1;
  nb2 = patch.nb2;
  if (n_f1 >= 0) nf1 = n_f1;
  if (n_f2 >= 0) nf2 = n_f2;
  nf = nf1 + nf2;
}


void ASMs3DmxLag::clear (bool retainGeometry)
{
  nx2 = ny2 = nz2 = 0;
  this->ASMs3DLag::clear(retainGeometry);
}


size_t ASMs3DmxLag::getNoNodes (int basis) const
{
  switch (basis)
    {
    case 1: return nb1;
    case 2: return nb2;
    }

  return MLGN.size();
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
  return inod <= nb1 || inod > nb1+nb2 ? nf1 : nf2;
}


char ASMs3DmxLag::getNodeType (size_t inod) const
{
  return inod <= nb1 ? 'D' : (inod <= nb1+nb2 ? 'P' : 'X');
}


void ASMs3DmxLag::initMADOF (const int* sysMadof)
{
  this->initMx(MLGN,sysMadof);
}


void ASMs3DmxLag::extractNodeVec (const Vector& globRes, Vector& nodeVec,
				  unsigned char, int basis) const
{
  this->extractNodeVecMx(globRes,nodeVec,basis);
}


bool ASMs3DmxLag::getSolution (Matrix& sField, const Vector& locSol,
			       const IntVec& nodes) const
{
  return this->getSolutionMx(sField,locSol,nodes);
}


bool ASMs3DmxLag::generateFEMTopology ()
{
  // Generate/check FE data for the geometry/basis1
  bool haveFEdata = !MLGN.empty();
  bool basis1IsOK = this->ASMs3DLag::generateFEMTopology();
  if ((haveFEdata && !shareFE) || !basis1IsOK) return basis1IsOK;

  // Order of 2nd basis in the three parametric directions (order = degree + 1)
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
  if (!this->getGridParameters(gpar[2],2,p3-1)) return false;

  // Number of nodes in each direction
  nx2 = gpar[0].size();
  ny2 = gpar[1].size();
  nz2 = gpar[2].size();

  nb1 = MLGN.size();
  nb2 = nx2*ny2*nz2;

  if (shareFE) return true;

  // Add nodes for second basis (coordinates are not needed)
  myMLGN.reserve(nb1+nb2);
  for (size_t i3 = 0; i3 < nz2; i3++)
    for (size_t i2 = 0; i2 < ny2; i2++)
      for (size_t i1 = 0; i1 < nx2; i1++)
	myMLGN.push_back(++gNod);

  // Number of elements in each direction
  const int nelx = (nx2-1)/(p1-1);
  const int nely = (ny2-1)/(p2-1);
  const int nelz = (nz2-1)/(p3-1);

  // Add connectivity for second basis: local --> global node relation
  int i, j, k, a, b, c, iel;
  for (k = iel = 0; k < nelz; k++)
    for (j = 0; j < nely; j++)
      for (i = 0; i < nelx; i++, iel++)
      {
	size_t nen1 = myMNPC[iel].size();
	myMNPC[iel].resize(nen1+p1*p2*p3);

	// First node in current element
	int corner = nb1 + (p3-1)*nx2*ny2*k + (p2-1)*nx2*j + (p1-1)*i;

	for (c = 0; c < p3; c++)
	{
	  int cornod = nen1 + p1*p2*c;
	  myMNPC[iel][cornod] = corner + nx2*ny2*c;
	  for (b = 1; b < p2; b++)
	  {
	    int facenod = cornod + p1*b;
	    myMNPC[iel][facenod] = myMNPC[iel][cornod] + nx2*b;
	    for (a = 1; a < p1; a++)
	    {
	      myMNPC[iel][facenod+a] = myMNPC[iel][facenod] + a;
	      myMNPC[iel][cornod+a]  = myMNPC[iel][cornod] + a;
	    }
	  }
	}
      }

  return true;
}


bool ASMs3DmxLag::connectPatch (int face, ASMs3D& neighbor,
				int nface, int norient)
{
  ASMs3DmxLag* neighMx = dynamic_cast<ASMs3DmxLag*>(&neighbor);
  if (!neighMx) return false;

  if (swapW && face > 4) // Account for swapped parameter direction
    face = 11-face;

  if (neighMx->swapW && face > 4) // Account for swapped parameter direction
    nface = 11-nface;

  return this->connectBasis(face,neighbor,nface,norient,1,0,0)
    &&   this->connectBasis(face,neighbor,nface,norient,2,nb1,neighMx->nb1);
}


void ASMs3DmxLag::closeFaces (int dir, int, int)
{
  this->ASMs3D::closeFaces(dir,1,1);
  this->ASMs3D::closeFaces(dir,2,nb1+1);
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
			     const TimeDomain& time)
{
  if (!svol) return true; // silently ignore empty patches

  // Get Gaussian quadrature points and weights
  const double* x = GaussQuadrature::getCoord(nGauss);
  const double* w = GaussQuadrature::getWeight(nGauss);
  if (!x || !w) return false;

  // Get parametric coordinates of the elements
  RealArray upar, vpar, wpar;
  this->getGridParameters(upar,0,1);
  this->getGridParameters(vpar,1,1);
  this->getGridParameters(wpar,2,1);

  // Number of elements in each direction
  const int nelx = upar.size() - 1;
  const int nely = vpar.size() - 1;

  // Order of basis in the two parametric directions (order = degree + 1)
  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);
  // The second basis is assumed one order lower in each direction
  const int q1 = p1 - 1;
  const int q2 = p2 - 1;
  const int q3 = p3 - 1;


  // === Assembly loop over all elements in the patch ==========================

  bool ok=true;
  for (size_t g=0;g<threadGroupsVol.size() && ok;++g) {
#pragma omp parallel for schedule(static)
    for (size_t t=0;t<threadGroupsVol[g].size();++t) {
      MxFiniteElement fe(p1*p2*p3,q1*q2*q3);
      Matrix dN1du, dN2du, Xnod, Jac;
      Vec4   X;
      for (size_t l = 0; l < threadGroupsVol[g][t].size() && ok; ++l)
      {
        int iel = threadGroupsVol[g][t][l];
        int i1  =  iel % nelx;
        int i2  = (iel / nelx) % nely;
        int i3  =  iel / (nelx*nely);

        // Set up control point coordinates for current element
        if (!this->getElementCoordinates(Xnod,++iel))
	{
          ok = false;
          break;
        }

        // Initialize element quantities
        fe.iel = MLGE[iel-1];
        IntVec::const_iterator f2start = MNPC[iel-1].begin() + fe.N1.size();
        LocalIntegral* A = integrand.getLocalIntegral(fe.N1.size(),fe.N2.size(),
                                                      fe.iel,false);
        if (!integrand.initElement(IntVec(MNPC[iel-1].begin(),f2start),
                                   IntVec(f2start,MNPC[iel-1].end()),nb1,*A))
        {
          A->destruct();
          ok = false;
          break;
        }


        // --- Integration loop over all Gauss points in each direction --------

        int jp = ((i3*nely + i2)*nelx + i1)*nGauss*nGauss*nGauss;
        fe.iGP = firstIp + jp; // Global integration point counter

        for (int k = 0; k < nGauss; k++)
          for (int j = 0; j < nGauss; j++)
            for (int i = 0; i < nGauss; i++, fe.iGP++)
            {
              // Parameter value of current integration point
              fe.u = 0.5*(upar[i1]*(1.0-x[i]) + upar[i1+1]*(1.0+x[i]));
              fe.v = 0.5*(vpar[i2]*(1.0-x[j]) + vpar[i2+1]*(1.0+x[j]));
              fe.w = 0.5*(wpar[i3]*(1.0-x[k]) + wpar[i3+1]*(1.0+x[k]));

              // Local coordinate of current integration point
              fe.xi   = x[i];
              fe.eta  = x[j];
              fe.zeta = x[k];

              // Compute basis function derivatives at current integration point
              // using tensor product of one-dimensional Lagrange polynomials
              if (!Lagrange::computeBasis(fe.N1,dN1du,p1,x[i],p2,x[j],p3,x[k]) ||
                  !Lagrange::computeBasis(fe.N2,dN2du,q1,x[i],q2,x[j],q3,x[k]))
                ok = false;

              // Compute Jacobian inverse of coordinate mapping and derivatives
              fe.detJxW = utl::Jacobian(Jac,fe.dN1dX,Xnod,dN1du);
              if (fe.detJxW == 0.0) continue; // skip singular points

              fe.dN2dX.multiply(dN2du,Jac); // dN2dX = dN2du * J^-1

              // Cartesian coordinates of current integration point
              X = Xnod * fe.N1;
              X.t = time.t;

              // Evaluate the integrand and accumulate element contributions
              fe.detJxW *= w[i]*w[j]*w[k];
              if (!integrand.evalIntMx(*A,fe,time,X))
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


bool ASMs3DmxLag::integrate (Integrand& integrand, int lIndex,
			     GlobalIntegral& glInt,
			     const TimeDomain& time)
{
  if (!svol) return true; // silently ignore empty patches

  std::map<char,ThreadGroups>::const_iterator tit;
  if ((tit = threadGroupsFace.find(lIndex)) == threadGroupsFace.end())
  {
    std::cerr <<" *** ASMs3DLag::integrate: No thread groups for face "<< lIndex
	      << std::endl;
    return false;
  }
  const ThreadGroups& threadGrp = tit->second;

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
  const int nel1 = (nx2-1)/(p1-1);
  const int nel2 = (ny2-1)/(p2-1);

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;


  // === Assembly loop over all elements on the patch face =====================

  bool ok = true;
  for (size_t g = 0; g < threadGrp.size() && ok; ++g) {
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < threadGrp[g].size(); ++t) {
      MxFiniteElement fe(p1*p2*p3,q1*q2*q3);
      Matrix dN1du, dN2du, Xnod, Jac;
      Vec4   X;
      Vec3   normal;
      double xi[3];

      for (size_t l = 0; l < threadGrp[g][t].size() && ok; ++l) {
        int iel = threadGrp[g][t][l];
        int i1  =  iel % nel1;
        int i2  = (iel / nel1) % nel2;
        int i3  =  iel / (nel1*nel2);

	// Set up control point coordinates for current element
	if (!this->getElementCoordinates(Xnod,++iel))
        {
          ok = false;
          break;
        }

	// Initialize element quantities
	fe.iel = MLGE[iel-1];
	IntVec::const_iterator f2start = MNPC[iel-1].begin() + fe.N1.size();
        LocalIntegral* A = integrand.getLocalIntegral(fe.N1.size(),fe.N2.size(),
                                                      fe.iel,true);
	if (!integrand.initElementBou(IntVec(MNPC[iel-1].begin(),f2start),
				      IntVec(f2start,MNPC[iel-1].end()),nb1,*A))
        {
          A->destruct();
          ok = false;
          break;
        }

        // Define some loop control variables depending on which face we are on
        int nf1, j1, j2;
        switch (abs(faceDir))
        {
          case 1: nf1 = nel2; j2 = i3; j1 = i2; break;
          case 2: nf1 = nel1; j2 = i3; j1 = i1; break;
          case 3: nf1 = nel1; j2 = i2; j1 = i1; break;
          default: nf1 = j1 = j2 = 0;
        }


	// --- Integration loop over all Gauss points in each direction --------

        int jp = (j2*nf1 + j1)*nGauss*nGauss;
        fe.iGP = firstp + jp; // Global integration point counter

	for (int j = 0; j < nGauss; j++)
	  for (int i = 0; i < nGauss; i++, fe.iGP++)
	  {
	    // Gauss point coordinates on the face
	    xi[t0-1] = faceDir < 0 ? -1.0 : 1.0;
	    xi[t1-1] = xg[i];
	    xi[t2-1] = xg[j];

	    // Compute the basis functions and their derivatives, using
	    // tensor product of one-dimensional Lagrange polynomials
	    if (!Lagrange::computeBasis(fe.N1,dN1du,p1,xi[0],p2,xi[1],p3,xi[2]) ||
		!Lagrange::computeBasis(fe.N2,dN2du,q1,xi[0],q2,xi[1],q3,xi[2]))
              ok = false;

	    // Compute basis function derivatives and the edge normal
	    fe.detJxW = utl::Jacobian(Jac,normal,fe.dN1dX,Xnod,dN1du,t1,t2);
	    if (fe.detJxW == 0.0) continue; // skip singular points

	    fe.dN2dX.multiply(dN2du,Jac); // dN2dX = dN2du * J^-1

	    if (faceDir < 0) normal *= -1.0;

	    // Cartesian coordinates of current integration point
	    X = Xnod * fe.N1;
	    X.t = time.t;

	    // Evaluate the integrand and accumulate element contributions
	    fe.detJxW *= wg[i]*wg[j];
	    if (!integrand.evalBouMx(*A,fe,time,X,normal))
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


bool ASMs3DmxLag::evalSolution (Matrix& sField, const Vector& locSol,
				const RealArray*, bool) const
{
  size_t nc1 = nf1;
  size_t nc2 = 0;
  if (nc1*nb1 < locSol.size())
    nc2 = (locSol.size() - nc1*nb1)/nb2;
  else
    nc1 = locSol.size()/nb1;

  if (nc1*nb1 + nc2*nb2 != locSol.size())
    return false;

  // TODO: Add evaluation of the second field at the nodes of the first field
  size_t nPoints = nb1;
  size_t nComp = nc1;
  size_t i, n, ip = 0;
  sField.resize(nComp,nPoints);
  for (n = 1; n <= nPoints; n++)
    for (i = 1; i <= nComp; i++)
      sField(i,n) = locSol(++ip);

  return true;
}


bool ASMs3DmxLag::evalSolution (Matrix& sField, const IntegrandBase& integrand,
				const RealArray*, bool) const
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
  const int nel = this->getNoElms(true);
  for (int iel = 1; iel <= nel; iel++)
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
	  if (!Lagrange::computeBasis(N1,dN1du,p1,xi,p2,eta,p3,zeta) ||
	      !Lagrange::computeBasis(N2,dN2du,q1,xi,q2,eta,q3,zeta))
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
