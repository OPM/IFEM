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
#include <numeric>


ASMs3DmxLag::ASMs3DmxLag (const CharVec& n_f)
  : ASMs3DLag(std::accumulate(n_f.begin(),n_f.end(),0)), ASMmxBase(n_f)
{
}


ASMs3DmxLag::ASMs3DmxLag (const ASMs3DmxLag& patch, const CharVec& n_f)
  : ASMs3DLag(patch), ASMmxBase(n_f)
{
  nxx = patch.nxx;
  nyx = patch.nyx;
  nzx = patch.nzx;
}


void ASMs3DmxLag::clear (bool retainGeometry)
{
  nxx.clear();
  nyx.clear();
  nzx.clear();
  this->ASMs3DLag::clear(retainGeometry);
}


size_t ASMs3DmxLag::getNoNodes (int basis) const
{
  if (basis > (int)nb.size())
    basis = 1;

  if (basis == 0)
    return MLGN.size();

  return nb[basis-1];
}


unsigned char ASMs3DmxLag::getNoFields (int basis) const
{
  if (basis == 0)
    return std::accumulate(nfx.begin(), nfx.end(), 0);

  return nfx[basis-1];
}


unsigned char ASMs3DmxLag::getNodalDOFs (size_t inod) const
{
  if (this->isLMn(inod)) return nLag;
  size_t nbc=0;
  for (size_t i=0;i<nb.size();++i)
    if (inod <= (nbc+=nb[i]))
      return nfx[i];

  return nfx[0];
}


char ASMs3DmxLag::getNodeType (size_t inod) const
{
  if (this->isLMn(inod)) return 'L';
  size_t nbc=0;
  for (size_t i=0;i<nb.size();++i)
    if (inod <= (nbc+=nb[i]))
      return i == 0 ? 'D' : 'P';

  return 'X';
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
  geoBasis = 1;
  if ((haveFEdata && !shareFE) || !basis1IsOK) return basis1IsOK;

  // Order of 2nd basis in the three parametric directions (order = degree + 1)
  const size_t p1 = svol->order(0)-1;
  const size_t p2 = svol->order(1)-1;
  const size_t p3 = svol->order(2)-1;
  if (p1 < 2 || p2 < 2 || p3 < 2)
  {
    std::cerr <<" *** ASMs3DmxLag::generateFEMTopology: Too low order "<< p1
	      <<","<< p2 <<","<< p3 <<" for the second basis."<< std::endl;
    return false;
  }

  // Evaluate the parametric values
  RealArray gpar1, gpar2, gpar3;
  if (!this->getGridParameters(gpar1,0,p1-1)) return false;
  if (!this->getGridParameters(gpar2,1,p2-1)) return false;
  if (!this->getGridParameters(gpar3,2,p3-1)) return false;

  // Number of nodes in each direction for each basis
  nxx.resize(2);
  nyx.resize(2);
  nzx.resize(2);
  nxx[0] = nx;
  nyx[0] = ny;
  nzx[0] = nz;
  nxx[1] = gpar1.size();
  nyx[1] = gpar2.size();
  nzx[1] = gpar3.size();

  // Number of elements in each direction for each basis
  elem_sizes.resize(2);
  elem_sizes[0][0] = p1+1;
  elem_sizes[0][1] = p2+1;
  elem_sizes[0][2] = p3+1;
  elem_sizes[1][0] = p1;
  elem_sizes[1][1] = p2;
  elem_sizes[1][2] = p3;

  // Total number of nodes for each basis
  nb.resize(2);
  nb[0] = MLGN.size();
  nb[1] = nxx[1]*nyx[1]*nzx[1];

  if (shareFE == 'F') return true;

  // Add nodes for second basis (coordinates are not needed)
  nnod += nb[1];
  myMLGN.reserve(nnod);
  for (size_t i3 = 0; i3 < nzx[1]; i3++)
    for (size_t i2 = 0; i2 < nyx[1]; i2++)
      for (size_t i1 = 0; i1 < nxx[1]; i1++)
	myMLGN.push_back(++gNod);

  // Number of elements in each direction
  const int nelx = (nxx[1]-1)/(p1-1);
  const int nely = (nyx[1]-1)/(p2-1);
  const int nelz = (nzx[1]-1)/(p3-1);

  // Add connectivity for second basis: local --> global node relation
  int i, j, k, iel;
  for (k = iel = 0; k < nelz; k++)
    for (j = 0; j < nely; j++)
      for (i = 0; i < nelx; i++, iel++)
      {
	size_t nen1 = myMNPC[iel].size();
	myMNPC[iel].resize(nen1+p1*p2*p3);

	// First node in current element
	int corner = nb[0] + (p3-1)*nxx[1]*nyx[1]*k + (p2-1)*nxx[1]*j + (p1-1)*i;

	for (size_t c = 0; c < p3; c++)
	{
	  int cornod = nen1 + p1*p2*c;
	  myMNPC[iel][cornod] = corner + nxx[1]*nyx[1]*c;
	  for (size_t b = 1; b < p2; b++)
	  {
	    int facenod = cornod + p1*b;
	    myMNPC[iel][facenod] = myMNPC[iel][cornod] + nxx[1]*b;
	    for (size_t a = 1; a < p1; a++)
	    {
	      myMNPC[iel][facenod+a] = myMNPC[iel][facenod] + a;
	      myMNPC[iel][cornod+a]  = myMNPC[iel][cornod] + a;
	    }
	  }
	}
      }

  return true;
}


bool ASMs3DmxLag::connectPatch (int face, ASM3D& neighbor, int nface,
                                int norient, int basis,
                                bool coordCheck, int thick)
{
  ASMs3DmxLag* neighMx = dynamic_cast<ASMs3DmxLag*>(&neighbor);
  if (!neighMx) return false;

  if (swapW && face > 4) // Account for swapped parameter direction
    face = 11-face;

  if (neighMx->swapW && face > 4) // Account for swapped parameter direction
    nface = 11-nface;

  size_t nb1=0;
  size_t nb2=0;
  for (size_t i = 1;i <= nxx.size(); ++i) {
    if (basis == 0 || i == (size_t)basis)
      if (!this->connectBasis(face,*neighMx,nface,norient,i,nb1,nb2,coordCheck,thick))
        return false;

    nb1 += nb[i-1];
    nb2 += neighMx->nb[i-1];
  }

  this->addNeighbor(neighMx);
  return true;
}


void ASMs3DmxLag::closeFaces (int dir, int, int)
{
  size_t nbi = 0;
  for (size_t i = 1;i <= nxx.size(); ++i) {
    this->ASMs3D::closeFaces(dir,i,nbi+1);
    nbi += nb[i-1];
  }
}


bool ASMs3DmxLag::getSize (int& n1, int& n2, int& n3, int basis) const
{
  if (basis <= 1)
    return this->ASMs3DLag::getSize(n1,n2,n3,1);

  n1 = nxx[1];
  n2 = nyx[1];
  n3 = nzx[1];

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

  std::vector<size_t> elem_size;
  for (size_t b = 0; b < nxx.size(); ++b)
    elem_size.push_back(elem_sizes[b][0]*elem_sizes[b][1]*elem_sizes[b][2]);

  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t g = 0; g < threadGroupsVol.size() && ok; g++)
  {
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < threadGroupsVol[g].size(); t++)
    {
      MxFiniteElement fe(elem_size);
      std::vector<Matrix> dNxdu;
      Matrix Xnod, Jac;
      Vec4   X;
      for (size_t l = 0; l < threadGroupsVol[g][t].size() && ok; l++)
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
        LocalIntegral* A = integrand.getLocalIntegral(elem_size,fe.iel,false);
        if (!integrand.initElement(MNPC[iel-1], elem_size, nb, *A))
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
              for (size_t b = 0; b < nxx.size(); ++b)
                if (!Lagrange::computeBasis(fe.basis(b+1),dNxdu[b],elem_sizes[b][0],
                                            x[i],elem_sizes[b][1],x[j],elem_sizes[b][2],x[k]))
                  ok = false;

              // Compute Jacobian inverse of coordinate mapping and derivatives
              fe.detJxW = utl::Jacobian(Jac,fe.grad(geoBasis),Xnod,dNxdu[geoBasis-1]);
              if (fe.detJxW == 0.0) continue; // skip singular points
              for (size_t b = 0; b < nxx.size(); ++b)
                if (b != (size_t)geoBasis-1)
                  fe.grad(b+1).multiply(dNxdu[b],Jac);

              // Cartesian coordinates of current integration point
              X = Xnod * fe.basis(geoBasis);
              X.t = time.t;

              // Evaluate the integrand and accumulate element contributions
              fe.detJxW *= w[i]*w[j]*w[k];
              if (!integrand.evalIntMx(*A,fe,time,X))
                ok = false;
            }

        // Finalize the element quantities
        if (ok && !integrand.finalizeElement(*A,time,firstIp+jp))
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


bool ASMs3DmxLag::integrate (Integrand& integrand, int lIndex,
			     GlobalIntegral& glInt,
			     const TimeDomain& time)
{
  if (!svol) return true; // silently ignore empty patches

  std::map<char,ThreadGroups>::const_iterator tit;
  if ((tit = threadGroupsFace.find(lIndex%10)) == threadGroupsFace.end())
  {
    std::cerr <<" *** ASMs3DmxLag::integrate: No thread groups for face "
              << lIndex%10 << std::endl;
    return false;
  }
  const ThreadGroups& threadGrp = tit->second;

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Find the parametric direction of the face normal {-3,-2,-1, 1, 2, 3}
  const int faceDir = (lIndex%10+1)/(lIndex%2 ? -2 : 2);

  const int t0 = abs(faceDir); // unsigned normal direction of the face
  const int t1 = 1 + t0%3; // first tangent direction of the face
  const int t2 = 1 + t1%3; // second tangent direction of the face

  // Extract the Neumann order flag (1 or higher) for the integrand
  integrand.setNeumannOrder(1 + lIndex/10);

  std::vector<size_t> elem_size;
  for (size_t b = 0; b < nxx.size(); ++b)
    elem_size.push_back(elem_sizes[b][0]*elem_sizes[b][1]*elem_sizes[b][2]);

  // Number of elements in each direction
  const int nel1 = (nxx[geoBasis-1]-1)/(elem_sizes[geoBasis-1][0]-1);
  const int nel2 = (nyx[geoBasis-1]-1)/(elem_sizes[geoBasis-1][1]-1);

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex%10);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;

  // === Assembly loop over all elements on the patch face =====================

  bool ok = true;
  for (size_t g = 0; g < threadGrp.size() && ok; g++)
  {
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < threadGrp[g].size(); t++)
    {
      MxFiniteElement fe(elem_size);
      std::vector<Matrix> dNxdu;
      Matrix Xnod, Jac;
      Vec4   X;
      Vec3   normal;
      double xi[3];

      for (size_t l = 0; l < threadGrp[g][t].size() && ok; l++)
      {
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
        LocalIntegral* A = integrand.getLocalIntegral(elem_size,fe.iel,true);
	if (!integrand.initElementBou(MNPC[iel-1],elem_size,nb,*A))
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
            for (size_t b = 0; b < nxx.size(); ++b)
              if (!Lagrange::computeBasis(fe.basis(b+1), dNxdu[b],elem_sizes[b][0],xi[0],
                                          elem_sizes[b][1],xi[1],elem_sizes[b][2],xi[2]))
                ok = false;

	    // Compute basis function derivatives and the edge normal
	    fe.detJxW = utl::Jacobian(Jac,normal,fe.grad(geoBasis),Xnod,dNxdu[geoBasis-1],t1,t2);
	    if (fe.detJxW == 0.0) continue; // skip singular points
            for (size_t b = 0; b < nxx.size(); ++b)
              if (b != (size_t)geoBasis-1)
                fe.grad(b+1).multiply(dNxdu[b],Jac);

	    if (faceDir < 0) normal *= -1.0;

	    // Cartesian coordinates of current integration point
	    X = Xnod * fe.basis(geoBasis);
	    X.t = time.t;

	    // Evaluate the integrand and accumulate element contributions
	    fe.detJxW *= wg[i]*wg[j];
	    if (!integrand.evalBouMx(*A,fe,time,X,normal))
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


bool ASMs3DmxLag::evalSolution (Matrix& sField, const Vector& locSol,
                                const RealArray*, bool, int) const
{
  size_t nc1 = nfx[0];
  size_t nc2 = 0;
  if (nc1*nb[0] < locSol.size())
    nc2 = (locSol.size() - nc1*nb[0])/nb[1];
  else
    nc1 = locSol.size()/nb[0];

  if (nc1*nb[0] + nc2*nb[1] != locSol.size())
    return false;

  // TODO: Add evaluation of the second field at the nodes of the first field
  size_t nPoints = nb[0];
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

  const size_t p1 = svol->order(0);
  const size_t p2 = svol->order(1);
  const size_t p3 = svol->order(2);
  const size_t q1 = p1 - 1;
  const size_t q2 = p2 - 1;
  const size_t q3 = p3 - 1;
  double incx = 2.0/double(q1);
  double incy = 2.0/double(q2);
  double incz = 2.0/double(q3);

  size_t nPoints = nb[0];
  IntVec check(nPoints,0);

  std::vector<size_t> elem_size;
  for (size_t b = 0; b < nxx.size(); ++b)
    elem_size.push_back(elem_sizes[b][0]*elem_sizes[b][1]*elem_sizes[b][2]);

  MxFiniteElement fe(elem_size);
  Vector          solPt;
  Vectors         globSolPt(nPoints);
  std::vector<Matrix> dNxdu(nxx.size());
  Matrix          Xnod, Jac;

  // Evaluate the secondary solution field at each point
  for (size_t iel = 1; iel <= nel; iel++)
  {
    IntVec::const_iterator f2start = geoBasis == 1? MNPC[iel-1].begin() :
                                     MNPC[iel-1].begin() +
                                     std::accumulate(elem_size.begin()+geoBasis-2,
                                                     elem_size.begin()+geoBasis-1, 0);
    IntVec::const_iterator f2end = f2start + elem_size[geoBasis-1];
    IntVec mnpc1(f2start,f2end);

    this->getElementCoordinates(Xnod,iel);

    size_t i, j, k, loc = 0;
    for (k = 0; k < p3; k++)
      for (j = 0; j < p2; j++)
	for (i = 0; i < p1; i++, loc++)
	{
	  fe.xi   = -1.0 + i*incx;
	  fe.eta  = -1.0 + j*incy;
	  fe.zeta = -1.0 + k*incz;
          for (size_t b = 0; b < nxx.size(); ++b)
            if (!Lagrange::computeBasis(fe.basis(b+1),dNxdu[b],elem_sizes[b][0],fe.xi,
                                        elem_sizes[b][1],fe.eta,elem_sizes[b][2],fe.zeta))
              return false;

	  // Compute the Jacobian inverse
          fe.detJxW = utl::Jacobian(Jac,fe.grad(geoBasis),Xnod,dNxdu[geoBasis-1]);

          for (size_t b = 1; b <= nxx.size(); b++)
            if (b != (size_t)geoBasis)
            {
              if (fe.detJxW == 0.0)
                fe.grad(b).clear();
              else
                fe.grad(b).multiply(dNxdu[b-1],Jac);
            }

	  // Now evaluate the solution field
	  if (!integrand.evalSol(solPt,fe,Xnod*fe.basis(geoBasis),MNPC[iel-1],elem_size,nb))
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
    sField.fillColumn(1+i,globSolPt[i] /= check[i]);

  return true;
}
