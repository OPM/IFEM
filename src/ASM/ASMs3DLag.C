// $Id$
//==============================================================================
//!
//! \file ASMs3DLag.C
//!
//! \date Feb 10 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for assembly of structured 3D Lagrange FE models.
//!
//==============================================================================

#include "GoTools/trivariate/SplineVolume.h"

#include "ASMs3DLag.h"
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
#include <array>


ASMs3DLag::ASMs3DLag (unsigned char n_f) : ASMs3D(n_f), coord(myCoord)
{
  nx = ny = nz = 0;
}


ASMs3DLag::ASMs3DLag (const ASMs3DLag& patch, unsigned char n_f)
  : ASMs3D(patch,n_f), coord(patch.myCoord)
{
  nx = patch.nx;
  ny = patch.ny;
  nz = patch.nz;
}


ASMs3DLag::ASMs3DLag (const ASMs3DLag& patch)
  : ASMs3D(patch), coord(myCoord)
{
  nx = patch.nx;
  ny = patch.ny;
  nz = patch.nz;
  myCoord = patch.coord;
}


void ASMs3DLag::clear (bool retainGeometry)
{
  myCoord.clear();
  nx = ny = nz = 0;

  this->ASMs3D::clear(retainGeometry);
}


bool ASMs3DLag::addXElms (short int dim, short int item, size_t nXn,
                          IntVec& nodes)
{
  if (!this->addXNodes(dim,nXn,nodes))
    return false;

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  const int nel1 = (nx-1)/(p1-1);
  const int nel2 = (ny-1)/(p2-1);
  const int nel3 = (nz-1)/(p3-1);

  int iel = 0;
  bool skipMe = false;
  for (int i3 = 0; i3 < nel3; i3++)
    for (int i2 = 0; i2 < nel2; i2++)
      for (int i1 = 0; i1 < nel1; i1++, iel++)
      {
        if (MLGE[iel] < 1) continue; // Skip zero-volume element

        // Skip elements that are not on current boundary face
        switch (item)
          {
          case 1: skipMe = i1 > 0;      break;
          case 2: skipMe = i1 < nel1-1; break;
          case 3: skipMe = i2 > 0;      break;
          case 4: skipMe = i2 < nel2-1; break;
          case 5: skipMe = i3 > 0;      break;
          case 6: skipMe = i3 < nel3-1; break;
          }
        if (skipMe) continue;

        IntVec& mnpc = myMNPC[nel+iel];
        if (!mnpc.empty())
        {
          std::cerr <<" *** ASMs3DLag::addXElms: Only one X-face allowed."
                    << std::endl;
          return false;
        }

        mnpc = MNPC[iel]; // Copy the ordinary element nodes

        // Negate node numbers that are not on the boundary face, to flag that
        // they shall not receive any tangent and/or residual contributions
        int lnod = 0;
        for (int j3 = 0; j3 < p3; j3++)
          for (int j2 = 0; j2 < p2; j2++)
            for (int j1 = 0; j1 < p1; j1++, lnod++)
            {
              switch (item)
                {
                case 1: skipMe = j1 > 0;    break;
                case 2: skipMe = j1 < p1-1; break;
                case 3: skipMe = j2 > 0;    break;
                case 4: skipMe = j2 < p2-1; break;
                case 5: skipMe = j3 > 0;    break;
                case 6: skipMe = j3 < p3-1; break;
                }
              if (skipMe) // Hack for node 0: Using -maxint as flag instead
                mnpc[lnod] = mnpc[lnod] == 0 ? -2147483648 : -mnpc[lnod];
            }

        // Add connectivity to the extra-ordinary nodes
        for (size_t i = 0; i < nXn; i++)
          mnpc.push_back(MLGN.size()-nXn+i);

        myMLGE[nel+iel] = -(++gEl); // Extra-ordinary element => negative sign
      }

  return true;
}


bool ASMs3DLag::generateFEMTopology ()
{
  if (!svol) return false;

  // Order of basis in the three parametric directions (order = degree + 1)
  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  // Evaluate the parametric values
  RealArray gpar1, gpar2, gpar3;
  if (!this->getGridParameters(gpar1,0,p1-1)) return false;
  if (!this->getGridParameters(gpar2,1,p2-1)) return false;
  if (!this->getGridParameters(gpar3,2,p3-1)) return false;

  // Number of nodes in each direction
  nx = gpar1.size();
  ny = gpar2.size();
  nz = gpar3.size();
  // Number of nodes in the patch
  nnod = nx*ny*nz;

  if (!coord.empty())
    return coord.size() == nnod;

  // Number of elements in each direction
  const int nelx = (nx-1)/(p1-1);
  const int nely = (ny-1)/(p2-1);
  const int nelz = (nz-1)/(p3-1);

  // Evaluate the nodal coordinates in the physical space
  RealArray XYZ(svol->dimension()*nnod);
  svol->gridEvaluator(gpar1,gpar2,gpar3,XYZ);

  size_t i1, j1;
  myMLGN.resize(nnod);
  myCoord.resize(nnod);
  for (i1 = j1 = 0; i1 < coord.size(); i1++)
  {
    myMLGN[i1] = ++gNod;
    myCoord[i1][0] = XYZ[j1++];
    myCoord[i1][1] = XYZ[j1++];
    myCoord[i1][2] = XYZ[j1++];
  }

  // Number of elements in patch
  nel = nelx*nely*nelz;
  // Number of nodes per element
  const int nen = p1*p2*p3;
  // Number of nodes in a xy-surface of an element
  const int ct  = p1*p2;

  // Connectivity array: local --> global node relation
  myMLGE.resize(nel);
  myMNPC.resize(nel);

  int i, j, k, a, b, c, iel = 0;
  for (k = 0; k < nelz; k++)
    for (j = 0; j < nely; j++)
      for (i = 0; i < nelx; i++, iel++)
      {
	myMLGE[iel] = ++gEl;
	myMNPC[iel].resize(nen);
	// First node in current element
	int corner = (p3-1)*(nx*ny)*k + (p2-1)*nx*j + (p1-1)*i;

	for (c = 0; c < p3; c++)
	{
	  int cornod = ct*c;
	  myMNPC[iel][cornod] = corner + c*nx*ny;
	  for (b = 1; b < p2; b++)
	  {
	    int facenod = cornod + b*p1;
	    myMNPC[iel][facenod] = myMNPC[iel][cornod] + b*nx;
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


Vec3 ASMs3DLag::getCoord (size_t inod) const
{
  if (inod < 1 || inod > coord.size()) return Vec3();

  return coord[inod-1];
}


bool ASMs3DLag::getElementCoordinates (Matrix& X, int iel) const
{
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs3DLag::getElementCoordinates: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return false;
  }

  // Number of nodes per element
  const size_t nen = svol->order(0)*svol->order(1)*svol->order(2);

  X.resize(3,nen);
  for (size_t i = 0; i < nen; i++)
    X.fillColumn(i+1,coord[MNPC[iel-1][i]].ptr());

  return true;
}


void ASMs3DLag::getNodalCoordinates (Matrix& X) const
{
  X.resize(3,coord.size());

  for (size_t inod = 0; inod < coord.size(); inod++)
    X.fillColumn(inod+1,coord[inod].ptr());
}


bool ASMs3DLag::updateCoords (const Vector& displ)
{
  if (shareFE) return true;

  if (displ.size() != 3*coord.size())
  {
    std::cerr <<" *** ASMs3DLag::updateCoords: Invalid dimension "
	      << displ.size() <<" on displ, should be "
	      << 3*coord.size() << std::endl;
    return false;
  }

  const double* u = displ.ptr();
  for (size_t inod = 0; inod < myCoord.size(); inod++, u += 3)
    myCoord[inod] += RealArray(u,u+3);

  return true;
}


bool ASMs3DLag::getSize (int& n1, int& n2, int& n3, int) const
{
  n1 = nx;
  n2 = ny;
  n3 = nz;

  return true;
}


size_t ASMs3DLag::getNoBoundaryElms (char lIndex, char ldim) const
{
  if (!svol) return 0;

  if (ldim < 1 && lIndex > 0)
    return 1;

  int nel[3]; // Number of elements in each direction
  nel[0] = (nx-1)/(svol->order(0)-1);
  nel[1] = (ny-1)/(svol->order(1)-1);
  nel[2] = (nz-1)/(svol->order(2)-1);

  if (ldim < 2 && lIndex > 0 && lIndex <= 12)
    return nel[(lIndex-1)/4];

  switch (lIndex)
    {
    case 1:
    case 2:
      return nel[1]*nel[2];
    case 3:
    case 4:
      return nel[0]*nel[2];
    case 5:
    case 6:
      return nel[0]*nel[1];
    }

  return 0;
}


bool ASMs3DLag::integrate (Integrand& integrand,
			   GlobalIntegral& glInt,
			   const TimeDomain& time)
{
  if (!svol) return true; // silently ignore empty patches

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Get the reduced integration quadrature points, if needed
  const double* xr = nullptr;
  const double* wr = nullptr;
  int nRed = integrand.getReducedIntegration(nGauss);
  if (nRed > 0)
  {
    xr = GaussQuadrature::getCoord(nRed);
    wr = GaussQuadrature::getWeight(nRed);
    if (!xr || !wr) return false;
  }
  else if (nRed < 0)
    nRed = nGauss; // The integrand needs to know nGauss

  // Get parametric coordinates of the elements
  RealArray upar, vpar, wpar;
  this->getGridParameters(upar,0,1);
  this->getGridParameters(vpar,1,1);
  this->getGridParameters(wpar,2,1);

  // Number of elements in each direction
  const int nel1 = upar.size() - 1;
  const int nel2 = vpar.size() - 1;

  // Order of basis in the three parametric directions (order = degree + 1)
  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);


  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t g = 0; g < threadGroupsVol.size() && ok; g++)
  {
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < threadGroupsVol[g].size(); t++)
    {
      FiniteElement fe(p1*p2*p3);
      Matrix dNdu, Xnod, Jac;
      Vec4   X;
      for (size_t l = 0; l < threadGroupsVol[g][t].size() && ok; l++)
      {
        int iel = threadGroupsVol[g][t][l];
        int i1  =  iel % nel1;
        int i2  = (iel / nel1) % nel2;
        int i3  =  iel / (nel1*nel2);

        // Set up nodal point coordinates for current element
        if (!this->getElementCoordinates(Xnod,++iel))
        {
          ok = false;
          break;
        }

        if (integrand.getIntegrandType() & Integrand::ELEMENT_CENTER)
        {
          // Compute the element "center" (average of element node coordinates)
          X = 0.0;
          for (size_t i = 1; i <= 3; i++)
            for (size_t j = 1; j <= Xnod.cols(); j++)
              X[i-1] += Xnod(i,j);

          X *= 1.0/(double)Xnod.cols();
        }

        // Initialize element quantities
        fe.iel = MLGE[iel-1];
        LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel);
        if (!integrand.initElement(MNPC[iel-1],fe,X,nRed*nRed*nRed,*A))
        {
          A->destruct();
          ok = false;
          break;
        }

        if (xr)
        {
          // --- Selective reduced integration loop ----------------------------

          for (int k = 0; k < nRed; k++)
            for (int j = 0; j < nRed; j++)
              for (int i = 0; i < nRed; i++)
              {
                // Local element coordinates of current integration point
                fe.xi   = xr[i];
                fe.eta  = xr[j];
                fe.zeta = xr[k];

                // Parameter value of current integration point
                fe.u = 0.5*(upar[i1]*(1.0-xr[i]) + upar[i1+1]*(1.0+xr[i]));
                fe.v = 0.5*(vpar[i2]*(1.0-xr[j]) + vpar[i2+1]*(1.0+xr[j]));
                fe.w = 0.5*(wpar[i3]*(1.0-xr[k]) + wpar[i3+1]*(1.0+xr[k]));

                // Compute basis function derivatives at current point
                // using tensor product of one-dimensional Lagrange polynomials
                if (!Lagrange::computeBasis(fe.N,dNdu,
                                            p1,xr[i],p2,xr[j],p3,xr[k]))
                {
                  ok = false;
                  break;
                }

                // Compute Jacobian inverse and derivatives
                fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

                // Cartesian coordinates of current integration point
                X = Xnod * fe.N;
                X.t = time.t;

                // Compute the reduced integration terms of the integrand
                fe.detJxW *= wr[i]*wr[j]*wr[k];
                if (!integrand.reducedInt(*A,fe,X))
                  ok = false;
              }
        }


        // --- Integration loop over all Gauss points in each direction --------

        int jp = ((i3*nel2 + i2)*nel1 + i1)*nGauss*nGauss*nGauss;
        fe.iGP = firstIp + jp; // Global integration point counter

        for (int k = 0; k < nGauss; k++)
          for (int j = 0; j < nGauss; j++)
            for (int i = 0; i < nGauss; i++, fe.iGP++)
            {
              // Local element coordinates of current integration point
              fe.xi   = xg[i];
              fe.eta  = xg[j];
              fe.zeta = xg[k];

              // Parameter value of current integration point
              fe.u = 0.5*(upar[i1]*(1.0-xg[i]) + upar[i1+1]*(1.0+xg[i]));
              fe.v = 0.5*(vpar[i2]*(1.0-xg[j]) + vpar[i2+1]*(1.0+xg[j]));
              fe.w = 0.5*(wpar[i3]*(1.0-xg[k]) + wpar[i3+1]*(1.0+xg[k]));

              // Compute basis function derivatives at current integration point
              // using tensor product of one-dimensional Lagrange polynomials
              if (!Lagrange::computeBasis(fe.N,dNdu,p1,xg[i],p2,xg[j],p3,xg[k]))
                ok = false;

              // Compute Jacobian inverse of coordinate mapping and derivatives
              fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);
              if (fe.detJxW == 0.0) continue; // skip singular points

              // Cartesian coordinates of current integration point
              X = Xnod * fe.N;
              X.t = time.t;

              // Evaluate the integrand and accumulate element contributions
              fe.detJxW *= wg[i]*wg[j]*wg[k];
              if (!integrand.evalInt(*A,fe,time,X))
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


bool ASMs3DLag::integrate (Integrand& integrand, int lIndex,
			   GlobalIntegral& glInt,
			   const TimeDomain& time)
{
  if (!svol) return true; // silently ignore empty patches

  std::map<char,ThreadGroups>::const_iterator tit;
  if ((tit = threadGroupsFace.find(lIndex%10)) == threadGroupsFace.end())
  {
    std::cerr <<" *** ASMs3DLag::integrate: No thread groups for face "
              << lIndex%10 << std::endl;
    return false;
  }
  const ThreadGroups& threadGrp = tit->second;

  // Get Gaussian quadrature points and weights
  int nGP = integrand.getBouIntegrationPoints(nGauss);
  const double* xg = GaussQuadrature::getCoord(nGP);
  const double* wg = GaussQuadrature::getWeight(nGP);
  if (!xg || !wg) return false;

  // Find the parametric direction of the face normal {-3,-2,-1, 1, 2, 3}
  const int faceDir = (lIndex%10+1)/(lIndex%2 ? -2 : 2);

  const int t0 = abs(faceDir); // unsigned normal direction of the face
  const int t1 = 1 + t0%3; // first tangent direction of the face
  const int t2 = 1 + t1%3; // second tangent direction of the face

  // Order of basis in the three parametric directions (order = degree + 1)
  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  // Number of elements in each direction
  const int nel1 = (nx-1)/(p1-1);
  const int nel2 = (ny-1)/(p2-1);
  const int nel3 = (nz-1)/(p3-1);

  // Get parametric coordinates of the elements
  RealArray upar, vpar, wpar;
  if (t0 == 1)
    upar.resize(1,faceDir < 0 ? svol->startparam(0) : svol->endparam(0));
  else if (t0 == 2)
    vpar.resize(1,faceDir < 0 ? svol->startparam(1) : svol->endparam(1));
  else if (t0 == 3)
    wpar.resize(1,faceDir < 0 ? svol->startparam(2) : svol->endparam(2));

  if (upar.empty()) this->getGridParameters(upar,0,1);
  if (vpar.empty()) this->getGridParameters(vpar,1,1);
  if (wpar.empty()) this->getGridParameters(wpar,2,1);

  // Extract the Neumann order flag (1 or higher) for the integrand
  integrand.setNeumannOrder(1 + lIndex/10);

  // Integrate the extraordinary elements?
  size_t doXelms = 0;
  if (integrand.getIntegrandType() & Integrand::XO_ELEMENTS)
    if ((doXelms = nel1*nel2*nel3)*2 > MNPC.size())
    {
      std::cerr <<" *** ASMs3DLag::integrate: Too few XO-elements "
                << MNPC.size() - doXelms << std::endl;
      return false;
    }

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex%10);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;


  // === Assembly loop over all elements on the patch face =====================

  bool ok = true;
  for (size_t g = 0; g < threadGrp.size() && ok; g++)
  {
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < threadGrp[g].size(); t++)
    {
      FiniteElement fe(p1*p2*p3);
      fe.u = upar.front();
      fe.v = vpar.front();
      fe.w = wpar.front();

      Matrix dNdu, Xnod, Jac;
      Vec4   X;
      Vec3   normal;
      double xi[3];

      for (size_t l = 0; l < threadGrp[g][t].size() && ok; l++)
      {
        int iel = threadGrp[g][t][l];
        int i1  =  iel % nel1;
        int i2  = (iel / nel1) % nel2;
        int i3  =  iel / (nel1*nel2);

        // Set up nodal point coordinates for current element
        if (!this->getElementCoordinates(Xnod,++iel))
        {
          ok = false;
          break;
        }

        // Initialize element quantities
        fe.iel = abs(MLGE[doXelms+iel-1]);
        LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel,true);
        if (!integrand.initElementBou(MNPC[doXelms+iel-1],*A))
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

        int k1, k2, k3;
        int jp = (j2*nf1 + j1)*nGP*nGP;
        fe.iGP = firstp + jp; // Global integration point counter

	for (int j = 0; j < nGP; j++)
	  for (int i = 0; i < nGP; i++, fe.iGP++)
	  {
	    // Local element coordinates of current integration point
	    xi[t0-1] = faceDir < 0 ? -1.0 : 1.0;
	    xi[t1-1] = xg[i];
	    xi[t2-1] = xg[j];
	    fe.xi   = xi[0];
	    fe.eta  = xi[1];
	    fe.zeta = xi[2];

	    // Local element coordinates and parameter values
	    // of current integration point
	    switch (abs(faceDir))
	    {
	      case 1: k2 = i; k3 = j; k1 = -1; break;
	      case 2: k1 = i; k3 = j; k2 = -1; break;
	      case 3: k1 = i; k2 = j; k3 = -1; break;
	      default: k1 = k2 = k3 = -1;
	    }
	    if (upar.size() > 1)
	      fe.u = 0.5*(upar[i1]*(1.0-xg[k1]) + upar[i1+1]*(1.0+xg[k1]));
	    if (vpar.size() > 1)
	      fe.v = 0.5*(vpar[i2]*(1.0-xg[k2]) + vpar[i2+1]*(1.0+xg[k2]));
	    if (wpar.size() > 1)
	      fe.w = 0.5*(wpar[i3]*(1.0-xg[k3]) + wpar[i3+1]*(1.0+xg[k3]));

	    // Compute the basis functions and their derivatives, using
	    // tensor product of one-dimensional Lagrange polynomials
            if (!Lagrange::computeBasis(fe.N,dNdu,p1,xi[0],p2,xi[1],p3,xi[2]))
              ok = false;

	    // Compute basis function derivatives and the face normal
	    fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,Xnod,dNdu,t1,t2);
	    if (fe.detJxW == 0.0) continue; // skip singular points

	    if (faceDir < 0) normal *= -1.0;

	    // Cartesian coordinates of current integration point
	    X = Xnod * fe.N;
	    X.t = time.t;

	    // Evaluate the integrand and accumulate element contributions
	    fe.detJxW *= wg[i]*wg[j];
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


bool ASMs3DLag::integrateEdge (Integrand& integrand, int lEdge,
			       GlobalIntegral& glInt,
			       const TimeDomain& time)
{
  if (!svol) return true; // silently ignore empty patches

  // Parametric direction of the edge {0, 1, 2}
  const int lDir = (lEdge-1)/4;

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Order of basis in the three parametric directions (order = degree + 1)
  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  // Number of elements in each direction
  const int nelx = (nx-1)/(p1-1);
  const int nely = (ny-1)/(p2-1);
  const int nelz = (nz-1)/(p3-1);

  FiniteElement fe(p1*p2*p3);
  Matrix dNdu, Xnod, Jac;
  Vec4   X;
  Vec3   tangent;
  double xi[3];

  switch (lEdge)
    {
    case  1: xi[1] = -1.0; xi[2] = -1.0; break;
    case  2: xi[1] =  1.0; xi[2] = -1.0; break;
    case  3: xi[1] = -1.0; xi[2] =  1.0; break;
    case  4: xi[1] =  1.0; xi[2] =  1.0; break;
    case  5: xi[0] = -1.0; xi[2] = -1.0; break;
    case  6: xi[0] =  1.0; xi[2] = -1.0; break;
    case  7: xi[0] = -1.0; xi[2] =  1.0; break;
    case  8: xi[0] =  1.0; xi[2] =  1.0; break;
    case  9: xi[0] = -1.0; xi[1] = -1.0; break;
    case 10: xi[0] =  1.0; xi[1] = -1.0; break;
    case 11: xi[0] = -1.0; xi[1] =  1.0; break;
    case 12: xi[0] =  1.0; xi[1] =  1.0; break;
    }

  std::map<char,size_t>::const_iterator iit = firstBp.find(lEdge);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;


  // === Assembly loop over all elements on the patch edge =====================

  int ip, iel = 1;
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

	if (lEdge < 5)
	  ip = i1*nGauss;
	else if (lEdge < 9)
	  ip = i2*nGauss;
	else
	  ip = i3*nGauss;

	// Set up nodal point coordinates for current element
	if (!this->getElementCoordinates(Xnod,iel)) return false;

	// Initialize element quantities
	fe.iel = MLGE[iel-1];
        LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel,true);
        bool ok = integrand.initElementBou(MNPC[iel-1],*A);


	// --- Integration loop over all Gauss points along the edge -----------

	fe.iGP = firstp + ip; // Global integration point counter

	for (int i = 0; i < nGauss && ok; i++, fe.iGP++)
	{
	  // Gauss point coordinates on the edge
	  xi[lDir] = xg[i];

	  // Compute the basis functions and their derivatives, using
	  // tensor product of one-dimensional Lagrange polynomials
	  if (!Lagrange::computeBasis(fe.N,dNdu,p1,xi[0],p2,xi[1],p3,xi[2]))
	    ok = false;

	  // Compute basis function derivatives and the edge tangent
	  fe.detJxW = utl::Jacobian(Jac,tangent,fe.dNdX,Xnod,dNdu,1+lDir);
	  if (fe.detJxW == 0.0) continue; // skip singular points

	  // Cartesian coordinates of current integration point
	  X = Xnod * fe.N;
	  X.t = time.t;

	  // Evaluate the integrand and accumulate element contributions
          if (!integrand.evalBou(*A,fe,time,X,tangent))
	    ok = false;
	}

        // Finalize the element quantities
        if (ok && !integrand.finalizeElementBou(*A,fe,time))
          ok = false;

	// Assembly of global system integral
	if (ok && !glInt.assemble(A->ref(),fe.iel))
	  ok = false;

	A->destruct();

	if (!ok) return false;
      }

  return true;
}


int ASMs3DLag::evalPoint (const double* xi, double* param, Vec3& X) const
{
  if (!svol) return -3;

  // Evaluate the parametric values of the point and nodes
  std::array<RealArray,3> u;
  for (int d = 0; d < 3; d++)
  {
    param[d] = (1.0-xi[d])*svol->startparam(d) + xi[d]*svol->endparam(d);
    if (!this->getGridParameters(u[d],d,svol->order(d)-1)) return -3;
  }

  // Search for the closest node
  size_t i = utl::find_closest(u[0],param[0]);
  size_t j = utl::find_closest(u[1],param[1]);
  size_t k = utl::find_closest(u[2],param[2]);
  size_t n = u[0].size()*(u[1].size()*k + j) + i;
  X = coord[n];

  return 1+n;
}


bool ASMs3DLag::tesselate (ElementBlock& grid, const int* npe) const
{
  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);
  if (p1 != npe[0] || p2 != npe[1] || p2 != npe[2])
  {
    int* newnpe = const_cast<int*>(npe);
    std::cout <<"\nLagrange elements: The number of visualization points are "
	      << p1 <<" "<< p2 <<" "<< p3 <<" by default\n"<< std::endl;
    newnpe[0] = p1;
    newnpe[1] = p2;
    newnpe[2] = p3;
  }

  if (!this->ASMs3D::tesselate(grid,npe))
    return false;

  // Adjust the element Id since each Lagrange element covers several knot-spans
  int i, ie, nse1 = p1-1;
  int j, je, nse2 = p2-1;
  int k, ke, nse3 = p3-1;
  int nelx = (nx-1)/nse1;
  int nely = (ny-1)/nse1;
  for (k = ke = 1; k < (int)nz; k++)
  {
    for (j = je = 1; j < (int)ny; j++)
    {
      for (i = ie = 1; i < (int)nx; i++)
      {
	grid.setElmId(((k-1)*(ny-1)+j-1)*(nx-1)+i,((ke-1)*nely+je-1)*nelx+ie);
	if (i%nse1 == 0) ie++;
      }
      if (j%nse2 == 0) je++;
    }
    if (k%nse3 == 0) ke++;
  }

  return true;
}


bool ASMs3DLag::evalSolution (Matrix& sField, const Vector& locSol,
			      const int*) const
{
  return this->evalSolution(sField,locSol,(const RealArray*)0,true);
}


bool ASMs3DLag::evalSolution (Matrix& sField, const Vector& locSol,
                              const RealArray*, bool, int) const
{
  size_t nPoints = coord.size();
  size_t nNodes = this->getNoNodes(-1);
  size_t nComp = locSol.size() / nNodes;
  if (nNodes < nPoints || nComp*nNodes != locSol.size())
    return false;

  sField.resize(nComp,nPoints);
  const double* u = locSol.ptr();
  for (size_t n = 1; n <= nPoints; n++, u += nComp)
    sField.fillColumn(n,u);

  return true;
}


bool ASMs3DLag::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			      const int*, bool) const
{
  return this->evalSolution(sField,integrand,(const RealArray*)0,true);
}


bool ASMs3DLag::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			      const RealArray*, bool) const
{
  sField.resize(0,0);

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);
  double incx = 2.0/double(p1-1);
  double incy = 2.0/double(p2-1);
  double incz = 2.0/double(p3-1);

  size_t nPoints = coord.size();
  IntVec check(nPoints,0);

  FiniteElement fe(p1*p2*p3);
  Vector        solPt;
  Vectors       globSolPt(nPoints);
  Matrix        dNdu, Xnod, Jac;

  // Evaluate the secondary solution field at each point
  const int nel = this->getNoElms(true);
  for (int iel = 1; iel <= nel; iel++)
  {
    const IntVec& mnpc = MNPC[iel-1];
    this->getElementCoordinates(Xnod,iel);

    int i, j, k, loc = 0;
    for (k = 0; k < p3; k++)
      for (j = 0; j < p2; j++)
	for (i = 0; i < p1; i++, loc++)
	{
	  fe.xi   = -1.0 + i*incx;
	  fe.eta  = -1.0 + j*incy;
	  fe.zeta = -1.0 + k*incz;
	  if (!Lagrange::computeBasis(fe.N,dNdu,p1,fe.xi,p2,fe.eta,p3,fe.zeta))
	    return false;

	  // Compute the Jacobian inverse
	  fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

	  // Now evaluate the solution field
	  if (!integrand.evalSol(solPt,fe,Xnod*fe.N,mnpc))
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


void ASMs3DLag::generateThreadGroups (const Integrand&, bool, bool)
{
  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  const int nel1 = (nx-1)/(p1-1);
  const int nel2 = (ny-1)/(p2-1);
  const int nel3 = (nz-1)/(p3-1);

  threadGroupsVol.calcGroups(nel1,nel2,nel3,1);
}


void ASMs3DLag::generateThreadGroups (char lIndex, bool, bool)
{
  if (threadGroupsFace.find(lIndex) != threadGroupsFace.end()) return;

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);
  const int n1 = (nx-1)/(p1-1);
  const int n2 = (ny-1)/(p2-1);
  const int n3 = (nz-1)/(p3-1);

  // Find elements that are on the boundary face 'lIndex'
  IntVec map; map.reserve(this->getNoBoundaryElms(lIndex,2));
  int d1, d2, iel = 0;
  for (int i3 = 1; i3 <= n3; i3++)
    for (int i2 = 1; i2 <= n2; i2++)
      for (int i1 = 1; i1 <= n1; i1++, iel++)
	switch (lIndex)
          {
	  case 1: if (i1 ==  1) map.push_back(iel); break;
	  case 2: if (i1 == n1) map.push_back(iel); break;
	  case 3: if (i2 ==  1) map.push_back(iel); break;
	  case 4: if (i2 == n2) map.push_back(iel); break;
	  case 5: if (i3 ==  1) map.push_back(iel); break;
	  case 6: if (i3 == n3) map.push_back(iel); break;
          }

  switch (lIndex)
    {
    case 1:
    case 2:
      d1 = n2;
      d2 = n3;
      break;
    case 3:
    case 4:
      d1 = n1;
      d2 = n3;
      break;
    default:
      d1 = n1;
      d2 = n2;
    }

  threadGroupsFace[lIndex].calcGroups(d1,d2,1);
  threadGroupsFace[lIndex].applyMap(map);
}
