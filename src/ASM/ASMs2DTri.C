// $Id$
//==============================================================================
//!
//! \file ASMs2DTri.C
//!
//! \date Feb 07 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of structured 2D triangle-based FE models.
//!
//==============================================================================

#include "GoTools/geometry/SplineSurface.h"

#include "ASMs2DTri.h"
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "TriangleQuadrature.h"
#include "GaussQuadrature.h"
#include "ElementBlock.h"
#include "Vec3Oper.h"
#include <numeric>


bool ASMs2DTri::generateFEMTopology ()
{
  // Generate the tensorial quadrilateral mesh first
  bool wasEmpty = MLGE.empty();
  if (!this->ASMs2DLag::generateFEMTopology())
    return false;
  else if (!wasEmpty)
    return true;

  if (surf->order_u() != 2 || surf->order_v() != 2)
  {
    std::cerr <<" *** ASMs2DTri::generateFEMTopology:"
              <<" Currently implemented for linear triangles only."<< std::endl;
    return false;
  }

  // Double the number of elements in the patch
  IntMat qMNPC(myMNPC);
  size_t nquad = nel;
  nel = 2*nquad;
  myMNPC.resize(nel);
  myMLGE.resize(nel);
  std::iota(myMLGE.begin()+nquad,myMLGE.end(),gEl+1);
  gEl += nquad;

  // Generate the triangular mesh topology
  size_t i, j, iel;
  for (iel = i = 0; iel < nel; i++, iel++)
  {
    myMNPC[iel].resize(3);
    myMNPC[iel][0] = qMNPC[i][0];
    myMNPC[iel][1] = qMNPC[i][1];
    myMNPC[iel][2] = qMNPC[i][3];
    myMNPC[++iel].resize(3);
    myMNPC[iel][0] = qMNPC[i][3];
    myMNPC[iel][1] = qMNPC[i][2];
    myMNPC[iel][2] = qMNPC[i][0];
  }

  // Convert to Union Jack pattern if even number of elements
  if (nx%2 == 1 && ny%2 == 1)
    for (iel = j = 0; j+1 < ny; j++)
      for (i = 0; i+1 < nx; i++, iel += 2)
        if (i%2 + j%2 == 1)
        {
          myMNPC[iel][2] = myMNPC[iel+1][1];
          myMNPC[iel+1][2] = myMNPC[iel][1];
        }

  return true;
}


void ASMs2DTri::getNoBouPoints (size_t& nPt, char ldim, char lindx)
{
  size_t nGp = 2; // Always two gauss points along an edge (assuming linear)

  firstBp[lindx] = nPt;

  nPt += this->getNoBoundaryElms(lindx,ldim)*nGp; // Includes 0-span elements
}


bool ASMs2DTri::addXElms (short int, short int, size_t, IntVec&)
{
  std::cerr <<" *** ASMs2DTri::addXElms: Not implemented yet."<< std::endl;
  return false;
}


/*!
  \brief Finds the parameter value of an integration point of a triangle.
*/

static void evalParam (double& u, double& v, double A1, double A2,
                       const double* upar, const double* vpar,
                       int upper, bool flippedDiag)
{
  if (flippedDiag)
  {
    // The diagonal is LR to UL
    u = (1.0-A2)*upar[upper] + A2*upar[1-upper];
    v = (1.0-A2)*vpar[upper] + A2*vpar[1-upper];
  }
  else
  {
    // The diagonal is LL to UR
    u = A1*upar[upper] + (1.0-A1)*upar[1-upper];
    v = A1*vpar[upper] + (1.0-A1)*vpar[1-upper];
  }
}


/*!
  \brief Establishes matrices with basis functions and their 1st derivatives.
*/

static void evalBasis (Vector& N, Matrix& dNdu, double A1, double A2)
{
  switch (N.size()) {
  case 3:
    N(1) = A1;
    N(2) = A2;
    N(3) = 1.0 - A1 - A2;
    dNdu(1,1) = dNdu(2,2) =  1.0;
    dNdu(3,1) = dNdu(3,2) = -1.0;
    break;
  }
}


bool ASMs2DTri::integrate (Integrand& integrand,
                           GlobalIntegral& glInt,
                           const TimeDomain& time)
{
  if (!surf) return true; // silently ignore empty patches

  // Order of basis in the two parametric directions (order = degree + 1)
  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  // Sanity check, the two parameter direction must have the same order
  if (p1 != p2)
  {
    std::cerr <<" *** ASMs2DTri::integrate: Unequal orders, "
              << p1 <<" != " << p2 << std::endl;
    return false;
  }

  // Get quadrature points and weights
  const double* xg = TriangleQuadrature::getCoord(nGauss);
  const double* wg = TriangleQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Get the reduced integration quadrature points, if needed
  const double* xr = nullptr;
  const double* wr = nullptr;
  int nRed = integrand.getReducedIntegration(nGauss);
  if (nRed > 0)
  {
    xr = TriangleQuadrature::getCoord(nRed);
    wr = TriangleQuadrature::getWeight(nRed);
    if (!xr || !wr) return false;
  }
  else if (nRed < 0)
    nRed = nGauss; // The integrand needs to know nGauss

  // Get parametric coordinates of the grid points
  RealArray upar, vpar;
  this->getGridParameters(upar,0,1);
  this->getGridParameters(vpar,1,1);

  // Number of elements in each direction
  const int nelx = upar.size() - 1;
  const int nely = vpar.size() - 1;
  bool UnionJack = nelx%2 == 0 && nely%2 == 0;

  // Number of element nodes
  const size_t nen = p1*(p2+1)/2;


  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t g = 0; g < threadGroups.size() && ok; g++)
  {
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < threadGroups[g].size(); t++)
    {
      FiniteElement fe(nen);
      Matrix        dNdu(nen,2), Xnod, Jac;
      Vec4          X;
      for (size_t i = 0; i < threadGroups[g][t].size() && ok; i++)
      {
        int iel = threadGroups[g][t][i];
        int i1  = (iel/2) % nelx;
        int i2  = (iel/2) / nelx;
        bool flipD = UnionJack && i1%2 + i2%2 == 1;

        // Set up nodal point coordinates for current element
        if (!this->getElementCoordinates(Xnod,1+iel))
        {
          ok = false;
          break;
        }

        if (integrand.getIntegrandType() & Integrand::ELEMENT_CENTER)
        {
          // Compute the element "center" (average of element node coordinates)
          X = 0.0;
          for (size_t d = 1; d <= nsd; d++)
            for (size_t j = 1; j <= Xnod.cols(); j++)
              X[d-1] += Xnod(d,j);

          X *= 1.0/(double)Xnod.cols();
        }

        // Initialize element quantities
        fe.iel = MLGE[iel];
        LocalIntegral* A = integrand.getLocalIntegral(nen,fe.iel);
        if (!integrand.initElement(MNPC[iel],fe,X,nRed*nRed,*A))
        {
          A->destruct();
          ok = false;
          break;
        }


        if (xr) // --- Selective reduced integration loop ----------------------

          for (int j = 0; j < nRed; j++)
          {
            // Area coordinates and parameter value of current integration point
            fe.xi  = xr[2*j];
            fe.eta = xr[2*j+1];
            evalParam(fe.u,fe.v,fe.xi,fe.eta,&upar[i1],&vpar[i2],iel%2,flipD);

            // Compute basis function values and derivatives at current point
            evalBasis(fe.N,dNdu,fe.xi,fe.eta);

            // Compute Jacobian inverse and derivatives
            fe.detJxW = 0.5*utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

            // Cartesian coordinates of current integration point
            X = Xnod * fe.N;
            X.t = time.t;

            // Compute the reduced integration terms of the integrand
            fe.detJxW *= wr[j];
            if (!integrand.reducedInt(*A,fe,X))
              ok = false;
          }


        // --- Integration loop over all Gauss points --------------------------

        int jp = (i2*nelx + i1)*2*nGauss;
        fe.iGP = firstIp + jp; // Global integration point counter

        for (int j = 0; j < nGauss; j++, fe.iGP++)
        {
          // Area coordinates and parameter value of current integration point
          fe.xi  = xg[2*j];
          fe.eta = xg[2*j+1];
          evalParam(fe.u,fe.v,fe.xi,fe.eta,&upar[i1],&vpar[i2],iel%2,flipD);

          // Compute basis function values and derivatives at current point
          evalBasis(fe.N,dNdu,fe.xi,fe.eta);

          // Compute Jacobian inverse of coordinate mapping and derivatives
          fe.detJxW = 0.5*utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);
          if (fe.detJxW == 0.0) continue; // skip singular points

          // Cartesian coordinates of current integration point
          X = Xnod * fe.N;
          X.t = time.t;

          // Evaluate the integrand and accumulate element contributions
          fe.detJxW *= wg[j];
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


bool ASMs2DTri::integrate (Integrand& integrand, int lIndex,
                           GlobalIntegral& glInt,
                           const TimeDomain& time)
{
  if (!surf) return true; // silently ignore empty patches

  // Get Gaussian quadrature points and weights
  int nGP = integrand.getBouIntegrationPoints(2);
  const double* xg = GaussQuadrature::getCoord(nGP);
  const double* wg = GaussQuadrature::getWeight(nGP);
  if (!xg || !wg) return false;

  // Find the parametric direction of the edge normal {-2,-1, 1, 2}
  const int edgeDir = (lIndex%10+1)/(lIndex%2 ? -2 : 2);

  const int t1 = abs(edgeDir); // tangent direction normal to the patch edge

  // Order of basis in the two parametric directions (order = degree + 1)
  const int p1 = surf->order_u();
  const int p2 = surf->order_v();

  // Number of elements in each direction
  const int nelx = (nx-1)/(p1-1);
  const int nely = (ny-1)/(p2-1);
  bool UnionJack = nelx%2 == 0 && nely%2 == 0;

  // Get parametric coordinates of the grid points
  FiniteElement fe(p1*(p2+1)/2);
  RealArray upar, vpar;
  if (t1 == 1)
  {
    fe.u = edgeDir < 0 ? surf->startparam_u() : surf->endparam_u();
    this->getGridParameters(vpar,1,1);
  }
  else if (t1 == 2)
  {
    this->getGridParameters(upar,0,1);
    fe.v = edgeDir < 0 ? surf->startparam_v() : surf->endparam_v();
  }

  // Extract the Neumann order flag (1 or higher) for the integrand
  integrand.setNeumannOrder(1 + lIndex/10);

  // Integrate the extraordinary elements?
  size_t doXelms = 0;
  if (integrand.getIntegrandType() & Integrand::XO_ELEMENTS)
    if ((doXelms = nelx*nely*2)*2 > MNPC.size())
    {
      std::cerr <<" *** ASMs2DTri::integrate: Too few XO-elements "
                << MNPC.size() - doXelms << std::endl;
      return false;
    }

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex%10);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;

  Matrix dNdu(fe.N.size(),2), Xnod, Jac;
  Vec4   X;
  Vec3   normal, v1, v2;


  // === Assembly loop over all elements on the patch edge =====================

  int iel = 0;
  for (int i2 = 0; i2 < nely; i2++)
    for (int i1 = 0; i1 < nelx; i1++, iel += 2)
    {
      // Skip elements that are not on current boundary edge
      bool flipDiag = UnionJack && i1%2 + i2%2 == 1;
      bool skipMe = false;
      int jel = iel;
      switch (edgeDir)
        {
        case -1:
          if (i1 > 0)
            skipMe = true;
          else if (!flipDiag)
            ++jel;
          break;
        case  1:
          if (i1 < nelx-1)
            skipMe = true;
          else if (flipDiag)
            ++jel;
          break;
        case -2:
          if (i2 > 0)
            skipMe = true;
          break;
        case  2:
          if (i2 < nely-1)
            skipMe = true;
          ++jel;
          break;
        }
      if (skipMe) continue;

      // Set up nodal point coordinates for current element
      if (!this->getElementCoordinates(Xnod,1+jel)) return false;

      // Initialize element quantities
      fe.iel = abs(MLGE[doXelms+jel]);
      LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel,true);
      bool ok = integrand.initElementBou(MNPC[doXelms+jel],*A);


      // --- Integration loop over all Gauss points along the edge -------------

      int jp = (t1 == 1 ? i2 : i1)*nGP;
      fe.iGP = firstp + jp; // Global integration point counter

      for (int i = 0; i < nGP && ok; i++, fe.iGP++)
      {
        // Convert from natural coordinates in bi-unit square [-1,1]x[-1,1]
        // to area coordinates [0,1] of the current triangle
        if (t1 == 2)
        {
          // On upper or lower edge; node 1 and 2 are in action
          fe.xi = 0.5-0.5*xg[i];
          fe.eta = 1.0 - fe.xi;
          if (edgeDir > 0) // upper edge
            std::swap(fe.xi,fe.eta);
        }
        else if (flipDiag)
        {
          // On left or right edge and flipped diagonal; node 1 and 3
          fe.xi = 0.5-0.5*xg[i];
          fe.eta = 0.0;
          if (edgeDir > 0) // right edge
            fe.xi = 1.0 - fe.xi;
        }
        else
        {
          // On left or right edge but not flipped diagonal; node 2 and 3
          fe.eta = 0.5-0.5*xg[i];
          fe.xi = 0.0;
          if (edgeDir < 0) // left edge
            fe.eta = 1.0 - fe.eta;
        }

        // Parameter values of current integration point
        if (upar.size() > 1)
          fe.u = 0.5*(upar[i1]*(1.0-xg[i]) + upar[i1+1]*(1.0+xg[i]));
        if (vpar.size() > 1)
          fe.v = 0.5*(vpar[i2]*(1.0-xg[i]) + vpar[i2+1]*(1.0+xg[i]));

        // Compute basis function values and derivatives at current point
        evalBasis(fe.N,dNdu,fe.xi,fe.eta);

        // Compute basis function derivatives
        utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

        // Compute the edge normal and the curve dilatation (dS = |v1|).
        // This calculation is valid for linear triangles only.
        if (t1 == 2)
        {
          v1 = Xnod.getColumn(2) - Xnod.getColumn(1);
          normal = Xnod.getColumn(3) - Xnod.getColumn(1);
        }
        else if (flipDiag)
        {
          v1 = Xnod.getColumn(1) - Xnod.getColumn(3);
          normal = Xnod.getColumn(2) - Xnod.getColumn(3);
        }
        else
        {
          v1 = Xnod.getColumn(3) - Xnod.getColumn(2);
          normal = Xnod.getColumn(1) - Xnod.getColumn(2);
        }
        v2.cross(v1,normal).normalize();
        fe.detJxW = 0.5*v1.normalize();
        normal.cross(v1,v2);

        // Cartesian coordinates of current integration point
        X = Xnod * fe.N;
        X.t = time.t;

        // Evaluate the integrand and accumulate element contributions
        fe.detJxW *= wg[i];
        if (ok && !integrand.evalBou(*A,fe,time,X,normal))
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


bool ASMs2DTri::tesselate (ElementBlock& grid, const int* npe) const
{
  if (npe[0] != 2 || npe[1] != 2)
  {
    int* newnpe = const_cast<int*>(npe);
    std::cout <<"\nTriangle elements: The number of visualization points are"
              <<" 2 2 by default\n"<< std::endl;
    newnpe[0] = newnpe[1] = 2;
  }

  // Establish the block grid coordinates
  size_t i, j, ip;
  grid.setNoElmNodes(3);
  grid.resize(nx,ny);
  for (i = 0; i < grid.getNoNodes(); i++)
    grid.setCoor(i,this->getCoord(1+i));

  // Establish the block grid topology
  for (i = ip = 0; i < MNPC.size(); i++)
    for (j = 0; j < 3; j++, ip++)
      grid.setNode(ip,MNPC[i][j]);

  return true;
}


bool ASMs2DTri::evalSolution (Matrix& sField, const IntegrandBase& integrand,
                             const RealArray*, bool) const
{
  sField.resize(0,0);

  FiniteElement fe(3);
  Vector        solPt;
  Vectors       globSolPt(nnod);
  IntVec        check(nnod,0);
  Matrix        dNdu(3,2), Xnod, Jac;

  // Evaluate the secondary solution field at each point
  const int nel = this->getNoElms(true);
  for (int iel = 1; iel <= nel; iel++)
  {
    const IntVec& mnpc = MNPC[iel-1];
    this->getElementCoordinates(Xnod,iel);

    for (int i = 0; i < 3; i++)
    {
      fe.xi = i == 0 ? 1.0 : 0.0;
      fe.eta = i == 1 ? 1.0 : 0.0;
      evalBasis(fe.N,dNdu,fe.xi,fe.eta);
      fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

      if (!integrand.evalSol(solPt,fe,Xnod*fe.N,mnpc))
        return false;
      else if (sField.empty())
        sField.resize(solPt.size(),nnod,true);

      if (++check[mnpc[i]] == 1)
        globSolPt[mnpc[i]] = solPt;
      else
        globSolPt[mnpc[i]] += solPt;
      }
  }

  for (size_t i = 0; i < nnod; i++)
    sField.fillColumn(1+i,globSolPt[i] /= check[i]);

  return true;
}


void ASMs2DTri::generateThreadGroups (const Integrand&, bool, bool)
{
  threadGroups.calcGroups(nx-1,ny-1,1);
#if defined(USE_OPENMP) && SP_DEBUG > 1
  std::cout <<"\nThreading groups after triangularization:"<< std::endl;
#endif
  for (size_t g = 0; g < threadGroups.size(); g++)
    for (size_t t = 0; t < threadGroups[g].size(); t++)
    {
      IntVec& tGroup = const_cast<IntMat&>(threadGroups[g])[t];
      size_t oldsize = tGroup.size();
      tGroup.resize(oldsize*2);
      for (int i = oldsize-1; i >= 0; i--)
      {
        int iel1 = 2*tGroup[i];
        tGroup[2*i]   = iel1;
        tGroup[2*i+1] = iel1+1;
      }
#if defined(USE_OPENMP) && SP_DEBUG > 1
      if (t == 0)
        std::cout <<"group "<< g << std::endl;
      std::cout <<"\tthread "<< t <<":";
      for (size_t k = 0; k < tGroup.size(); k++)
        std::cout <<" "<< tGroup[k];
      std::cout << std::endl;
#endif
    }
}
