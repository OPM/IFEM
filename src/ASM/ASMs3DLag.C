// $Id$
//==============================================================================
//!
//! \file ASMs3DLag.C
//!
//! \date Feb 10 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for assembly of structured 3D %Lagrange FE models.
//!
//==============================================================================

#include "GoTools/trivariate/SplineVolume.h"

#include "ASMs3DLag.h"
#include "Lagrange.h"
#include "SparseMatrix.h"
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlbL2projector.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "ElementBlock.h"
#include "Utilities.h"
#include <array>
#include <numeric>


ASMs3DLag::ASMs3DLag (unsigned char n_f) : ASMs3D(n_f)
{
  nx = ny = nz = 0;
  p1 = p2 = p3 = 0;
}


ASMs3DLag::ASMs3DLag (const ASMs3DLag& patch, unsigned char n_f)
  : ASMs3D(patch,n_f), ASMLagBase(patch,false)
{
  nx = patch.nx;
  ny = patch.ny;
  nz = patch.nz;
  p1 = patch.p1;
  p2 = patch.p2;
  p3 = patch.p3;
}


ASMs3DLag::ASMs3DLag (const ASMs3DLag& patch)
  : ASMs3D(patch), ASMLagBase(patch)
{
  nx = patch.nx;
  ny = patch.ny;
  nz = patch.nz;
  p1 = patch.p1;
  p2 = patch.p2;
  p3 = patch.p3;
}


void ASMs3DLag::clear (bool retainGeometry)
{
  myCoord.clear();
  nx = ny = nz = 0;
  p1 = p2 = p3 = 0;

  this->ASMs3D::clear(retainGeometry);
}


bool ASMs3DLag::addXElms (short int dim, short int item, size_t nXn,
                          IntVec& nodes)
{
  if (!this->addXNodes(dim,nXn,nodes))
    return false;
  else if (p1 < 2 || p2 < 2 || p3 < 2)
    return false;

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
  if (!projB) projB = svol;

  // Order of basis in the three parametric directions (order = degree + 1)
  p1 = svol->order(0);
  p2 = svol->order(1);
  p3 = svol->order(2);

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
  // Number of elements in patch
  nel = nelx*nely*nelz;

  // Global node numbers
  myMLGN.resize(nnod);
  std::iota(myMLGN.begin(), myMLGN.end(), gNod+1);
  gNod += nnod;

  // Evaluate the nodal coordinates in the physical space
  size_t i1, j1;
  myCoord.resize(nnod);
  RealArray XYZ(svol->dimension()*nnod);
  svol->gridEvaluator(gpar1,gpar2,gpar3,XYZ);
  for (i1 = j1 = 0; i1 < myCoord.size(); i1++)
    for (int d = 0; d < 3; d++, j1++)
      myCoord[i1][d] = XYZ[j1];

  // Global element numbers
  myMLGE.resize(nel);
  std::iota(myMLGE.begin(), myMLGE.end(), gEl+1);
  gEl += nel;

  // Connectivity array: local --> global node relation
  ASMs3DLag::createMNPC(nx,ny,nz,p1,p2,p3,myMNPC);

  return true;
}


Vec3 ASMs3DLag::getCoord (size_t inod) const
{
  if (inod < 1 || inod > coord.size())
    return Vec3();

  return coord[inod-1];
}


void ASMs3DLag::setCoord (size_t inod, const Vec3& Xnod)
{
  if (inod < 1)
    return;

  if (inod > nnod)
    myCoord.resize(nnod = inod);

  myCoord[inod-1] = Xnod;
}


Vec3 ASMs3DLag::getElementCenter (int iel) const
{
  if (iel < 1 || static_cast<size_t>(iel) > MNPC.size())
  {
    std::cerr <<" *** ASMs3DLag::getElementCenter: Element index "<< iel
              <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return Vec3();
  }

  return this->getGeometricCenter(MNPC[--iel]);
}


bool ASMs3DLag::getElementCoordinates (Matrix& X, int iel, bool) const
{
  if (iel < 1 || static_cast<size_t>(iel) > MNPC.size())
  {
    std::cerr <<" *** ASMs3DLag::getElementCoordinates: Element index "<< iel
              <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return false;
  }

  // Number of nodes per element
  size_t nen = std::min(static_cast<size_t>(p1*p2*p3),MNPC[--iel].size());

  X.resize(3,nen);
  for (size_t i = 0; i < nen; i++)
    X.fillColumn(i+1,coord[MNPC[iel][i]].ptr());

  return true;
}


void ASMs3DLag::getNodalCoordinates (Matrix& X, bool) const
{
  X.resize(3,coord.size());

  for (size_t inod = 0; inod < coord.size(); inod++)
    X.fillColumn(inod+1,coord[inod].ptr());
}


bool ASMs3DLag::updateCoords (const Vector& displ)
{
  return shareFE ? true : this->ASMLagBase::updateCoords(displ,3);
}


bool ASMs3DLag::getOrder (int& pu, int& pv, int& pw) const
{
  pu = p1;
  pv = p2;
  pw = p3;

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
  if (ldim < 1 && lIndex > 0)
    return 1;

  size_t nel[3] = { // Number of elements in each direction
    nx > 1 && p1 > 1 ? (nx-1)/(p1-1) : 0 ,
    ny > 1 && p2 > 1 ? (ny-1)/(p2-1) : 0 ,
    nz > 1 && p3 > 1 ? (nz-1)/(p3-1) : 0 };

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
  if (this->empty()) return true; // silently ignore empty patches

  if (myCache.empty())
    myCache.emplace_back(std::make_unique<BasisFunctionCache>(*this));

  ASMs3D::BasisFunctionCache& cache = *myCache.front();
  cache.setIntegrand(&integrand);
  cache.init(integrand.getIntegrandType() & Integrand::NO_DERIVATIVES ? 0 : 1);

  // Get Gaussian quadrature points and weights
  const std::array<int,3>& ng = cache.nGauss();
  const std::array<const double*,3>& xg = cache.coord();
  const std::array<const double*,3>& wg = cache.weight();

  // Get the reduced integration quadrature points, if needed
  const double* xr = cache.coord(true)[0];
  const double* wr = cache.weight(true)[0];

  // Number of elements in each direction
  const int nel1 = cache.noElms()[0];
  const int nel2 = cache.noElms()[1];


  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t g = 0; g < threadGroupsVol.size() && ok; g++)
  {
#pragma omp parallel for schedule(static)
    for (const IntVec& group : threadGroupsVol[g])
    {
      FiniteElement fe;
      Matrix Jac;
      Vec4 X(nullptr,time.t);
      for (size_t e = 0; e < group.size() && ok; e++)
      {
        int iel = group[e];
        if (iel < 0 || iel >= static_cast<int>(nel))
        {
          ok = false;
          break;
        }

        fe.iel = MLGE[iel];
        if (!this->isElementActive(fe.iel)) continue; // zero-volume element

        // Set up nodal point coordinates for current element
        this->getElementCoordinates(fe.Xn,1+iel);

        if (integrand.getIntegrandType() & Integrand::ELEMENT_CENTER)
        {
          // Compute the element "center" (average of element node coordinates)
          X = 0.0;
          for (size_t i = 1; i <= 3; i++)
            for (size_t j = 1; j <= fe.Xn.cols(); j++)
              X[i-1] += fe.Xn(i,j);

          X *= 1.0 / static_cast<double>(fe.Xn.cols());
        }

        // Initialize element quantities
        LocalIntegral* A = integrand.getLocalIntegral(fe.Xn.cols(),fe.iel);
        const int nRed = cache.nGauss(true).front();
        if (!integrand.initElement(MNPC[iel],fe,X,nRed*nRed*nRed,*A))
        {
          A->destruct();
          ok = false;
          break;
        }

        int i1 = nel1*nel2 > 0 ?  iel % nel1         : 0;
        int i2 = nel1*nel2 > 0 ? (iel / nel1) % nel2 : 0;
        int i3 = nel1*nel2 > 0 ?  iel / (nel1*nel2)  : 0;

        if (xr)
        {
          // --- Selective reduced integration loop ----------------------------

          size_t ip = 0;
          for (int k = 0; k < nRed; k++)
            for (int j = 0; j < nRed; j++)
              for (int i = 0; i < nRed; i++, ++ip)
              {
                // Local element coordinates of current integration point
                fe.xi   = xr[i];
                fe.eta  = xr[j];
                fe.zeta = xr[k];

                // Parameter value of current integration point
                fe.u = cache.getParam(0,i1,i,true);
                fe.v = cache.getParam(1,i2,j,true);
                fe.w = cache.getParam(2,i3,k,true);

                // Compute basis function derivatives at current point
                if (integrand.getIntegrandType() & Integrand::NO_DERIVATIVES)
                  fe.N = cache.getVals(iel,ip,true).N;
                else
                {
                  const BasisFunctionVals& bfs = cache.getVals(iel,ip,true);
                  fe.N = bfs.N;

                  // Compute Jacobian inverse and derivatives
                  fe.detJxW = utl::Jacobian(Jac,fe.dNdX,fe.Xn,bfs.dNdu);
                }

                // Cartesian coordinates of current integration point
                X.assign(fe.Xn * fe.N);

                // Compute the reduced integration terms of the integrand
                fe.detJxW *= wr[i]*wr[j]*wr[k];
                if (!integrand.reducedInt(*A,fe,X))
                  ok = false;
              }
        }


        // --- Integration loop over all Gauss points in each direction --------

        size_t ip = 0;
        int jp = iel*ng[0]*ng[1]*ng[2];
        fe.iGP = firstIp + jp; // Global integration point counter

        for (int k = 0; k < ng[2]; k++)
          for (int j = 0; j < ng[1]; j++)
            for (int i = 0; i < ng[0]; i++, fe.iGP++, ++ip)
            {
              // Local element coordinates of current integration point
              fe.xi   = xg[0][i];
              fe.eta  = xg[1][j];
              fe.zeta = xg[2][k];

              // Parameter value of current integration point
              fe.u = cache.getParam(0,i1,i);
              fe.v = cache.getParam(1,i2,j);
              fe.w = cache.getParam(2,i3,k);

              if (integrand.getIntegrandType() & Integrand::NO_DERIVATIVES)
                fe.N = cache.getVals(iel,ip).N;
              else
              {
                const BasisFunctionVals& bfs = cache.getVals(iel,ip);
                fe.N = bfs.N;

                // Compute Jacobian inverse of coordinate mapping + derivatives
                fe.detJxW = utl::Jacobian(Jac,fe.dNdX,fe.Xn,bfs.dNdu);
                if (fe.detJxW == 0.0) continue; // skip singular points
              }

              // Cartesian coordinates of current integration point
              X.assign(fe.Xn * fe.N);

              // Evaluate the integrand and accumulate element contributions
              fe.detJxW *= wg[0][i]*wg[1][j]*wg[2][k];
              if (!integrand.evalInt(*A,fe,time,X))
                ok = false;
            }

        // Finalize the element quantities
        if (ok && !integrand.finalizeElement(*A,fe,time,firstIp+jp))
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
  if (this->empty()) return true; // silently ignore empty patches

  std::map<char,ThreadGroups>::const_iterator tit;
  if ((tit = threadGroupsFace.find(lIndex%10)) == threadGroupsFace.end())
  {
    std::cerr <<" *** ASMs3DLag::integrate: No thread groups for face "
              << lIndex%10 << std::endl;
    return false;
  }
  const ThreadGroups& threadGrp = tit->second;

  // Find the parametric direction of the face normal {-3,-2,-1, 1, 2, 3}
  const int faceDir = (lIndex%10+1)/((lIndex%2) ? -2 : 2);

  const int t0 = abs(faceDir); // unsigned normal direction of the face
  const int t1 = 1 + t0%3; // first tangent direction of the face
  const int t2 = 1 + t1%3; // second tangent direction of the face

  // Get Gaussian quadrature points and weights
  // and compute parameter values of the Gauss points over the whole patch face
  std::array<int,3> ng;
  std::array<const double*,3> xg, wg;
  for (int d = 0; d < 3; d++)
    if (-1-d == faceDir)
    {
      ng[d] = 1;
      xg[d] = nullptr;
      wg[d] = nullptr;
    }
    else if (1+d == faceDir)
    {
      ng[d] = 1;
      xg[d] = nullptr;
      wg[d] = nullptr;
    }
    else
    {
      int n = svol ? this->getNoGaussPt(svol->order(d),true) : 0;
      ng[d] = integrand.getBouIntegrationPoints(n);
      xg[d] = GaussQuadrature::getCoord(ng[d]);
      wg[d] = GaussQuadrature::getWeight(ng[d]);
    }

  const int tt0 = t0-1;
  const int tt1 = t1 > t2 ? t2-1 : t1-1;
  const int tt2 = t1 > t2 ? t1-1 : t2-1;
  const int nG1 = ng[tt1];
  const int nG2 = ng[tt2];

  // Number of elements in each direction
  const int nel1 = (nx-1)/(p1-1);
  const int nel2 = (ny-1)/(p2-1);
  const int nel3 = (nz-1)/(p3-1);

  // Get parametric coordinates of the elements
  RealArray upar, vpar, wpar;
  if (svol)
  {
    if (t0 == 1)
      upar = { faceDir < 0 ? svol->startparam(0) : svol->endparam(0) };
    else if (t0 == 2)
      vpar = { faceDir < 0 ? svol->startparam(1) : svol->endparam(1) };
    else if (t0 == 3)
      wpar = { faceDir < 0 ? svol->startparam(2) : svol->endparam(2) };

    if (upar.empty()) this->getGridParameters(upar,0,1);
    if (vpar.empty()) this->getGridParameters(vpar,1,1);
    if (wpar.empty()) this->getGridParameters(wpar,2,1);
  }

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

  // Lambda function for linear interpolation between two values
  auto&& linearInt = [](double xi, double v1, double v2)
  {
    return 0.5*(v1*(1.0-xi) + v2*(1.0+xi));
  };


  // === Assembly loop over all elements on the patch face =====================

  bool ok = true;
  for (size_t g = 0; g < threadGrp.size() && ok; g++)
  {
#pragma omp parallel for schedule(static)
    for (const IntVec& group : threadGrp[g])
    {
      FiniteElement fe;
      fe.u = upar.front();
      fe.v = vpar.front();
      fe.w = wpar.front();

      Matrix dNdu, Jac;
      double param[3] = { fe.u, fe.v, fe.w };
      Vec4   X(param,time.t);
      Vec3   normal;
      double xi[3];

      for (size_t e = 0; e < group.size() && ok; e++)
      {
        int iel = group[e];
        if (iel < 0 || iel >= static_cast<int>(nel))
        {
          ok = false;
          break;
        }

        fe.iel = abs(MLGE[doXelms+iel]);
        if (!this->isElementActive(fe.iel)) continue; // zero-volume element

        // Set up nodal point coordinates for current element
        this->getElementCoordinates(fe.Xn,1+iel);

        // Initialize element quantities
        LocalIntegral* A = integrand.getLocalIntegral(fe.Xn.cols(),fe.iel,true);
        if (!integrand.initElementBou(MNPC[doXelms+iel],*A))
        {
          A->destruct();
          ok = false;
          break;
        }

        int i1 =  iel % nel1;
        int i2 = (iel / nel1) % nel2;
        int i3 =  iel / (nel1*nel2);

        // Define some loop control variables depending on which face we are on
        int nf1 = 0, j1 = 0, j2 = 0;
        switch (abs(faceDir))
        {
          case 1: nf1 = nel2; j2 = i3; j1 = i2; break;
          case 2: nf1 = nel1; j2 = i3; j1 = i1; break;
          case 3: nf1 = nel1; j2 = i2; j1 = i1; break;
        }


        // --- Integration loop over all Gauss points in each direction --------

        int jp = (j2*nf1 + j1)*nG1*nG2;
        fe.iGP = firstp + jp; // Global integration point counter

        for (int j = 0; j < nG2 && ok; j++)
          for (int i = 0; i < nG1 && ok; i++, fe.iGP++)
          {
            // Local element coordinates of current integration point
            xi[tt0] = faceDir < 0 ? -1.0 : 1.0;
            xi[tt1] = xg[tt1][i];
            xi[tt2] = xg[tt2][j];
            fe.xi   = xi[0];
            fe.eta  = xi[1];
            fe.zeta = xi[2];

            // Local element coordinates and parameter values
            // of current integration point
            int k1 = -1, k2 = -1, k3 = -1;
            switch (abs(faceDir))
            {
              case 1: k2 = i; k3 = j; break;
              case 2: k1 = i; k3 = j; break;
              case 3: k1 = i; k2 = j; break;
            }
            if (xg[0]) fe.u = linearInt(xg[0][k1],upar[i1],upar[i1+1]);
            if (xg[1]) fe.v = linearInt(xg[1][k2],vpar[i2],vpar[i2+1]);
            if (xg[2]) fe.w = linearInt(xg[2][k3],wpar[i3],wpar[i3+1]);

            // Compute the basis functions and their derivatives, using
            // tensor product of one-dimensional Lagrange polynomials
            if (integrand.getIntegrandType() & Integrand::NO_DERIVATIVES)
              ok = Lagrange::computeBasis(fe.N,p1,xi[0],p2,xi[1],p3,xi[2]);
            else
            {
              ok = Lagrange::computeBasis(fe.N,dNdu,p1,xi[0],p2,xi[1],p3,xi[2]);

              // Compute basis function derivatives and the face normal
              fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,fe.Xn,dNdu,t1,t2);
              if (fe.detJxW == 0.0) continue; // skip singular points

              if (faceDir < 0) normal *= -1.0;
            }

            // Cartesian coordinates of current integration point
            X.assign(fe.Xn * fe.N);

            // Evaluate the integrand and accumulate element contributions
            fe.detJxW *= wg[tt1][i]*wg[tt2][j];
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
  if (this->empty()) return true; // silently ignore empty patches

  // Parametric direction of the edge {0, 1, 2}
  const int lDir = (lEdge-1)/4;

  // Get Gaussian quadrature points and weights
  int ng = svol ? this->getNoGaussPt(svol->order(lDir),true) : 0;
  const double* xg = GaussQuadrature::getCoord(ng);
  const double* wg = GaussQuadrature::getWeight(ng);
  if (!xg || !wg) return false;

  // Number of elements in each direction
  const int nelx = (nx-1)/(p1-1);
  const int nely = (ny-1)/(p2-1);
  const int nelz = (nz-1)/(p3-1);

  FiniteElement fe;
  Matrix dNdu, Jac;
  Vec4   X(nullptr,time.t);
  Vec3   tangent;
  double xi[3] = {0.0, 0.0, 0.0};

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
        fe.iel = MLGE[iel-1];
        if (!this->isElementActive(fe.iel)) continue; // zero-volume element

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
          ip = i1*ng;
        else if (lEdge < 9)
          ip = i2*ng;
        else
          ip = i3*ng;

        // Set up nodal point coordinates for current element
        this->getElementCoordinates(fe.Xn,iel);

        // Initialize element quantities
        LocalIntegral* A = integrand.getLocalIntegral(fe.Xn.cols(),fe.iel,true);
        bool ok = integrand.initElementBou(MNPC[iel-1],*A);


        // --- Integration loop over all Gauss points along the edge -----------

        fe.iGP = firstp + ip; // Global integration point counter

        for (int i = 0; i < ng && ok; i++, fe.iGP++)
        {
          // Gauss point coordinates on the edge
          xi[lDir] = xg[i];

          // Compute the basis functions and their derivatives, using
          // tensor product of one-dimensional Lagrange polynomials
          if (integrand.getIntegrandType() & Integrand::NO_DERIVATIVES)
            ok = Lagrange::computeBasis(fe.N,p1,xi[0],p2,xi[1],p3,xi[2]);
          else
          {
            ok = Lagrange::computeBasis(fe.N,dNdu,p1,xi[0],p2,xi[1],p3,xi[2]);

            // Compute basis function derivatives and the edge tangent
            fe.detJxW = utl::Jacobian(Jac,tangent,fe.dNdX,fe.Xn,dNdu,1+lDir);
            if (fe.detJxW == 0.0) continue; // skip singular points
          }

          // Cartesian coordinates of current integration point
          X.assign(fe.Xn * fe.N);

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
  // Evaluate the parametric values of the point and nodes
  std::array<RealArray,3> u;
  for (int d = 0; d < 3; d++)
  {
    if (svol)
      param[d] = (1.0-xi[d])*svol->startparam(d) + xi[d]*svol->endparam(d);
    else
      param[d] = xi[d];
    if (svol && !this->getGridParameters(u[d],d,svol->order(d)-1)) return -3;
  }

  // Search for the closest node
  size_t i = utl::find_closest(u[0],param[0]);
  size_t j = utl::find_closest(u[1],param[1]);
  size_t k = utl::find_closest(u[2],param[2]);
  size_t n = u[0].size()*(u[1].size()*k + j) + i;
  X = coord[n];

  return 1+n;
}


bool ASMs3DLag::tesselate (ElementBlock& grid, const int*) const
{
  const int npe[3] = { p1, p2, p3 };
  if (!this->ASMs3D::tesselate(grid,npe))
    return false;

  // Adjust the element Id since each Lagrange element covers several knot-spans
  int ie, nse1 = p1-1;
  int je, nse2 = p2-1;
  int ke, nse3 = p3-1;
  int nelx = (nx-1)/nse1;
  int nely = (ny-1)/nse1;
  for (size_t k = ke = 1; k < nz; k++)
  {
    for (size_t j = je = 1; j < ny; j++)
    {
      for (size_t i = ie = 1; i < nx; i++)
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
                              const int*, int n_f, bool) const
{
  return this->evalSolution(sField,locSol,nullptr,false,0,n_f);
}


int ASMs3DLag::findElement (double u, double v, double w,
                            double* xi, double* eta, double* zeta) const
{
  if (!svol)
  {
    std::cerr <<" *** ASMs3DLag::findElement: No spline geometry"<< std::endl;
    return -1;
  }

  const std::array<std::pair<double,int>,3> knot {{
    {u, svol->basis(0).knotInterval(u)},
    {v, svol->basis(1).knotInterval(v)},
    {w, svol->basis(2).knotInterval(w)}
  }};

  const std::array<int,3> elm {
    knot[0].second - (p1 - 1),
    knot[1].second - (p2 - 1),
    knot[2].second - (p3 - 1)
  };

  const std::array<int,2> nel {
    svol->numCoefs(0) - svol->order(0) + 1,
    svol->numCoefs(1) - svol->order(1) + 1
  };

  auto getParam = [this,&knot](int dir)
  {
    const double knot_1 = *(svol->basis(dir).begin() + knot[dir].second);
    const double knot_2 = *(svol->basis(dir).begin() + knot[dir].second + 1);
    return -1.0 + 2.0 * (knot[dir].first - knot_1) / (knot_2 - knot_1);
  };

  if (xi)
    *xi = getParam(0);
  if (eta)
    *eta = getParam(1);
  if (zeta)
    *zeta = getParam(2);

  return 1 + elm[0] + (elm[1] + elm[2] * nel[1]) * nel[0];
}


bool ASMs3DLag::evalSolution (Matrix& sField, const Vector& locSol,
                              const RealArray* gpar, bool regular,
                              int, int) const
{
  if (!gpar && !regular) // Direct nodal evaluation
    return this->nodalField(sField,locSol,this->getNoNodes(-1));

  size_t nCmp = locSol.size() / this->getNoProjectionNodes();
  size_t ni   = gpar ? gpar[0].size() : nel;
  size_t nj   = gpar ? gpar[1].size() : 1;
  size_t nk   = gpar ? gpar[2].size() : 1;
  size_t nen  = p1*p2*p3;

  sField.resize(nCmp,ni*nj*nk);
  Matrix elmSol(nCmp,nen);
  RealArray N(nen), val;

  double xi = 0.0, eta = 0.0, zeta = 0.0;
  if (!gpar && !Lagrange::computeBasis(N,p1,xi,p2,eta,p3,zeta))
    return false;

  size_t ip = 1;
  int iel = 0;
  for (size_t k = 0; k < nk; k++)
    for (size_t j = 0; j < nj; j++)
      for (size_t i = 0; i < ni; i++, ip++)
      {
        if (gpar)
        {
          iel = this->findElement(gpar[0][i], gpar[1][j], gpar[2][k],
                                  &xi, &eta, &zeta);
          if (iel < 1 || iel > static_cast<int>(nel))
            return false;
          if (!Lagrange::computeBasis(N,p1,xi,p2,eta,p3,zeta))
            return false;
        }
        else
          iel++;

        for (size_t a = 1; a <= nen; a++)
          elmSol.fillColumn(a, locSol.ptr() + nCmp*MNPC[iel-1][a-1]);

        elmSol.multiply(N,val);
        sField.fillColumn(ip,val);
      }

  return true;
}


bool ASMs3DLag::evalSolution (Matrix& sField, const IntegrandBase& integrand,
                              const int*, char) const
{
  sField.resize(0,0);

  double incx = 2.0/double(p1-1);
  double incy = 2.0/double(p2-1);
  double incz = 2.0/double(p3-1);

  size_t nPoints = coord.size();
  IntVec check(nPoints,0);

  FiniteElement fe(p1*p2*p3);
  Vector        solPt;
  Vectors       globSolPt(nPoints);
  Matrix        dNdu, Jac;

  // Evaluate the secondary solution field at each point
  for (size_t iel = 0; iel < nel; iel++)
  {
    fe.iel = MLGE[iel];
    if (fe.iel < 1) continue; // zero-volume element

    const IntVec& mnpc = MNPC[iel];
    this->getElementCoordinates(fe.Xn,1+iel);

    int i, j, k, loc = 0;
    for (k = 0; k < p3; k++)
      for (j = 0; j < p2; j++)
        for (i = 0; i < p1; i++, loc++)
        {
          fe.xi   = -1.0 + i*incx;
          fe.eta  = -1.0 + j*incy;
          fe.zeta = -1.0 + k*incz;
          if (integrand.getIntegrandType() & Integrand::NO_DERIVATIVES)
            Lagrange::computeBasis(fe.N,p1,fe.xi,p2,fe.eta,p3,fe.zeta);
          else
          {
            Lagrange::computeBasis(fe.N,dNdu,p1,fe.xi,p2,fe.eta,p3,fe.zeta);

            // Compute the Jacobian inverse
            fe.detJxW = utl::Jacobian(Jac,fe.dNdX,fe.Xn,dNdu);
            if (fe.detJxW == 0.0) continue; // skip singular points
          }

          // Now evaluate the solution field
          utl::Point X4(fe.Xn*fe.N,{fe.u,fe.v,fe.w});
          if (!integrand.evalSol(solPt,fe,X4,mnpc))
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
    if (check[i] == 1)
      sField.fillColumn(1+i, globSolPt[i]);
    else if (check[i] > 1)
      sField.fillColumn(1+i, globSolPt[i] /= check[i]);

  return true;
}


bool ASMs3DLag::evalSolution (Matrix& sField, const IntegrandBase& integrand,
                              const RealArray* gpar, bool regular) const
{
  sField.resize(0,0);

  bool elCenters = gpar ? regular && gpar->empty() : true;
  size_t nPoints = elCenters ? nel : gpar->size();

  FiniteElement fe(p1*p2*p3);
  Vector        solPt;
  Vectors       globSolPt(nPoints);
  Matrix        dNdu, Jac;

  // Evaluate the secondary solution field at each point or element center
  int iel = 0;
  for (size_t i = 0; i < nPoints; i++)
  {
    if (elCenters)
      iel++;
    else
    {
      iel = this->findElement(gpar[0][i], gpar[1][i], gpar[2][i],
                              &fe.xi, &fe.eta, &fe.zeta);
      if (iel < 1 || iel > static_cast<int>(nel))
        return false;
    }

    fe.iel = MLGE[iel-1];
    if (fe.iel < 1) continue; // zero-volume element

    this->getElementCoordinates(fe.Xn,iel);

    if (!Lagrange::computeBasis(fe.N,dNdu,p1,fe.xi,p2,fe.eta,p3,fe.zeta))
      return false;

    // Compute the Jacobian inverse
    fe.detJxW = utl::Jacobian(Jac,fe.dNdX,fe.Xn,dNdu);
    if (fe.detJxW == 0.0) continue; // skip singular points

    // Now evaluate the solution field
    utl::Point X4(fe.Xn*fe.N,{fe.u,fe.v,fe.w});
    if (!integrand.evalSol(solPt,fe,X4,MNPC[iel-1]))
      return false;
    else if (sField.empty())
      sField.resize(solPt.size(),nPoints,true);

    sField.fillColumn(1+i, solPt);
  }

  return true;
}


void ASMs3DLag::generateThreadGroups (const Integrand&, bool, bool)
{
  threadGroupsVol.calcGroups((nx-1)/(p1-1),(ny-1)/(p2-1),(nz-1)/(p3-1),1);
}


void ASMs3DLag::generateThreadGroups (char lIndex, bool, bool)
{
  if (threadGroupsFace.find(lIndex) != threadGroupsFace.end()) return;

  ThreadGroups& fGrp = threadGroupsFace[lIndex];
  switch (lIndex)
  {
    case 1:
    case 2:
      fGrp.calcGroups((ny-1)/(p2-1), (nz-1)/(p3-1), 1);
      break;
    case 3:
    case 4:
      fGrp.calcGroups((nx-1)/(p1-1), (nz-1)/(p3-1), 1);
      break;
    default:
      fGrp.calcGroups((nx-1)/(p1-1), (ny-1)/(p2-1), 1);
  }

  // Find elements that are on the boundary face 'lIndex'
  IntVec map;
  this->findBoundaryElms(map,lIndex);

  fGrp.applyMap(map);
}


IntMat ASMs3DLag::getElmNodes (int basis) const
{
  const Go::SplineVolume* sv = this->getBasis(basis);

  // Find number of nodes in each parameter direction,
  // accounting for possible zero-length knot-spans
  // (see ASMs3D::getGridParameters())
  int nx[3] = { 0, 0, 3 };
  RealArray::const_iterator uit, uend;
  for (int d = 0; d < 3; d++)
  {
    uit  = sv->basis(d).begin() + sv->basis(d).order()-1;
    uend = sv->basis(d).begin() + sv->basis(d).numCoefs()+1;
    for (double uprev = *(uit++); uit != uend; ++uit)
    {
      if (*uit > uprev)
        nx[d] += sv->basis(d).order()-1;
      uprev = *uit;
    }
    if (sv->basis(d).order() > 2)
      nx[d]++;
  }

  IntMat result;
  ASMs3DLag::createMNPC(nx[0],nx[1],nx[2],
                        sv->order(0),sv->order(1),sv->order(2), result);
  return result;
}


void ASMs3DLag::createMNPC (size_t nx, size_t ny, size_t nz,
                            int p1, int p2, int p3, IntMat& result)
{
  const int nelx = (nx-1) / (p1-1);
  const int nely = (ny-1) / (p2-1);
  const int nelz = (nz-1) / (p3-1);
  result.resize(nelx*nely*nelz);

  // Number of nodes per element
  const int nen = p1*p2*p3;
  // Number of nodes in a xy-surface of an element
  const int ct  = p1*p2;

  int iel = 0;
  for (int k = 0; k < nelz; k++)
    for (int j = 0; j < nely; j++)
      for (int i = 0; i < nelx; i++, iel++)
      {
        result[iel].resize(nen);
        // First node in current element
        int corner = (p3-1)*(nx*ny)*k + (p2-1)*nx*j + (p1-1)*i;

        for (int c = 0; c < p3; c++)
        {
          int cornod = ct*c;
          result[iel][cornod] = corner + c*nx*ny;
          for (int b = 1; b < p2; b++)
          {
            int facenod = cornod + b*p1;
            result[iel][facenod] = result[iel][cornod] + b*nx;
            for (int a = 1; a < p1; a++)
            {
              result[iel][facenod+a] = result[iel][facenod] + a;
              result[iel][cornod+a]  = result[iel][cornod] + a;
            }
          }
        }
      }
}


void ASMs3DLag::findBoundaryElms (IntVec& elms, int lIndex, int) const
{
  const int N1m = (nx-1)/(p1-1);
  const int N2m = (ny-1)/(p2-1);
  const int N3m = (nz-1)/(p3-1);

  elms.clear();
  switch (lIndex) {
  case 1:
  case 2:
    elms.reserve(N2m*N3m);
    for (int k = 0; k < N3m; ++k)
      for (int j = 0; j < N2m; ++j)
        elms.push_back(j*N1m + k*N1m*N2m + (lIndex-1)*(N1m-1));
    break;
  case 3:
  case 4:
    elms.reserve(N1m*N3m);
    for (int k = 0; k < N3m; ++k)
      for (int i = 0; i < N1m; ++i)
        elms.push_back(i + k*N1m*N2m + (lIndex-3)*(N1m*(N2m-1)));
    break;
  case 5:
  case 6:
    elms.reserve(N1m*N2m);
    for (int j = 0; j < N2m; ++j)
      for (int i = 0; i < N1m; ++i)
        elms.push_back(i + j*N1m + (lIndex-5)*(N1m*N2m*(N3m-1)));
  }
}


bool ASMs3DLag::getGridParameters (RealArray& prm, int dir, int nSegPerSpan) const
{
  if (svol)
    return this->ASMs3D::getGridParameters(prm, dir, nSegPerSpan);

  if (nSegPerSpan < 1)
  {
    std::cerr <<" *** ASMs3DLag::getGridParameters: Too few knot-span points "
              << nSegPerSpan+1 <<" in direction "<< dir << std::endl;
    return false;
  }

  int nel1 = dir == 0 ? (nx-1)/(p1-1) : (dir == 1 ? (ny-1)/(p2-1) : (nz-1)/(p3-1));

  double ucurr = 0.0, uprev = 0.0, du = 1.0 / nel;
  for (int i = 0; i < nel1; ++i)
  {
    ucurr += du;
    if (ucurr > uprev)
      if (nSegPerSpan == 1)
        prm.push_back(uprev);
      else for (int j = 0; j < nSegPerSpan; j++)
      {
        double xg = (double)(2*j-nSegPerSpan)/(double)nSegPerSpan;
        prm.push_back(0.5*(ucurr*(1.0+xg) + uprev*(1.0-xg)));
      }
    uprev = ucurr;
  }

  if (ucurr > prm.back())
    prm.push_back(ucurr);

  return true;
}


bool ASMs3DLag::write (std::ostream& os, int) const
{
  return this->writeLagBasis(os,"hexahedron");
}


bool ASMs3DLag::evaluate (const ASMbase* basis, const Vector& locVec,
                          RealArray& vec, int basisNum) const
{
  // Compute parameter values of the result sampling points
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++) {
    int nel1 = dir == 0 ? (nx-1)/(p1-1) : (dir == 1 ? (ny-1)/(p2-1) : (nz-1)/(p3-1));
    gpar[dir].resize(nel1+1);
    double du = 1.0 / nel1;
    for (int i = 0; i <= nel1; ++i)
      gpar[dir][i] = i*du;
  }

  // Evaluate the result field at all sampling points.
  // Note: it is here assumed that *basis and *this have spline bases
  // defined over the same parameter domain.
  Matrix sValues;
  if (!basis->evalSolution(sValues,locVec,gpar.data()))
    return false;

  vec = sValues;
  return true;
}


void ASMs3DLag::constrainFace (int dir, bool open, int dof,
                               int code, char basis)
{
  this->ASMs3D::constrainFace(dir, open, dof, code > 0 ? -code : code, basis);
}


void ASMs3DLag::constrainEdge (int lEdge, bool open, int dof,
                               int code, char basis)
{
  this->ASMs3D::constrainEdge(lEdge, open, dof, code > 0 ? -code : code, basis);
}


size_t ASMs3DLag::getNoProjectionNodes () const
{
  return this->getNoNodes(1);
}


bool ASMs3DLag::assembleL2matrices (SystemMatrix& A, SystemVector& B,
                                    const L2Integrand& integrand,
                                    bool continuous) const
{
  BasisFunctionCache& cache = static_cast<BasisFunctionCache&>(*myCache.front());
  cache.init(1);

  const std::array<const double*,3>& wg = cache.weight();
  const std::array<int,3> nGP = cache.nGauss();

  Matrix dNdX, Xnod, J;

  const size_t nnod = this->getNoProjectionNodes();

  // === Assembly loop over all elements in the patch ==========================

  int iel = 0;
  for (size_t i3 = 0; i3 < cache.noElms()[2]; ++i3)
    for (size_t i2 = 0; i2 < cache.noElms()[1]; ++i2)
      for (size_t i1 = 0; i1 < cache.noElms()[0]; ++i1, ++iel)
      {
        if (!this->getElementCoordinates(Xnod,1+iel))
          return false;

        std::array<RealArray,3> GP;
        GP[0].reserve(nGP[0]*nGP[1]*nGP[2]);
        GP[1].reserve(nGP[0]*nGP[1]*nGP[2]);
        GP[2].reserve(nGP[0]*nGP[1]*nGP[2]);
        for (int k = 0; k < nGP[2]; ++k)
          for (int j = 0; j < nGP[1]; ++j)
            for (int i = 0; i < nGP[0]; ++i) {
              GP[0].push_back(cache.getParam(0, i1, i));
              GP[1].push_back(cache.getParam(1, i2, j));
              GP[2].push_back(cache.getParam(2, i3, k));
            }

        Matrix sField;
        integrand.evaluate(sField, GP.data());

        // --- Integration loop over all Gauss points in each direction --------

        Matrix eA(p1*p2*p3, p1*p2*p3);
        Vectors eB(sField.rows(), Vector(p1*p2*p3));
        size_t ip = 0;
        for (int k = 0; k < nGP[2]; ++k)
          for (int j = 0; j < nGP[1]; ++j)
            for (int i = 0; i < nGP[0]; ++i, ++ip) {
              const BasisFunctionVals& bfs = cache.getVals(iel,ip);

              double dJw = (wg[0][i]*wg[1][j]*wg[2][k] *
                            utl::Jacobian(J,dNdX,Xnod,bfs.dNdu));

              // Integrate the mass matrix
              eA.outer_product(bfs.N, bfs.N, true, dJw);

              // Integrate the rhs vector B
              for (size_t r = 1; r <= sField.rows(); r++)
                eB[r-1].add(bfs.N,sField(r,ip+1)*dJw);
            }

        const IntVec& mnpc = MNPC[iel];
        A.assemble(eA, mnpc);
        B.assemble(eB, mnpc, nnod);
      }

  return true;
}


ASMs3DLag::BasisFunctionCache::BasisFunctionCache (const ASMs3DLag& pch)
  : ASMs3D::BasisFunctionCache(pch)
{
  nel[0] = (pch.nx-1) / (pch.p1-1);
  nel[1] = (pch.ny-1) / (pch.p2-1);
  nel[2] = (pch.nz-1) / (pch.p3-1);
}


ASMs3DLag::BasisFunctionCache::BasisFunctionCache (const ASMs3D::BasisFunctionCache& cache,
                                                   int b) :
  ASMs3D::BasisFunctionCache(cache,b)
{
}


bool ASMs3DLag::BasisFunctionCache::internalInit ()
{
  if (!mainQ->xg[0])
    this->setupQuadrature();

  nTotal = mainQ->ng[0]*mainQ->ng[1]*mainQ->ng[2];
  if (reducedQ->xg[0])
    nTotalRed = reducedQ->ng[0]*reducedQ->ng[1]*reducedQ->ng[2];

  return true;
}


void ASMs3DLag::BasisFunctionCache::setupParameters ()
{
  // Compute parameter values of the Gauss points over the whole patch
  for (int d = 0; d < 3; d++) {
    patch.getGridParameters(mainQ->gpar[d],d,1);
    if (reducedQ->xg[0])
      reducedQ->gpar[d] = mainQ->gpar[d];
  }
}


double ASMs3DLag::BasisFunctionCache::getParam (int dir, size_t el,
                                                size_t gp, bool reduced) const
{
  const Quadrature& q = reduced ? *reducedQ : *mainQ;
  return 0.5*(q.gpar[dir][el  ]*(1.0-q.xg[dir][gp]) +
              q.gpar[dir][el+1]*(1.0+q.xg[dir][gp]));
}


BasisFunctionVals ASMs3DLag::BasisFunctionCache::calculatePt (size_t, size_t gp,
                                                              bool red) const
{
  const Quadrature& q = red ? *reducedQ : *mainQ;
  std::array<size_t,3> gpIdx = this->gpIndex(gp,red);
  const size_t i = gpIdx[0];
  const size_t j = gpIdx[1];
  const size_t k = gpIdx[2];

  int p1, p2, p3;
  patch.getOrder(p1,p2,p3);

  BasisFunctionVals result;
  if (integrand->getIntegrandType() & Integrand::NO_DERIVATIVES)
    Lagrange::computeBasis(result.N,
                           p1,q.xg[0][i], p2,q.xg[1][j], p3,q.xg[2][k]);
  else
    Lagrange::computeBasis(result.N, result.dNdu,
                           p1,q.xg[0][i], p2,q.xg[1][j], p3,q.xg[2][k]);

  return result;
}


void ASMs3DLag::BasisFunctionCache::calculateAll ()
{
  int i, j, k, p1, p2, p3;
  patch.getOrder(p1,p2,p3);
  bool noDerivs = integrand->getIntegrandType() & Integrand::NO_DERIVATIVES;

  // Evaluate basis function values and derivatives at all integration points.
  // They will be the same for all elements in the patch,
  // so no need to repeat for all.
  std::vector<BasisFunctionVals>::iterator it = values.begin();
  for (k = 0; k < mainQ->ng[2]; k++)
    for (j = 0; j < mainQ->ng[1]; j++)
      for (i = 0; i < mainQ->ng[0]; i++, ++it)
        if (noDerivs)
          Lagrange::computeBasis(it->N,
                                 p1, mainQ->xg[0][i],
                                 p2, mainQ->xg[1][j],
                                 p3, mainQ->xg[2][k]);
        else
          Lagrange::computeBasis(it->N, it->dNdu,
                                 p1, mainQ->xg[0][i],
                                 p2, mainQ->xg[1][j],
                                 p3, mainQ->xg[2][k]);

  if (reducedQ->xg[0])
    for (k = 0, it = valuesRed.begin(); k < reducedQ->ng[2]; k++)
      for (j = 0; j < reducedQ->ng[1]; j++)
        for (i = 0; i < reducedQ->ng[0]; i++, ++it)
          if (noDerivs)
            Lagrange::computeBasis(it->N,
                                   p1, reducedQ->xg[0][i],
                                   p2, reducedQ->xg[1][j],
                                   p3, reducedQ->xg[2][k]);
          else
            Lagrange::computeBasis(it->N, it->dNdu,
                                   p1, reducedQ->xg[0][i],
                                   p2, reducedQ->xg[1][j],
                                   p3, reducedQ->xg[2][k]);
}
