// $Id$
//==============================================================================
//!
//! \file ASMs2DLag.C
//!
//! \date Mar 20 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for assembly of structured 2D %Lagrange FE models.
//!
//==============================================================================

#include "GoTools/geometry/SplineSurface.h"

#include "ASMs2DLag.h"
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


ASMs2DLag::ASMs2DLag (unsigned char n_s, unsigned char n_f) : ASMs2D(n_s,n_f)
{
  nx = ny = 0;
  p1 = p2 = 0;
}


ASMs2DLag::ASMs2DLag (const ASMs2DLag& patch, unsigned char n_f)
  : ASMs2D(patch,n_f), ASMLagBase(patch,false)
{
  nx = patch.nx;
  ny = patch.ny;
  p1 = patch.p1;
  p2 = patch.p2;
}


ASMs2DLag::ASMs2DLag (const ASMs2DLag& patch)
  : ASMs2D(patch), ASMLagBase(patch)
{
  nx = patch.nx;
  ny = patch.ny;
  p1 = patch.p1;
  p2 = patch.p2;
}


void ASMs2DLag::clear (bool retainGeometry)
{
  myCoord.clear();
  nx = ny = 0;
  p1 = p2 = 0;

  this->ASMs2D::clear(retainGeometry);
}


bool ASMs2DLag::addXElms (short int dim, short int item, size_t nXn,
                          IntVec& nodes)
{
  if (!this->addXNodes(dim,nXn,nodes))
    return false;
  else if (p1 < 2 || p2 < 2)
    return false;

  const int nelx = (nx-1)/(p1-1);
  const int nely = (ny-1)/(p2-1);

  int iel = 0;
  bool skipMe = false;
  for (int i2 = 0; i2 < nely; i2++)
    for (int i1 = 0; i1 < nelx; i1++, iel++)
    {
      if (MLGE[iel] < 1) continue; // Skip zero-area element

      // Skip elements that are not on current boundary edge
      switch (item)
      {
        case 1: skipMe = i1 > 0;      break;
        case 2: skipMe = i1 < nelx-1; break;
        case 3: skipMe = i2 > 0;      break;
        case 4: skipMe = i2 < nely-1; break;
      }
      if (skipMe) continue;

      IntVec& mnpc = myMNPC[nel+iel];
      if (!mnpc.empty())
      {
        std::cerr <<" *** ASMs2DLag::addXElms: Only one X-edge allowed."
                  << std::endl;
        return false;
      }

      mnpc = MNPC[iel]; // Copy the ordinary element nodes

      // Negate node numbers that are not on the boundary edge, to flag that
      // they shall not receive any tangent and/or residual contributions
      int lnod = 0;
      for (int j2 = 0; j2 < p2; j2++)
        for (int j1 = 0; j1 < p1; j1++, lnod++)
        {
          switch (item)
          {
            case 1: skipMe = j1 > 0;    break;
            case 2: skipMe = j1 < p1-1; break;
            case 3: skipMe = j2 > 0;    break;
            case 4: skipMe = j2 < p2-1; break;
          }
          if (skipMe) // Hack for node 0: Using -maxint as flag instead
            mnpc[lnod] = mnpc[lnod] == 0 ? -2147483648 : -mnpc[lnod];
        }

      // Add connectivity to the extra-ordinary nodes
      for (size_t i = 0; i < nXn; i++)
        mnpc.push_back(MLGN.size()-nXn+i);

      myMLGE[nel+iel] = -(++gEl); // Flag extraordinary element by negative sign
    }

  return true;
}


bool ASMs2DLag::generateFEMTopology ()
{
  if (!surf) return false;
  if (!projB) projB = surf;

  // Order of basis in the two parametric directions (order = degree + 1)
  p1 = surf->order_u();
  p2 = surf->order_v();

  // Evaluate the parametric values
  RealArray gpar1, gpar2;
  if (!this->getGridParameters(gpar1,0,p1-1)) return false;
  if (!this->getGridParameters(gpar2,1,p2-1)) return false;

  // Number of nodes in each direction
  nx = gpar1.size();
  ny = gpar2.size();
  // Number of nodes in the patch
  nnod = nx*ny;

  if (!coord.empty())
    return coord.size() == nnod;

  // Number of elements in each direction
  const int nelx = (nx-1)/(p1-1);
  const int nely = (ny-1)/(p2-1);

  // Evaluate the nodal coordinates in the physical space
  RealArray XYZ(surf->dimension()*nnod);
  surf->gridEvaluator(XYZ,gpar1,gpar2);

  size_t i1, j1;
  myMLGN.resize(nnod);
  myCoord.resize(nnod);
  for (i1 = j1 = 0; i1 < myCoord.size(); i1++)
  {
    myMLGN[i1] = ++gNod;
    for (size_t d = 0; d < nsd; d++)
      myCoord[i1][d] = XYZ[j1+d];
    j1 += surf->dimension();
  }

  // Number of elements in patch
  nel = nelx*nely;
  // Number of nodes per element
  const int nen = p1*p2;

  // Connectivity array: local --> global node relation
  myMLGE.resize(nel);
  myMNPC.resize(nel);

  int i, j, a, b, iel = 0;
  for (j = 0; j < nely; j++)
    for (i = 0; i < nelx; i++, iel++)
    {
      myMLGE[iel] = ++gEl;
      myMNPC[iel].resize(nen);
      // First node in current element
      int corner = (p2-1)*nx*j + (p1-1)*i;

      for (b = 0; b < p2; b++)
      {
        int facenod = b*p1;
        myMNPC[iel][facenod] = corner + b*nx;
        for (a = 1; a < p1; a++)
          myMNPC[iel][facenod+a] = myMNPC[iel][facenod] + a;
      }
    }

  return true;
}


Vec3 ASMs2DLag::getCoord (size_t inod) const
{
  if (inod < 1 || inod > coord.size())
    return Vec3();

  return coord[inod-1];
}


void ASMs2DLag::setCoord (size_t inod, const Vec3& Xnod)
{
  if (inod < 1)
    return;

  if (inod > nnod)
    myCoord.resize(nnod = inod);

  myCoord[inod-1] = Xnod;
}


Vec3 ASMs2DLag::getElementCenter (int iel) const
{
  if (iel < 1 || static_cast<size_t>(iel) > MNPC.size())
  {
    std::cerr <<" *** ASMs2DLag::getElementCenter: Element index "<< iel
              <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return Vec3();
  }

  return this->getGeometricCenter(MNPC[--iel]);
}


bool ASMs2DLag::getElementCoordinates (Matrix& X, int iel, bool) const
{
  if (iel < 1 || static_cast<size_t>(iel) > MNPC.size())
  {
    std::cerr <<" *** ASMs2DLag::getElementCoordinates: Element index "<< iel
              <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return false;
  }

  // Number of nodes per element
  size_t nen = std::min(static_cast<size_t>(p1*p2),MNPC[--iel].size());

  X.resize(nsd,nen);
  for (size_t i = 0; i < nen; i++)
    X.fillColumn(i+1,coord[MNPC[iel][i]].ptr());

  return true;
}


void ASMs2DLag::getNodalCoordinates (Matrix& X, bool) const
{
  X.resize(nsd,coord.size());

  for (size_t inod = 0; inod < coord.size(); inod++)
    X.fillColumn(inod+1,coord[inod].ptr());
}


bool ASMs2DLag::updateCoords (const Vector& displ)
{
  if (shareFE) return true;

  if (displ.size() != nsd*coord.size())
  {
    std::cerr <<" *** ASMs2DLag::updateCoords: Invalid dimension "
              << displ.size() <<" on displ, should be "
              << nsd*coord.size() << std::endl;
    return false;
  }

  const double* u = displ.ptr();
  for (size_t inod = 0; inod < coord.size(); inod++, u += nsd)
    myCoord[inod] += RealArray(u,u+nsd);

  return true;
}


bool ASMs2DLag::getOrder (int& pu, int& pv, int& pw) const
{
  pu = p1;
  pv = p2;
  pw = 0;

  return true;
}


bool ASMs2DLag::getSize (int& n1, int& n2, int) const
{
  n1 = nx;
  n2 = ny;

  return true;
}


size_t ASMs2DLag::getNoBoundaryElms (char lIndex, char ldim) const
{
  if (ldim < 1 && lIndex > 0)
    return 1;

  switch (lIndex)
  {
    case 1:
    case 2:
      return ny > 1 && p2 > 1 ? (ny-1)/(p2-1) : 0;
    case 3:
    case 4:
      return nx > 1 && p1 > 1 ? (nx-1)/(p1-1) : 0;
  }

  return 0;
}


bool ASMs2DLag::integrate (Integrand& integrand,
                           GlobalIntegral& glInt,
                           const TimeDomain& time)
{
  if (this->empty()) return true; // silently ignore empty patches

  if (myCache.empty())
    myCache.emplace_back(std::make_unique<BasisFunctionCache>(*this));

  ASMs2D::BasisFunctionCache& cache = *myCache.front();
  cache.setIntegrand(&integrand);
  if (!cache.init(1))
    return false;

  // Get Gaussian quadrature points and weights
  const std::array<const double*,2>& xg = cache.coord();
  const std::array<const double*,2>& wg = cache.weight();

  // Get the reduced integration quadrature points, if needed
  const double* xr = cache.coord(true)[0];
  const double* wr = cache.weight(true)[0];

  // Number of elements in first parameter direction
  const int nelx = (nx-1)/(p1-1);


  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t g = 0; g < threadGroups.size() && ok; g++)
  {
#pragma omp parallel for schedule(static)
    for (const IntVec& group : threadGroups[g])
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
        if (!this->isElementActive(fe.iel)) continue; // zero-area element

        // Set up nodal point coordinates for current element
        this->getElementCoordinates(fe.Xn,1+iel);

        if (integrand.getIntegrandType() & Integrand::ELEMENT_CENTER)
        {
          // Compute the element "center" (average of element node coordinates)
          X = 0.0;
          for (size_t i = 1; i <= nsd; i++)
            for (size_t j = 1; j <= fe.Xn.cols(); j++)
              X[i-1] += fe.Xn(i,j);

          X *= 1.0 / static_cast<double>(fe.Xn.cols());
        }

        // Initialize element quantities
        LocalIntegral* A = integrand.getLocalIntegral(fe.Xn.cols(),fe.iel);
        const int nRed = fe.Xn.cols() < 4 ? 0 : cache.nGauss(true).front();
        if (!integrand.initElement(MNPC[iel],fe,X,nRed*nRed,*A))
        {
          A->destruct();
          ok = false;
          break;
        }

        int i1 = nelx > 0 ? iel % nelx : 0;
        int i2 = nelx > 0 ? iel / nelx : 0;

        if (xr)
        {
          // --- Selective reduced integration loop ----------------------------

          size_t ip = 0;
          for (int j = 0; j < nRed; j++)
            for (int i = 0; i < nRed; i++, ++ip)
            {
              // Local element coordinates of current integration point
              fe.xi  = xr[i];
              fe.eta = xr[j];

              // Parameter value of current integration point
              fe.u = cache.getParam(0,i1,i,true);
              fe.v = cache.getParam(1,i2,j,true);

              const BasisFunctionVals& bfs = cache.getVals(iel,ip,true);
              fe.N = bfs.N;

              // Compute Jacobian inverse and derivatives
              fe.detJxW = utl::Jacobian(Jac,fe.dNdX,fe.Xn,bfs.dNdu);
              if (fe.detJxW == 0.0) continue; // skip singular points

              // Store tangent vectors in fe.G for shells
              if (nsd > 2) fe.G = Jac;

              // Cartesian coordinates of current integration point
              X.assign(fe.Xn * fe.N);

              // Compute the reduced integration terms of the integrand
              fe.detJxW *= wr[i]*wr[j];
              if (!integrand.reducedInt(*A,fe,X))
                ok = false;
            }
        }


        // --- Integration loop over all Gauss points in each direction --------

        const int ng1 = fe.Xn.cols() < 4 ? 0 : cache.nGauss().front();
        const int ng2 = fe.Xn.cols() < 4 ? 0 : cache.nGauss().back();

        size_t ip = 0;
        int jp = iel*ng1*ng2;
        fe.iGP = firstIp + jp; // Global integration point counter

        for (int j = 0; j < ng2; j++)
          for (int i = 0; i < ng1; i++, fe.iGP++, ++ip)
          {
            // Local element coordinates of current integration point
            fe.xi  = xg[0][i];
            fe.eta = xg[1][j];

            // Parameter value of current integration point
            fe.u = cache.getParam(0,i1,i);
            fe.v = cache.getParam(1,i2,j);

            const BasisFunctionVals& bfs = cache.getVals(iel,ip);
            fe.N = bfs.N;

            // Compute Jacobian inverse of coordinate mapping and derivatives
            fe.detJxW = utl::Jacobian(Jac,fe.dNdX,fe.Xn,bfs.dNdu);
            if (fe.detJxW == 0.0) continue; // skip singular points

            // Store tangent vectors in fe.G for shells
            if (nsd > 2) fe.G = Jac;

            // Cartesian coordinates of current integration point
            X.assign(fe.Xn * fe.N);

            // Evaluate the integrand and accumulate element contributions
            fe.detJxW *= wg[0][i]*wg[1][j];
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

  cache.finalizeAssembly();
  return ok;
}


bool ASMs2DLag::integrate (Integrand& integrand, int lIndex,
                           GlobalIntegral& glInt,
                           const TimeDomain& time)
{
  if (this->empty()) return true; // silently ignore empty patches

  // Get Gaussian quadrature points and weights
  int nG1 = this->getNoGaussPt(lIndex%10 < 3 ? p2 : p1, true);
  int nGP = integrand.getBouIntegrationPoints(nG1);
  const double* xg = GaussQuadrature::getCoord(nGP);
  const double* wg = GaussQuadrature::getWeight(nGP);
  if (!xg || !wg) return false;

  // Find the parametric direction of the edge normal {-2,-1, 1, 2}
  const int edgeDir = (lIndex%10+1)/((lIndex%2) ? -2 : 2);

  const int t1 = abs(edgeDir); // tangent direction normal to the patch edge
  const int t2 = 3-t1;         // tangent direction along the patch edge

  // Number of elements in each direction
  const int nelx = (nx-1)/(p1-1);
  const int nely = (ny-1)/(p2-1);

  FiniteElement fe;
  RealArray upar, vpar;
  if (surf)
  {
    // Get parametric coordinates of the elements
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
  }

  // Extract the Neumann order flag (1 or higher) for the integrand
  integrand.setNeumannOrder(1 + lIndex/10);

  // Integrate the extraordinary elements?
  size_t doXelms = 0;
  if (integrand.getIntegrandType() & Integrand::XO_ELEMENTS)
    if ((doXelms = nelx*nely)*2 > MNPC.size())
    {
      std::cerr <<" *** ASMs2DLag::integrate: Too few XO-elements "
                << MNPC.size() - doXelms << std::endl;
      return false;
    }

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex%10);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;

  Matrix dNdu, Jac;
  Vec4   X(nullptr,time.t);
  Vec3   normal;
  double xi[2];


  // === Assembly loop over all elements on the patch edge =====================

  int iel = 1;
  for (int i2 = 0; i2 < nely; i2++)
    for (int i1 = 0; i1 < nelx; i1++, iel++)
    {
      fe.iel = abs(MLGE[doXelms+iel-1]);
      if (!this->isElementActive(fe.iel)) continue; // zero-area element

      // Skip elements that are not on current boundary edge
      bool skipMe = false;
      switch (edgeDir)
      {
        case -1: if (i1 > 0)      skipMe = true; break;
        case  1: if (i1 < nelx-1) skipMe = true; break;
        case -2: if (i2 > 0)      skipMe = true; break;
        case  2: if (i2 < nely-1) skipMe = true; break;
      }
      if (skipMe) continue;

      // Set up nodal point coordinates for current element
      this->getElementCoordinates(fe.Xn,iel);

      // Initialize element quantities
      LocalIntegral* A = integrand.getLocalIntegral(fe.Xn.cols(),fe.iel,true);
      bool ok = integrand.initElementBou(MNPC[doXelms+iel-1],*A);


      // --- Integration loop over all Gauss points along the edge -------------

      int jp = (t1 == 1 ? i2 : i1)*nGP;
      fe.iGP = firstp + jp; // Global integration point counter

      for (int i = 0; i < nGP && ok; i++, fe.iGP++)
      {
        // Local element coordinates of current integration point
        xi[t1-1] = edgeDir < 0 ? -1.0 : 1.0;
        xi[t2-1] = xg[i];
        fe.xi  = xi[0];
        fe.eta = xi[1];

        // Parameter values of current integration point
        if (upar.size() > 1)
          fe.u = 0.5*(upar[i1]*(1.0-xg[i]) + upar[i1+1]*(1.0+xg[i]));
        if (vpar.size() > 1)
          fe.v = 0.5*(vpar[i2]*(1.0-xg[i]) + vpar[i2+1]*(1.0+xg[i]));

        // Compute the basis functions and their derivatives, using
        // tensor product of one-dimensional Lagrange polynomials
        if (!Lagrange::computeBasis(fe.N,dNdu,p1,xi[0],p2,xi[1]))
          ok = false;

        // Compute basis function derivatives and the edge normal
        fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,fe.Xn,dNdu,t1,t2);
        if (fe.detJxW == 0.0) continue; // skip singular points

        if (edgeDir < 0) normal *= -1.0;

        // Store tangent vectors in fe.G for shells
        if (nsd > 2) fe.G = std::move(Jac);

        // Cartesian coordinates of current integration point
        X.assign(fe.Xn * fe.N);

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


int ASMs2DLag::evalPoint (const double* xi, double* param, Vec3& X) const
{
  if (surf)
  {
    param[0] = (1.0-xi[0])*surf->startparam_u() + xi[0]*surf->endparam_u();
    param[1] = (1.0-xi[1])*surf->startparam_v() + xi[1]*surf->endparam_v();
  }
  else
    memcpy(param,xi,2*sizeof(double));

  // Evaluate the parametric values of the nodes
  RealArray u, v;
  if (!this->getGridParameters(u,0,p1-1)) return -2;
  if (!this->getGridParameters(v,1,p2-1)) return -2;

  // Search for the closest node
  size_t i = utl::find_closest(u,param[0]);
  size_t j = utl::find_closest(v,param[1]);
  size_t n = u.size()*j + i;
  X = coord[n];

  return 1+n;
}


bool ASMs2DLag::tesselate (ElementBlock& grid, const int* npe) const
{
  if (p1 != npe[0] || p2 != npe[1])
  {
    int* newnpe = const_cast<int*>(npe);
    std::cout <<"\nLagrange elements: The number of visualization points are "
              << p1 <<" "<< p2 <<" by default\n"<< std::endl;
    newnpe[0] = p1;
    newnpe[1] = p2;
  }

  if (!this->ASMs2D::tesselate(grid,npe))
    return false;

  // Adjust the element Id since each Lagrange element covers several knot-spans
  int i, ie, nse1 = p1-1;
  int j, je, nse2 = p2-1;
  int nelx = (nx-1)/nse1;
  for (j = je = 1; j < (int)ny; j++)
  {
    for (i = ie = 1; i < (int)nx; i++)
    {
      grid.setElmId((j-1)*(nx-1)+i,(je-1)*nelx+ie);
      if (i%nse1 == 0) ie++;
    }
    if (j%nse2 == 0) je++;
  }

  return true;
}


void ASMs2DLag::constrainEdge (int dir, bool open, int dof, int code, char basis)
{
  this->ASMs2D::constrainEdge(dir, open, dof, code > 0 ? -code : code, basis);
}


bool ASMs2DLag::evalSolution (Matrix& sField, const Vector& locSol,
                              const int*, int n_f, bool) const
{
  return this->evalSolution(sField,locSol,nullptr,false,0,n_f);
}


bool ASMs2DLag::evalSolution (Matrix& sField, const Vector& locSol,
                              const RealArray* gpar, bool regular,
                              int, int) const
{
  if (!gpar && !regular) // Direct nodal evaluation
    return this->nodalField(sField,locSol,this->getNoNodes(-1));

  size_t nCmp = locSol.size() / this->getNoProjectionNodes();
  size_t ni   = gpar ? gpar[0].size() : nel;
  size_t nj   = gpar ? gpar[1].size() : 1;
  size_t nen  = p1*p2;

  sField.resize(nCmp,ni*nj);
  Matrix elmSol(nCmp,nen);
  RealArray N(nen), val;

  double xi = 0.0, eta = 0.0;
  if (!gpar && !Lagrange::computeBasis(N,p1,xi,p2,eta))
    return false;

  size_t ip = 1;
  int iel = 0;
  for (size_t j = 0; j < nj; j++)
    for (size_t i = 0; i < ni; i++, ip++)
    {
      if (gpar)
      {
        iel = this->findElement(gpar[0][i], gpar[1][j], &xi, &eta);
        if (iel < 1 || iel > static_cast<int>(nel))
          return false;
        if (!Lagrange::computeBasis(N,p1,xi,p2,eta))
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


bool ASMs2DLag::evalSolution (Matrix& sField, const IntegrandBase& integrand,
                              const int*, char) const
{
  sField.resize(0,0);

  double incx = 2.0/double(p1-1);
  double incy = 2.0/double(p2-1);

  size_t nPoints = coord.size();
  IntVec check(nPoints,0);

  FiniteElement fe(p1*p2);
  Vector        solPt;
  Vectors       globSolPt(nPoints);
  Matrix        dNdu, Jac;

  // Evaluate the secondary solution field at each point
  for (size_t iel = 0; iel < nel; iel++)
  {
    fe.iel = MLGE[iel];
    if (fe.iel < 1) continue; // zero-area element

    const IntVec& mnpc = MNPC[iel];
    this->getElementCoordinates(fe.Xn,1+iel);

    int i, j, loc = 0;
    for (j = 0; j < p2; j++)
      for (i = 0; i < p1; i++, loc++)
      {
        fe.xi  = -1.0 + i*incx;
        fe.eta = -1.0 + j*incy;
        if (!Lagrange::computeBasis(fe.N,dNdu,p1,fe.xi,p2,fe.eta))
          return false;

        // Compute the Jacobian inverse
        fe.detJxW = utl::Jacobian(Jac,fe.dNdX,fe.Xn,dNdu);
        if (fe.detJxW == 0.0) continue; // skip singular points

        // Now evaluate the solution field
        utl::Point X4(fe.Xn*fe.N,{fe.u,fe.v});
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


bool ASMs2DLag::evalSolution (Matrix& sField, const IntegrandBase& integrand,
                              const RealArray* gpar, bool regular) const
{
  sField.resize(0,0);

  bool elCenters = gpar ? regular && gpar->empty() : true;
  size_t nPoints = elCenters ? nel : gpar->size();

  FiniteElement fe(p1*p2);
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
      iel = this->findElement(gpar[0][i], gpar[1][i], &fe.xi, &fe.eta);
      if (iel < 1 || iel > static_cast<int>(nel))
        return false;
    }

    fe.iel = MLGE[iel-1];
    if (fe.iel < 1) continue; // zero-area element

    this->getElementCoordinates(fe.Xn,iel);

    if (!Lagrange::computeBasis(fe.N,dNdu,p1,fe.xi,p2,fe.eta))
      return false;

    // Compute the Jacobian inverse
    fe.detJxW = utl::Jacobian(Jac,fe.dNdX,fe.Xn,dNdu);
    if (fe.detJxW == 0.0) continue; // skip singular points

    // Store tangent vectors in fe.G for shells
    if (nsd > 2) fe.G = std::move(Jac);

    // Now evaluate the solution field
    utl::Point X4(fe.Xn*fe.N,{fe.u,fe.v});
    if (!integrand.evalSol(solPt,fe,X4,MNPC[iel-1]))
      return false;
    else if (sField.empty())
      sField.resize(solPt.size(),nPoints,true);

    sField.fillColumn(1+i,solPt);
  }

  return true;
}


void ASMs2DLag::generateThreadGroups (const Integrand&, bool, bool)
{
  threadGroups.calcGroups((nx-1)/(p1-1),(ny-1)/(p2-1),1);
}


void ASMs2DLag::findBoundaryElms (IntVec& elms, int lIndex, int) const
{
  const int N1m = (nx-1)/(p1-1);
  const int N2m = (ny-1)/(p2-1);

  elms.clear();
  switch (lIndex) {
  case 1:
  case 2:
    elms.reserve(N2m);
    for (int i = 0; i < N2m; ++i)
      elms.push_back(i*N1m + (lIndex-1)*(N1m-1));
    break;
  case 3:
  case 4:
    elms.reserve(N1m);
    for (int i = 0; i < N1m; ++i)
      elms.push_back(i + (lIndex-3)*N1m*(N2m-1));
  }
}


bool ASMs2DLag::write (std::ostream& os, int) const
{
  return this->writeLagBasis(os, "quad");
}


int ASMs2DLag::findElement (double u, double v, double* xi, double* eta) const
{
  if (!surf)
  {
    std::cerr <<" *** ASMs2DLag::findElement: No spline geometry"<< std::endl;
    return -1;
  }

  int ku = surf->basis(0).knotInterval(u);
  int kv = surf->basis(1).knotInterval(v);
  int elmx = ku - (p1 - 1);
  int elmy = kv - (p2 - 1);

  if (xi) {
    const double knot_1 = *(surf->basis(0).begin() + ku);
    const double knot_2 = *(surf->basis(0).begin() + ku+1);
    *xi = -1.0 + 2.0 * (u - knot_1) / (knot_2 - knot_1);
  }
  if (eta) {
    const double knot_1 = *(surf->basis(1).begin() + kv);
    const double knot_2 = *(surf->basis(1).begin() + kv+1);
    *eta = -1.0 + 2.0 * (v - knot_1) / (knot_2 - knot_1);
  }

  int nel1 = surf->numCoefs_u() - p1 + 1;

  return 1 + elmx + elmy*nel1;
}


bool ASMs2DLag::evaluate (const ASMbase* basis, const Vector& locVec,
                          RealArray& vec, int basisNum) const
{
  // Compute parameter values of the result sampling points
  std::array<RealArray,2> gpar;
  for (int dir = 0; dir < 2; dir++) {
    int nel1 = dir == 0 ? (nx-1)/(p1-1) : (ny-1)/(p2-1);
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


size_t ASMs2DLag::getNoProjectionNodes () const
{
  return this->getNoNodes(1);
}


bool ASMs2DLag::assembleL2matrices (SparseMatrix& A, StdVector& B,
                                    const L2Integrand& integrand,
                                    bool continuous) const
{
  ASMs2D::BasisFunctionCache& cache = *myCache.front();
  if (!cache.init(1))
    return false;

  const std::array<const double*,2>& wg = cache.weight();
  const std::array<int,2> nGP = cache.nGauss();

  Matrix dNdX, Xnod, J;

  const size_t nnod = this->getNoProjectionNodes();

  // === Assembly loop over all elements in the patch ==========================

  int iel = 0;
  for (size_t i2 = 0; i2 < cache.noElms()[1]; ++i2)
    for (size_t i1 = 0; i1 < cache.noElms()[0]; ++i1, ++iel)
    {
      if (!this->getElementCoordinates(Xnod,1+iel))
        return false;

      std::array<RealArray,2> GP;
      GP[0].reserve(nGP[0]*nGP[1]);
      GP[1].reserve(nGP[0]*nGP[1]);
      for (int j = 0; j < nGP[1]; ++j)
        for (int i = 0; i < nGP[0]; ++i) {
          GP[0].push_back(cache.getParam(0, i1, i));
          GP[1].push_back(cache.getParam(1, i2, j));
        }

      Matrix sField;
      integrand.evaluate(sField, GP.data());

      // --- Integration loop over all Gauss points in each direction ----------

      Matrix eA(p1*p2, p1*p2);
      Vectors eB(sField.rows(), Vector(p1*p2));
      size_t ip = 0;
      for (int j = 0; j < nGP[1]; ++j)
        for (int i = 0; i < nGP[0]; ++i, ++ip) {
          const BasisFunctionVals& bfs = cache.getVals(iel,ip);

          double dJw = wg[0][i]*wg[1][j]*utl::Jacobian(J,dNdX,Xnod,bfs.dNdu);

          // Integrate the mass matrix
          eA.outer_product(bfs.N, bfs.N, true, dJw);

          // Integrate the rhs vector B
          for (size_t r = 1; r <= sField.rows(); r++)
            eB[r-1].add(bfs.N,sField(r,ip+1)*dJw);
        }

      const IntVec& mnpc = MNPC[iel];

      for (int i = 0; i < p1*p2; ++i) {
        for (int j = 0; j < p1*p2; ++j)
          A(mnpc[i]+1, mnpc[j]+1) += eA(i+1, j+1);

        int jp = mnpc[i]+1;
        for (size_t r = 0; r < sField.rows(); r++, jp += nnod)
          B(jp) += eB[r](1+i);
      }
    }

  return true;
}


ASMs2DLag::BasisFunctionCache::BasisFunctionCache (const ASMs2DLag& pch)
  : ASMs2D::BasisFunctionCache(pch)
{
  nel[0] = (pch.nx-1) / (pch.p1-1);
  nel[1] = (pch.ny-1) / (pch.p2-1);
}


ASMs2DLag::BasisFunctionCache::BasisFunctionCache (const ASMs2D::BasisFunctionCache& cache,
                                                   int b) :
  ASMs2D::BasisFunctionCache(cache,b)
{
}


bool ASMs2DLag::BasisFunctionCache::internalInit ()
{
  if (!mainQ->xg[0])
    this->setupQuadrature();

  nTotal = mainQ->ng[0]*mainQ->ng[1];
  if (reducedQ->xg[0])
    nTotalRed = reducedQ->ng[0]*reducedQ->ng[1];

  return true;
}


void ASMs2DLag::BasisFunctionCache::setupParameters ()
{
  // Compute parameter values of the Gauss points over the whole patch
  for (int d = 0; d < 2; d++) {
    RealArray par;
    patch.getGridParameters(par,d,1);
    mainQ->gpar[d].resize(par.size(), 1);
    mainQ->gpar[d].fillColumn(1,par.data());
    if (reducedQ->xg[0])
      reducedQ->gpar[d] = mainQ->gpar[d];
  }
}


double ASMs2DLag::BasisFunctionCache::getParam (int dir, size_t el,
                                                size_t gp, bool reduced) const
{
  const Quadrature& q = reduced ? *reducedQ : *mainQ;
  return 0.5*(q.gpar[dir](el+1,1)*(1.0-q.xg[dir][gp]) +
              q.gpar[dir](el+2,1)*(1.0+q.xg[dir][gp]));
}


BasisFunctionVals ASMs2DLag::BasisFunctionCache::calculatePt (size_t, size_t gp,
                                                              bool red) const
{
  std::array<size_t,2> gpIdx = this->gpIndex(gp,red);
  const Quadrature& q = red ? *reducedQ : *mainQ;

  int p1, p2, p3;
  patch.getOrder(p1,p2,p3);

  BasisFunctionVals result;
  Lagrange::computeBasis(result.N,result.dNdu,
                         p1,q.xg[0][gpIdx[0]],
                         p2,q.xg[1][gpIdx[1]]);
  return result;
}


void ASMs2DLag::BasisFunctionCache::calculateAll ()
{
  int i, j, p1, p2;
  patch.getOrder(p1,p2,i);

  // Evaluate basis function values and derivatives at all integration points.
  // They will be the same for all elements in the patch,
  // so no need to repeat for all.
  std::vector<BasisFunctionVals>::iterator it = values.begin();
  for (j = 0; j < mainQ->ng[1]; j++)
    for (i = 0; i < mainQ->ng[0]; i++, ++it)
      Lagrange::computeBasis(it->N, it->dNdu,
                             p1, mainQ->xg[0][i], p2, mainQ->xg[1][j]);

  if (reducedQ->xg[0])
    for (j = 0, it = valuesRed.begin(); j < reducedQ->ng[1]; j++)
      for (i = 0; i < reducedQ->ng[0]; i++, ++it)
        Lagrange::computeBasis(it->N, it->dNdu,
                               p1, reducedQ->xg[0][i],  p2, reducedQ->xg[1][j]);
}
