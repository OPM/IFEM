// $Id$
//==============================================================================
//!
//! \file ASMu2Drecovery.C
//!
//! \date February 2012
//!
//! \author Kjetil A. Johannessen and Knut Morten Okstad / SINTEF
//!
//! \brief Recovery techniques for unstructured LR B-splines.
//!
//==============================================================================

#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Element.h"

#include "ASMu2D.h"
#include "FiniteElement.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "SparseMatrix.h"
#include "DenseMatrix.h"
#include "Function.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Profiler.h"
#include <array>
#include <fstream>


bool ASMu2D::getGrevilleParameters (RealArray& prm, int dir) const
{
  if (!lrspline || dir < 0 || dir > 1) return false;

  prm.clear();
  prm.reserve(lrspline->nBasisFunctions());
  for(LR::Basisfunction *b : lrspline->getAllBasisfunctions())
    prm.push_back(b->getGrevilleParameter()[dir]);

  return true;
}


/*!
  \brief Expands a tensor parametrization point to an unstructured one.
  \details Takes as input a tensor mesh, for instance
     in[0] = {0,1,2}
     in[1] = {2,3,5}
   and expands this to an unstructred representation, i.e.,
     out[0] = {0,1,2,0,1,2,0,1,2}
     out[1] = {2,2,2,3,3,3,5,5,5}
*/

static void expandTensorGrid (const RealArray* in, RealArray* out)
{
  out[0].resize(in[0].size()*in[1].size());
  out[1].resize(in[0].size()*in[1].size());

  size_t i, j, ip = 0;
  for (j = 0; j < in[1].size(); j++)
    for (i = 0; i < in[0].size(); i++, ip++) {
      out[0][ip] = in[0][i];
      out[1][ip] = in[1][j];
    }
}


LR::LRSplineSurface* ASMu2D::projectSolution (const IntegrandBase& integr) const
{
  PROFILE2("ASMu2D::projectSolution");

  // Compute parameter values of the result sampling points (Greville points)
  std::array<RealArray,2> gpar;
  for (int dir = 0; dir < 2; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir))
      return nullptr;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integr,gpar.data()))
    return nullptr;

  // Project the results onto the spline basis to find control point
  // values based on the result values evaluated at the Greville points.
  // Note that we here implicitly assume that the number of Greville points
  // equals the number of control points such that we don't have to resize
  // the result array. Think that is always the case, but beware if trying
  // other projection schemes later.

  return this->regularInterpolation(gpar[0],gpar[1],sValues);
}


LR::LRSpline* ASMu2D::evalSolution (const IntegrandBase& integrand) const
{
  return this->projectSolution(integrand);
}


bool ASMu2D::assembleL2matrices (SparseMatrix& A, StdVector& B,
                                 const IntegrandBase& integrand,
                                 bool continuous) const
{
  const size_t nnod = this->getNoNodes();

  const int p1 = lrspline->order(0);
  const int p2 = lrspline->order(1);

  // Get Gaussian quadrature points
  const int ng1 = continuous ? nGauss : p1 - 1;
  const int ng2 = continuous ? nGauss : p2 - 1;
  const double* xg = GaussQuadrature::getCoord(ng1);
  const double* yg = GaussQuadrature::getCoord(ng2);
  const double* wg = continuous ? GaussQuadrature::getWeight(nGauss) : nullptr;
  if (!xg || !yg) return false;
  if (continuous && !wg) return false;

  double dA = 0.0;
  Vector phi;
  Matrix dNdu, Xnod, Jac;
  Go::BasisDerivsSf spl1;
  Go::BasisPtsSf    spl0;


  // === Assembly loop over all elements in the patch ==========================

  std::vector<LR::Element*>::iterator el = lrspline->elementBegin();
  for (int iel = 1; el != lrspline->elementEnd(); el++, iel++)
  {
    if (continuous)
    {
      // Set up control point (nodal) coordinates for current element
      if (!this->getElementCoordinates(Xnod,iel))
        return false;
      else if ((dA = 0.25*this->getParametricArea(iel)) < 0.0)
        return false; // topology error (probably logic error)
    }

    // Compute parameter values of the Gauss points over this element
    std::array<RealArray,2> gpar, unstrGpar;
    this->getGaussPointParameters(gpar[0],0,ng1,iel,xg);
    this->getGaussPointParameters(gpar[1],1,ng2,iel,yg);

    // convert to unstructred mesh representation
    expandTensorGrid(gpar.data(),unstrGpar.data());

    // Evaluate the secondary solution at all integration points
    Matrix sField;
    if (!this->evalSolution(sField,integrand,unstrGpar.data()))
      return false;

    // set up basis function size (for extractBasis subroutine)
    phi.resize((**el).nBasisFunctions());

    // --- Integration loop over all Gauss points in each direction ----------
    int ip = 0;
    for (int j = 0; j < ng2; j++)
      for (int i = 0; i < ng1; i++, ip++)
      {
        if (continuous)
        {
          lrspline->computeBasis(gpar[0][i],gpar[1][j],spl1,iel-1);
          SplineUtils::extractBasis(spl1,phi,dNdu);
        }
        else
        {
          lrspline->computeBasis(gpar[0][i],gpar[1][j],spl0,iel-1);
          phi = spl0.basisValues;
        }

        // Compute the Jacobian inverse and derivatives
        double dJw = 1.0;
        if (continuous)
        {
          dJw = dA*wg[i]*wg[j]*utl::Jacobian(Jac,dNdu,Xnod,dNdu,false);
          if (dJw == 0.0) continue; // skip singular points
        }

        // Integrate the linear system A*x=B
        for (size_t ii = 0; ii < phi.size(); ii++)
        {
          int inod = MNPC[iel-1][ii]+1;
          for (size_t jj = 0; jj < phi.size(); jj++)
          {
            int jnod = MNPC[iel-1][jj]+1;
            A(inod,jnod) += phi[ii]*phi[jj]*dJw;
          }
          for (size_t r = 1; r <= sField.rows(); r++)
            B(inod+(r-1)*nnod) += phi[ii]*sField(r,ip+1)*dJw;
        }
      }
  }

  return true;
}


/*!
  \brief Evaluates monomials in 2D up to the specified order.
*/

static void evalMonomials (int p1, int p2, double x, double y, Vector& P)
{
  P.resize(p1*p2);
  P.fill(1.0);

  int i, j, k, ip = 1;
  for (j = 0; j < p2; j++)
    for (i = 0; i < p1; i++, ip++)
    {
      for (k = 0; k < i; k++) P(ip) *= x;
      for (k = 0; k < j; k++) P(ip) *= y;
    }
}


LR::LRSplineSurface* ASMu2D::scRecovery (const IntegrandBase& integrand) const
{
  PROFILE2("ASMu2D::scRecovery");

  const int m = integrand.derivativeOrder();
  const int p1 = lrspline->order(0);
  const int p2 = lrspline->order(1);

  // Get Gaussian quadrature point coordinates
  const int ng1 = p1 - m;
  const int ng2 = p2 - m;
  const double* xg = GaussQuadrature::getCoord(ng1);
  const double* yg = GaussQuadrature::getCoord(ng2);
  if (!xg || !yg) return nullptr;

  // Compute parameter values of the Greville points
  std::array<RealArray,2> gpar;
  if (!this->getGrevilleParameters(gpar[0],0)) return nullptr;
  if (!this->getGrevilleParameters(gpar[1],1)) return nullptr;

  const int n1 = p1 - m + 1; // Patch size in first parameter direction
  const int n2 = p2 - m + 1; // Patch size in second parameter direction

  const size_t nCmp = integrand.getNoFields(); // Number of result components
  const size_t nPol = n1*n2; // Number of terms in polynomial expansion

  Matrix sValues(nCmp,gpar[0].size());
  Vector P(nPol);
  Go::Point X, G;

  // Loop over all Greville points (one for each basis function)
  size_t k, l, ip = 0;
  std::vector<LR::Element*>::const_iterator elStart, elEnd, el;
  std::vector<LR::Element*> supportElements;
  for (LR::Basisfunction *b : lrspline->getAllBasisfunctions())
  {
#if SP_DEBUG > 2
    std::cout <<"Basis: "<< *b <<"\n  ng1 ="<< ng1 <<"\n  ng2 ="<< ng2
              <<"\n  nPol="<< nPol << std::endl;
#endif

    // Special case for basis functions with too many zero knot spans by using
    // the extended support
    // if (nel*ng1*ng2 < nPol)
    if (true)
    {
      // I am unsure as to the consequence of going back to previous if-statement
      // here so we keep if (true) for now. This was introduced mainly when considering
      // functions that live on the boundary and have support on few elements;
      // corner functions have support on one element. Using i.e. 2x2 points
      // for every element is not enough to fit 1,x,x^2,x^3,y,xy,...x^3y^3 when
      // we only have one element. The solution is getExtendedSupport, which is the
      // union of support from all functions that overlap *b.
      supportElements = b->getExtendedSupport();
      elStart = supportElements.begin();
      elEnd   = supportElements.end();
#if SP_DEBUG > 2
      std::cout <<"Extended basis:";
      for (el = elStart; el != elEnd; el++)
        std::cout <<"\n  " << **el;
      std::cout << std::endl;
#endif
    }
    else
    {
      elStart = b->supportedElementBegin();
      elEnd   = b->supportedElementEnd();
    }

    // Physical coordinates of current Greville point
    lrspline->point(G,gpar[0][ip],gpar[1][ip]);

    // Set up the local projection matrices
    DenseMatrix A(nPol,nPol);
    Matrix B(nPol,nCmp);

    // Loop over all non-zero knot-spans in the support of
    // the basis function associated with current Greville point
    for (el = elStart; el != elEnd; el++)
    {
      int iel = (**el).getId()+1;

      // evaluate all gauss points for this element
      std::array<RealArray,2> gaussPt, unstrGauss;
      this->getGaussPointParameters(gaussPt[0],0,ng1,iel,xg);
      this->getGaussPointParameters(gaussPt[1],1,ng2,iel,yg);

#if SP_DEBUG > 2
      std::cout << "Element " << **el << std::endl;
#endif

      // convert to unstructred mesh representation
      expandTensorGrid(gaussPt.data(),unstrGauss.data());

      // Evaluate the secondary solution at all Gauss points
      Matrix sField;
      if (!this->evalSolution(sField,integrand,unstrGauss.data()))
        return nullptr;

      // Loop over the Gauss points in current knot-span
      int i, j, ig = 1;
      for (j = 0; j < ng2; j++)
        for (i = 0; i < ng1; i++, ig++)
        {
          // Evaluate the polynomial expansion at current Gauss point
          lrspline->point(X,gaussPt[0][i],gaussPt[1][j]);
          evalMonomials(n1,n2,X[0]-G[0],X[1]-G[1],P);
#if SP_DEBUG > 2
          std::cout <<"Greville     point: "<< G
                    <<"\nGauss        point: "<< gaussPt[0][i] <<", "<< gaussPt[0][j]
                    <<"\nMapped gauss point: "<< X
                    <<"\nP-matrix:"<< P <<"--------------------\n"<< std::endl;
#endif

          for (k = 1; k <= nPol; k++)
          {
            // Accumulate the projection matrix, A += P^t * P
            for (l = 1; l <= nPol; l++)
              A(k,l) += P(k)*P(l);

            // Accumulate the right-hand-side matrix, B += P^t * sigma
            for (l = 1; l <= nCmp; l++)
              B(k,l) += P(k)*sField(l,ig);
          }
        }
    }

#if SP_DEBUG > 2
    std::cout <<"---- Matrix A -----"<< A
              <<"-------------------"<< std::endl;
    std::cout <<"---- Vector B -----"<< B
              <<"-------------------"<< std::endl;
#endif

    // Solve the local equation system
    if (!A.solve(B)) return nullptr;

    // Evaluate the projected field at current Greville point (first row of B)
    ip++;
    for (l = 1; l <= nCmp; l++)
      sValues(l,ip) = B(1,l);

#if SP_DEBUG > 1
    std::cout <<"Greville point "<< ip <<" (u,v) = "
              << gpar[0][ip-1] <<" "<< gpar[1][ip-1] <<" :";
    for (k = 1; k <= nCmp; k++)
      std::cout <<" "<< sValues(k,ip);
    std::cout << std::endl;
#endif
  }

  // Project the Greville point results onto the spline basis
  // to find the control point values

  return this->regularInterpolation(gpar[0],gpar[1],sValues);
}


LR::LRSplineSurface* ASMu2D::regularInterpolation (const RealArray& upar,
                                                   const RealArray& vpar,
                                                   const Matrix& points) const
{
  if (lrspline->rational())
  {
    std::cerr <<" *** ASMu2D::regularInterpolation:"
              <<" Rational LR B-splines not supported yet."<< std::endl;
    return nullptr;
  }

  // sanity check on input parameters
  const size_t nBasis = lrspline->nBasisFunctions();
  if (upar.size() != nBasis || vpar.size() != nBasis || points.cols() != nBasis)
  {
    std::cerr <<" *** ASMu2D::regularInterpolation:"
              <<" Mismatching input array sizes.\n"
              <<"     size(upar)="<< upar.size() <<" size(vpar)="<< vpar.size()
              <<" size(points)="<< points.cols() <<" nBasis="<< nBasis
              << std::endl;
    return nullptr;
  }

  SparseMatrix   A(SparseMatrix::SUPERLU);
  A.resize(nBasis, nBasis);
  Matrix B2(points,true); // transpose to get one vector per field
  StdVector B(B2);
  Go::BasisPtsSf splineValues;

  // Evaluate all basis functions at all points, stored in the A-matrix
  // (same row = same evaluation point)
  for (size_t i = 0; i < nBasis; i++)
  {
    int iel = lrspline->getElementContaining(upar[i], vpar[i]);
    lrspline->computeBasis(upar[i],vpar[i],splineValues, iel);
    size_t k = 0;
    for (auto& function : lrspline->getElement(iel)->support()) {
      int j = function->getId();
      A(i+1,j+1) = splineValues.basisValues[k++];
    }
  }

  // Solve for all solution components - one right-hand-side for each
  if (!A.solve(B))
    return nullptr;

  // Copy all basis functions and mesh
  LR::LRSplineSurface* ans = lrspline->copy();
  ans->rebuildDimension(points.rows());

  // Back to interleaved data
  std::vector<double> interleave;
  interleave.reserve(B.dim());
  for (size_t i = 0; i < nBasis; ++i)
    for (size_t j = 0; j < points.rows(); j++) {
        interleave.push_back(B(1+j*points.cols()+i));
  }

  ans->setControlPoints(interleave);

  return ans;
}


bool ASMu2D::edgeL2projection (const DirichletEdge& edge,
                               const RealFunc& values,
                               RealArray& result,
                               double time) const
{
  size_t n = edge.MLGN.size();
  SparseMatrix A(SparseMatrix::SUPERLU);
  StdVector    B(n);
  A.resize(n,n);

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // find the normal and tangent direction for the edge
  int edgeDir, t1, t2;
  switch (edge.edg)
  {       //          id          normal  tangent
    case LR::WEST:  edgeDir = -1;  t1=1;    t2=2; break;
    case LR::EAST:  edgeDir = +1;  t1=1;    t2=2; break;
    case LR::SOUTH: edgeDir = -2;  t1=2;    t2=1; break;
    case LR::NORTH: edgeDir = +2;  t1=2;    t2=1; break;
    default:        return false;
  }

  std::array<Vector,2> gpar;
  for (int d = 0; d < 2; d++)
    if (-1-d == edgeDir)
    {
      gpar[d].resize(nGauss);
      gpar[d].fill(d == 0 ? lrspline->startparam(0) : lrspline->startparam(1));
    }
    else if (1+d == edgeDir)
    {
      gpar[d].resize(nGauss);
      gpar[d].fill(d == 0 ? lrspline->endparam(0) : lrspline->endparam(1));
    }


  Matrix dNdu, Xnod, Jac;
  Vec4   X;
  Vec3   normal;
  int    iGp = 0;

  // === Assembly loop over all elements on the patch edge =====================
  for (size_t i=0; i<edge.MLGE.size(); i++) // for all edge elements
  {
    FiniteElement fe(edge.MNPC[i].size());
    fe.iel = edge.MLGE[i];

    // Get element edge length in the parameter space
    double dS = this->getParametricLength(fe.iel+1,t1);
    if (dS < 0.0) return false; // topology error (probably logic error)

    // Set up control point coordinates for current element
    if (!this->getElementCoordinates(Xnod,fe.iel+1)) return false;

    // Initialize element quantities
    fe.xi = fe.eta = edgeDir < 0 ? -1.0 : 1.0;

    // Get integration gauss points over this element
    this->getGaussPointParameters(gpar[t2-1],t2-1,nGauss,fe.iel+1,xg);

    // --- Integration loop over all Gauss points along the edge -------------
    for (int j = 0; j < nGauss; j++)
    {
      fe.iGP = iGp++; // Global integration point counter

      // Local element coordinates and parameter values
      // of current integration point
      fe.xi  = xg[j];
      fe.eta = xg[j];
      fe.u = gpar[0][j];
      fe.v = gpar[1][j];

      // Evaluate basis function (geometry) derivatives at current integration points
      Go::BasisDerivsSf spline;
      lrspline->computeBasis(fe.u, fe.v, spline, fe.iel);

      // Fetch basis function derivatives at current integration point
      SplineUtils::extractBasis(spline,fe.N,dNdu);

      // Compute basis function derivatives and the edge normal
      fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,Xnod,dNdu,t1,t2);
      if (fe.detJxW == 0.0) continue; // skip singular points

      if (edgeDir < 0) normal *= -1.0;

      // Cartesian coordinates of current integration point
      X = Xnod * fe.N;
      X.t = time;

      // Evaluate the integrand and accumulate element contributions
      fe.detJxW *= 0.5*dS*wg[j];

      // For mixed basis, we need to compute functions separate from geometry
      if (edge.lr != lrspline.get())
      {
        // different lrspline instances enumerate elements differently
        int basis_el =  edge.lr->getElementContaining(fe.u, fe.v);
        edge.lr->computeBasis(fe.u, fe.v, spline, basis_el);
        SplineUtils::extractBasis(spline,fe.N,dNdu);
      }

      // Assemble into matrix A
      for (size_t il=0; il<edge.MNPC[i].size(); il++)     // local i-index
        if (edge.MNPC[i][il] != -1)
        {
          size_t ig = edge.MNPC[i][il]+1;                 // global i-index
          for (size_t jl=0; jl<edge.MNPC[i].size(); jl++) // local j-index
            if (edge.MNPC[i][jl] != -1)
            {
              size_t jg = edge.MNPC[i][jl]+1;             // global j-index
              A(ig,jg) = A(ig,jg) + fe.N[il]*fe.N[jl]*fe.detJxW;
          }
          B(ig) = B(ig) + fe.N[il]*values(X)*fe.detJxW;
      } // end basis-function loop
    } // end gauss-point loop
  } // end element loop

#if SP_DEBUG > 2
  std::cout <<"---- Matrix A -----\n"<< A
            <<"-------------------"<< std::endl;
  std::cout <<"---- Vector B -----\n"<< B
            <<"-------------------"<< std::endl;

  // dump mesh enumerations to file
  std::ofstream out("mesh.eps");
  lrspline->writePostscriptFunctionSpace(out);

  std::cout <<"---- Edge-nodes (g2l-mapping) -----\n";
  int i=-1;
  for (auto d : edge.MLGN) {
    i++;
    std::cout << d.first << ": " << d.second << std::endl;
  }
  std::cout <<"-------------------"<< std::endl;
  std::cout <<"---- Element-nodes (-1 is interior element nodes) -----\n";
  for (auto d : edge.MNPC) {
    for (int c : d)
      std::cout << c << " ";
    std::cout << std::endl;
  }
  std::cout <<"-------------------"<< std::endl;
#endif

  // Solve the edge-global equation system
  if (!A.solve(B)) return false;

#if SP_DEBUG > 2
  std::cout <<"---- SOLUTION -----\n"<< B
            <<"-------------------"<< std::endl;
#endif

  // Store the control-point values of the projected field
  result.resize(n);
  for (size_t i = 0; i < n; i++)
    result[i] = B[i];

  return true;
}


bool ASMu2D::edgeL2projection (const DirichletEdge& edge,
                               const VecFunc& values,
                               Real2DMat& result,
                               double time) const
{
  size_t n = edge.MLGN.size();
  SparseMatrix A(SparseMatrix::SUPERLU);
  std::vector<StdVector> B(nsd, StdVector(n));
  A.resize(n,n);

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // find the normal and tangent direction for the edge
  int edgeDir, t1, t2;
  switch (edge.edg)
  {       //          id          normal  tangent
    case LR::WEST:  edgeDir = -1;  t1=1;    t2=2; break;
    case LR::EAST:  edgeDir = +1;  t1=1;    t2=2; break;
    case LR::SOUTH: edgeDir = -2;  t1=2;    t2=1; break;
    case LR::NORTH: edgeDir = +2;  t1=2;    t2=1; break;
    default:        return false;
  }

  std::array<Vector,2> gpar;
  for (int d = 0; d < 2; d++)
    if (-1-d == edgeDir)
    {
      gpar[d].resize(nGauss);
      gpar[d].fill(d == 0 ? lrspline->startparam(0) : lrspline->startparam(1));
    }
    else if (1+d == edgeDir)
    {
      gpar[d].resize(nGauss);
      gpar[d].fill(d == 0 ? lrspline->endparam(0) : lrspline->endparam(1));
    }


  Matrix dNdu, Xnod, Jac;
  Vec4   X;
  Vec3   normal;
  int    iGp = 0;

  // === Assembly loop over all elements on the patch edge =====================
  for (size_t i=0; i<edge.MLGE.size(); i++) // for all edge elements
  {
    FiniteElement fe(edge.MNPC[i].size());
    fe.iel = edge.MLGE[i];

    // Get element edge length in the parameter space
    double dS = this->getParametricLength(fe.iel+1,t1);
    if (dS < 0.0) return false; // topology error (probably logic error)

    // Set up control point coordinates for current element
    if (!this->getElementCoordinates(Xnod,fe.iel+1)) return false;

    // Initialize element quantities
    fe.xi = fe.eta = edgeDir < 0 ? -1.0 : 1.0;

    // Get integration gauss points over this element
    this->getGaussPointParameters(gpar[t2-1],t2-1,nGauss,fe.iel+1,xg);

    // --- Integration loop over all Gauss points along the edge -------------
    for (int j = 0; j < nGauss; j++)
    {
      fe.iGP = iGp++; // Global integration point counter

      // Local element coordinates and parameter values
      // of current integration point
      fe.xi  = xg[j];
      fe.eta = xg[j];
      fe.u = gpar[0][j];
      fe.v = gpar[1][j];

      // Evaluate basis function derivatives at current integration points
      Go::BasisDerivsSf spline;
      lrspline->computeBasis(fe.u, fe.v, spline, fe.iel);

      // Fetch basis function derivatives at current integration point
      SplineUtils::extractBasis(spline,fe.N,dNdu);

      // Compute basis function derivatives and the edge normal
      fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,Xnod,dNdu,t1,t2);
      if (fe.detJxW == 0.0) continue; // skip singular points

      if (edgeDir < 0) normal *= -1.0;

      // Cartesian coordinates of current integration point
      X = Xnod * fe.N;
      X.t = time;

      // Evaluate the integrand and accumulate element contributions
      fe.detJxW *= 0.5*dS*wg[j];

      // For mixed basis, we need to compute functions separate from geometry
      if (edge.lr != lrspline.get())
      {
        // different lrspline instances enumerate elements differently
        int basis_el =  edge.lr->getElementContaining(fe.u, fe.v);
        edge.lr->computeBasis(fe.u, fe.v, spline, basis_el);
        SplineUtils::extractBasis(spline,fe.N,dNdu);
      }

      // Assemble into matrix A
      for (size_t il=0; il<edge.MNPC[i].size(); il++)     // local i-index
        if (edge.MNPC[i][il] != -1)
        {
          size_t ig = edge.MNPC[i][il]+1;                 // global i-index
          for (size_t jl=0; jl<edge.MNPC[i].size(); jl++) // local j-index
            if (edge.MNPC[i][jl] != -1)
            {
              size_t jg = edge.MNPC[i][jl]+1;             // global j-index
              A(ig,jg) = A(ig,jg) + fe.N[il]*fe.N[jl]*fe.detJxW;
          }
          Vec3 val = values(X);
          for(size_t k=0; k<nsd; k++)
            B[k](ig) = B[k](ig) + fe.N[il]*val[k]*fe.detJxW;
      } // end basis-function loop
    } // end gauss-point loop
  } // end element loop

#if SP_DEBUG > 2
  std::cout <<"---- Matrix A -----\n"<< A
            <<"-------------------"<< std::endl;
  for(size_t k=0; k<nsd; k++)
  {
    std::cout <<"--- Vector B(" << (k+1) << ") ----\n" << B[k]
              <<"-------------------"<< std::endl;
  }

  // dump mesh enumerations to file
  std::ofstream out("mesh.eps");
  lrspline->writePostscriptFunctionSpace(out);

  std::cout <<"---- Edge-nodes (g2l-mapping) -----\n";
  int i=-1;
  for (auto d : edge.MLGN) {
    i++;
    std::cout << d.first << ": " << d.second << std::endl;
  }
  std::cout <<"-------------------"<< std::endl;
  std::cout <<"---- Element-nodes (-1 is interior element nodes) -----\n";
  for (auto d : edge.MNPC) {
    for (int c : d)
      std::cout << c << " ";
    std::cout << std::endl;
  }
  std::cout <<"-------------------"<< std::endl;
#endif

  // Solve the edge-global equation system
  if (!A.solve(B[0], true)) return false;
  // Solve the system for the rest of the right-hand-side components (re-use LU factorization)
  for (size_t k=1; k<nsd; k++)
    if (!A.solve(B[k], false)) return false;

#if SP_DEBUG > 2
  for(size_t k=0; k<nsd; k++)
  {
    std::cout <<"--- Solution(" << (k+1) << ") ----\n" << B[k]
              <<"-------------------"<< std::endl;
  }
#endif
  result.resize(nsd);
  for(size_t i=0; i<nsd; i++)
  {
    result[i].resize(n);
    std::copy(B[i].begin(), B[i].end(), result[i].begin());
  }

  return true;
}
