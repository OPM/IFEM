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
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "SparseMatrix.h"
#include "DenseMatrix.h"
#include "SplineUtils.h"
#include "Profiler.h"
#include <array>
#include <fstream>


bool ASMu2D::getGrevilleParameters (RealArray& prm, int dir, int basisNum) const
{
  if (!this->getBasis(basisNum) || dir < 0 || dir > 1) return false;

  const LR::LRSpline* lrspline = this->getBasis(basisNum);

  prm.clear();
  prm.reserve(lrspline->nBasisFunctions());

  for (const LR::Basisfunction* b : lrspline->getAllBasisfunctions())
    prm.push_back(b->getGrevilleParameter()[dir]);

  return true;
}


/*!
  \brief Expands a tensor parametrization point to an unstructured one.
  \details Takes as input a tensor mesh, for instance
     in[0] = {0,1,2}
     in[1] = {2,3,5}
   and expands this to an unstructured representation, i.e.,
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
    if (!this->getGrevilleParameters(gpar[dir],dir,1))
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

  return this->regularInterpolation(gpar[0],gpar[1],sValues,1);
}


LR::LRSpline* ASMu2D::evalSolution (const IntegrandBase& integrand) const
{
  return this->projectSolution(integrand);
}


bool ASMu2D::assembleL2matrices (SparseMatrix& A, StdVector& B,
                                 const IntegrandBase& integrand,
                                 bool continuous) const
{
  size_t nnod = this->getNoProjectionNodes();

  const int p1 = projBasis->order(0);
  const int p2 = projBasis->order(1);

  // Get Gaussian quadrature points
  const int ng1 = continuous ? nGauss : p1 - 1;
  const int ng2 = continuous ? nGauss : p2 - 1;
  const double* xg = GaussQuadrature::getCoord(ng1);
  const double* yg = GaussQuadrature::getCoord(ng2);
  const double* wg = continuous ? GaussQuadrature::getWeight(nGauss) : nullptr;
  if (!xg || !yg) return false;
  if (continuous && !wg) return false;

  double dA = 0.0;
  Vector phi, phi2;
  Matrix dNdu, Xnod, Jac;
  Go::BasisPtsSf    spl0;
  Go::BasisDerivsSf spl1, spl2;


  // === Assembly loop over all elements in the patch ==========================

  for (const LR::Element* el1 : lrspline->getAllElements())
  {
    double uh = (el1->umin()+el1->umax())/2.0;
    double vh = (el1->vmin()+el1->vmax())/2.0;
    int ielp = projBasis->getElementContaining(uh,vh);
    int iel = lrspline->getElementContaining(uh,vh)+1;

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
    expandTensorGrid(gpar.data(),unstrGpar.data());

    // Evaluate the secondary solution at all integration points
    Matrix sField;
    if (!this->evalSolution(sField,integrand,unstrGpar.data()))
      return false;

    // Set up basis function size (for extractBasis subroutine)
    const LR::Element* elm = projBasis->getElement(ielp);
    size_t nbf = elm->nBasisFunctions();

    IntVec lmnpc;
    if (projBasis != lrspline)
    {
      lmnpc.reserve(nbf);
      for (const LR::Basisfunction* f : elm->support())
        lmnpc.push_back(f->getId());
    }
    const IntVec& mnpc = projBasis == lrspline ? MNPC[iel-1] : lmnpc;

    // --- Integration loop over all Gauss points in each direction ----------

    Matrix eA(nbf, nbf);
    Vectors eB(sField.rows(), Vector(nbf));
    int ip = 0;
    for (int j = 0; j < ng2; j++)
      for (int i = 0; i < ng1; i++, ip++)
      {
        if (continuous)
        {
          this->computeBasis(gpar[0][i],gpar[1][j],spl1,ielp,projBasis.get());
          SplineUtils::extractBasis(spl1,phi,dNdu);
          this->computeBasis(gpar[0][i],gpar[1][j],spl2,iel-1);
          SplineUtils::extractBasis(spl2,phi2,dNdu);
        }
        else
        {
          this->computeBasis(gpar[0][i],gpar[1][j],spl0,ielp,projBasis.get());
          phi = spl0.basisValues;
        }

        // Compute the Jacobian inverse and derivatives
        double dJw = 1.0;
        if (continuous)
        {
          dJw = dA*wg[i]*wg[j]*utl::Jacobian(Jac,dNdu,Xnod,dNdu,false);
          if (dJw == 0.0) continue; // skip singular points
        }

        // Integrate the mass matrix
        eA.outer_product(phi, phi, true, dJw);

        // Integrate the rhs vector B
        for (size_t r = 1; r <= sField.rows(); r++)
          eB[r-1].add(phi,sField(r,ip+1)*dJw);
      }

    for (size_t i = 0; i < eA.rows(); ++i) {
      for (size_t j = 0; j < eA.cols(); ++j)
        A(mnpc[i]+1, mnpc[j]+1) += eA(i+1,j+1);

      int jp = mnpc[i]+1;
      for (size_t r = 0; r < sField.rows(); r++, jp += nnod)
        B(jp) += eB[r](1+i);
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
  if (!this->getGrevilleParameters(gpar[0],0,1)) return nullptr;
  if (!this->getGrevilleParameters(gpar[1],1,1)) return nullptr;

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
  for (LR::Basisfunction* b : lrspline->getAllBasisfunctions())
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
    for (el = elStart; el != elEnd; ++el)
    {
      int iel = (**el).getId()+1;
#if SP_DEBUG > 2
      std::cout <<"Element "<< **el << std::endl;
#endif

      // Compute parameter values of the Gauss points over this element
      std::array<RealArray,2> gaussPt, unstrGauss;
      this->getGaussPointParameters(gaussPt[0],0,ng1,iel,xg);
      this->getGaussPointParameters(gaussPt[1],1,ng2,iel,yg);
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

  return this->regularInterpolation(gpar[0],gpar[1],sValues,1);
}


LR::LRSplineSurface* ASMu2D::regularInterpolation (const RealArray& upar,
                                                   const RealArray& vpar,
                                                   const Matrix& points,
                                                   int basis) const
{
  const LR::LRSplineSurface* lrspline = this->getBasis(basis);
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

  SparseMatrix A(SparseMatrix::SUPERLU);
  A.resize(nBasis, nBasis);
  Matrix B2(points,true); // transpose to get one vector per field
  StdVector B(B2);
  Go::BasisPtsSf splineValues;

  // Evaluate all basis functions at all points, stored in the A-matrix
  // (same row = same evaluation point)
  for (size_t i = 0; i < nBasis; i++)
  {
    size_t k = 0;
    int iel = lrspline->getElementContaining(upar[i], vpar[i]);
    this->computeBasis(upar[i],vpar[i],splineValues,iel);
    for (LR::Basisfunction* function : lrspline->getElement(iel)->support())
      A(i+1,function->getId()+1) = splineValues.basisValues[k++];
  }

  // Solve for all solution components - one right-hand-side for each
  if (!A.solve(B)) return nullptr;

  // Copy all basis functions and mesh
  LR::LRSplineSurface* ans = lrspline->copy();
  ans->rebuildDimension(points.rows());

  // Back to interleaved data
  RealArray interleave;
  interleave.reserve(B.dim());
  for (size_t i = 0; i < nBasis; i++)
    for (size_t j = 0; j < points.rows(); j++)
      interleave.push_back(B(1+j*points.cols()+i));

  ans->setControlPoints(interleave);

  return ans;
}


bool ASMu2D::edgeL2projection (const DirichletEdge& edge,
                               const FunctionBase& values,
                               Real2DMat& result,
                               double time) const
{
  size_t n = edge.MLGN.size();
  size_t m = values.dim();
  SparseMatrix A(SparseMatrix::SUPERLU);
  StdVector B(n*m);
  A.resize(n,n);

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // find the normal and tangent direction for the edge
  int edgeDir, t1, t2;
  switch (edge.edg)
  { //              id            normal  tangent
    case LR::WEST:  edgeDir = -1; t1 = 1; t2 = 2; break;
    case LR::EAST:  edgeDir =  1; t1 = 1; t2 = 2; break;
    case LR::SOUTH: edgeDir = -2; t1 = 2; t2 = 1; break;
    case LR::NORTH: edgeDir =  2; t1 = 2; t2 = 1; break;
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

  Vector N;
  Matrix dNdu, dNdX, Xnod, Jac;
  double param[3] = { 0.0, 0.0, 0.0 };
  Vec4   X(param);

  // === Assembly loop over all elements on the patch edge =====================

  for (size_t i = 0; i < edge.MLGE.size(); i++) // for all edge elements
  {
    int iel = 1 + edge.MLGE[i];

    // Get element edge length in the parameter space
    double dS = 0.5*this->getParametricLength(iel,t2);
    if (dS < 0.0) return false; // topology error (probably logic error)

    // Set up control point coordinates for current element
    if (!this->getElementCoordinates(Xnod,iel)) return false;

    // Get integration gauss points over this element
    this->getGaussPointParameters(gpar[t2-1],t2-1,nGauss,iel,xg);

    // --- Integration loop over all Gauss points along the edge ---------------

    for (int j = 0; j < nGauss; j++)
    {
      // Parameter values of current integration point
      double u = param[0] = gpar[0][j];
      double v = param[1] = gpar[1][j];

      // Evaluate basis function derivatives at current integration points
      Go::BasisDerivsSf spline;
      this->computeBasis(u,v,spline,iel-1);

      // Fetch basis function derivatives at current integration point
      SplineUtils::extractBasis(spline,N,dNdu);

      // Compute basis function derivatives
      double detJxW = dS*utl::Jacobian(Jac,X,dNdX,Xnod,dNdu,t1,t2)*wg[j];
      if (detJxW == 0.0) continue; // skip singular points

      // Cartesian coordinates of current integration point
      X.assign(Xnod * N);
      X.t = time;

      // For mixed basis, we need to compute functions separate from geometry
      if (edge.lr != lrspline.get())
      {
        // different lrspline instances enumerate elements differently
        edge.lr->computeBasis(u,v,spline,edge.lr->getElementContaining(u,v));
        SplineUtils::extractBasis(spline,N,dNdu);
      }

      // Assemble into matrix A and vector B
      for (size_t il = 0; il < edge.MNPC[i].size(); il++) { // local i-index
        int ig;
        if ((ig = 1+edge.MNPC[i][il]) > 0)         // global i-index
        {
          for (size_t jl = 0; jl < edge.MNPC[i].size(); jl++) { // local j-index
            int jg;
            if ((jg = 1+edge.MNPC[i][jl]) > 0)         // global j-index
              A(ig,jg) += N[il]*N[jl]*detJxW;
          }

          RealArray val = values.getValue(X);
          for (size_t k = 0; k < m; k++)
            B(ig+k*n) += N[il]*val[k]*detJxW;
        } // end basis-function loop
      }
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

  std::cout <<"---- Edge-nodes (g2l-mapping) -----";
  for (size_t i = 0; i < edge.MLGN.size(); i++)
    std::cout <<"\n     "<< i <<": "<< edge.MLGN[i];
  std::cout <<"\n-------------------";
  std::cout <<"\n---- Element-nodes (-1 is interior element nodes) -----";
  for (const IntVec& d : edge.MNPC) {
    std::cout <<"\n    ";
    for (int c : d)
      std::cout <<" "<< c;
  }
  std::cout <<"\n-------------------"<< std::endl;
#endif

  // Solve the edge-global equation system
  if (!A.solve(B)) return false;

#if SP_DEBUG > 2
  std::cout <<"---- SOLUTION -----\n"<< B
            <<"-------------------"<< std::endl;
#endif

  // Store the control-point values of the projected field
  result.resize(m,RealArray(n));
  RealArray::const_iterator it = B.begin();
  for (size_t i = 0; i < m; i++, it += n)
    std::copy(it, it+n, result[i].begin());

  return true;
}
