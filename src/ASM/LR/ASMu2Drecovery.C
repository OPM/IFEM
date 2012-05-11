// $Id$
//==============================================================================
//!
//! \file ASMu2D.C
//!
//! \date February 2012 
//!
//! \author Kjetil A. Johannessen and Knut Morten Okstad / SINTEF
//!
//! \brief Recovery techniques for unstructured LR B-splines
//!
//==============================================================================

#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Element.h"

#include "ASMu2D.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "SparseMatrix.h"
#include "DenseMatrix.h"
#include "Utilities.h"
#include "Profiler.h"
#include "IntegrandBase.h"


bool ASMu2D::getGrevilleParameters (RealArray& prm, int dir) const
{
  if (!lrspline || dir < 0 || dir > 1) return false;

  int i=0;
  prm.resize(lrspline->nBasisFunctions());
  std::vector<LR::Basisfunction*>::iterator bit;
  for(bit = lrspline->basisBegin(); bit != lrspline->basisEnd(); bit++, i++)
  	prm[i] = (**bit).getGrevilleParameter()[dir];

  return true;
}

void ASMu2D::expandTensorGrid(RealArray *in, RealArray *out) const {
  out[0].resize((in[0].size()) * (in[1].size()));
  out[1].resize((in[0].size()) * (in[1].size()));

  int ip=0;
  for(size_t j=0; j<in[1].size(); j++) {
    for(size_t i=0; i<in[0].size(); i++, ip++) {
      out[0][ip] = in[0][i];
      out[1][ip] = in[1][j];
    }
  }

}


LR::LRSplineSurface* ASMu2D::projectSolution (const IntegrandBase& integrnd) const
{
  PROFILE2("ASMu2D::projectSolution");

  // Compute parameter values of the result sampling points (Greville points)
  RealArray gpar[2];
  for (int dir = 0; dir < 2; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir))
      return 0;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrnd,gpar))
    return 0;

  // Project the results onto the spline basis to find control point
  // values based on the result values evaluated at the Greville points.
  // Note that we here implicitly assume that the number of Greville points
  // equals the number of control points such that we don't have to resize
  // the result array. Think that is always the case, but beware if trying
  // other projection schemes later.

  if (lrspline->rational()) {
    std::cerr << "Error: rational LR B-splines not supported yet\n";
    return NULL;
  }

  return regularInterpolation(lrspline,
			      gpar[0], gpar[1],
			      sValues,
			      sValues.rows());
}


LR::LRSplineSurface* ASMu2D::evalSolution (const IntegrandBase& integrand) const
{
  return this->projectSolution(integrand);
}


bool ASMu2D::globalL2projection (Matrix& sField,
				 const IntegrandBase& integrand,
				 bool continuous) const
{
  if (!lrspline) return true; // silently ignore empty patches

  PROFILE2("ASMu2D::globalL2");


/*
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();
  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;
*/
  const int p1 = lrspline->order_u();
  const int p2 = lrspline->order_v();

  // Get Gaussian quadrature points
  const int ng1 = continuous ? nGauss : p1 - 1;
  const int ng2 = continuous ? nGauss : p2 - 1;
  const double* xg = GaussQuadrature::getCoord(ng1);
  const double* yg = GaussQuadrature::getCoord(ng2);
  const double* wg = continuous ? GaussQuadrature::getWeight(nGauss) : 0;
  if (!xg || !yg) return false;
  if (continuous && !wg) return false;


  // Set up the projection matrices
  const size_t nnod = this->getNoNodes();
  const size_t ncomp = integrand.getNoFields(); 
  SparseMatrix A(SparseMatrix::SUPERLU);
  StdVector B(nnod*ncomp);
  A.redim(nnod,nnod);

  double dA = 0.0;
  Go::BasisDerivsSf spl1;
  Go::BasisPtsSf    spl0;
  Vector phi;
  Matrix dNdu, Xnod, Jac;
  int ip;


  // === Assembly loop over all elements in the patch ==========================

  std::vector<LR::Element*>::iterator el = lrspline->elementBegin();
  for (int iel = 1; el != lrspline->elementEnd(); el++, iel++)
  {
    // workingEl = iel-1;
    if (continuous)
    {
      // Set up control point (nodal) coordinates for current element
      if (!this->getElementCoordinates(Xnod,iel))
        return false;
      else if ((dA = 0.25*this->getParametricArea(iel)) < 0.0)
        return false; // topology error (probably logic error)
    }

    // Compute parameter values of the Gauss points over this element
    RealArray gpar[2], unstrGpar[2];
    for (int d = 0; d < 2; d++)
      this->getGaussPointParameters(gpar[d],d,(d==0)?ng1:ng2,iel,(d==0)?xg:yg);

    // convert to unstructred mesh representation
    expandTensorGrid(gpar, unstrGpar);

    // Evaluate the secondary solution at all integration points
    if (!this->evalSolution(sField,integrand,unstrGpar))
      return false;

    // set up basis function size (for extractBasis subroutine)
    phi.resize((**el).nBasisFunctions());

    // --- Integration loop over all Gauss points in each direction ----------
    ip = 0;
    for (int j = 0; j < ng2; j++)
      for (int i = 0; i < ng1; i++, ip++)
      {
        
        if (continuous)
	{
	  lrspline->computeBasis(gpar[0][i], gpar[1][j], spl1, iel-1);
	  extractBasis(spl1,phi,dNdu);
	}
	else
	{
	  lrspline->computeBasis(gpar[0][i], gpar[1][j], spl0, iel-1);
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
          for (size_t r = 1; r <= ncomp; r++)
            B(inod+(r-1)*nnod) += phi[ii]*sField(r,ip+1)*dJw;
        }
      }
  }
  workingEl = -1;
#if SP_DEBUG > 2
  std::cout << " ---- Matrix A -----\n";
  std::cout << A << std::endl;
  std::cout << " -------------------\n";
  std::cout << " ---- Vector B -----\n";
  std::cout << B << std::endl;
  std::cout << " -------------------\n";
#endif


  // Solve the patch-global equation system
  if (!A.solve(B)) return false;

  // Store the control-point values of the projected field
  sField.resize(ncomp,nnod);
  for (size_t i = 1; i <= nnod; i++)
    for (size_t j = 1; j <= ncomp; j++)
      sField(j,i) = B(i+(j-1)*nnod);

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

  const int p1 = lrspline->order_u();
  const int p2 = lrspline->order_v();
  // const int nel1 = surf->numCoefs_u() - p1 + 1;

  // Get Gaussian quadrature point coordinates
  const int ng1 = p1 - 1;
  const int ng2 = p2 - 1;
  const double* xg = GaussQuadrature::getCoord(ng1);
  const double* yg = GaussQuadrature::getCoord(ng2);
  if (!xg || !yg) return 0;
		  

  // Compute parameter values of the Gauss points over the whole patch
  RealArray gpar[2];

  // Compute parameter values of the Greville points
  if (!this->getGrevilleParameters(gpar[0],0)) return 0;
  if (!this->getGrevilleParameters(gpar[1],1)) return 0;

  const size_t nCmp = integrand.getNoFields(); // Number of result components
  const size_t nPol = p1*p2;                   // Number of terms in polynomial expansion

  Matrix sValues(nCmp,gpar[0].size());
  Vector P(nPol);
  Go::Point X, G;

  // Loop over all Greville points (one for each basis function)
  size_t k, l, ip = 0;
  LR::Basisfunction *b;
  std::vector<LR::Basisfunction*>::iterator bit;
  std::vector<LR::Element*>::iterator elStart, elEnd, el;
  std::vector<LR::Element*> supportElements;
  for(bit = lrspline->basisBegin(); bit != lrspline->basisEnd(); bit++, ip++)
  {
    b = *bit;
    
#if SP_DEBUG > 2
    std::cout << "Basis: " << *b  << std::endl;
    std::cout << "  nel =" << nel  << std::endl;
    std::cout << "  ng1 =" << ng1  << std::endl;
    std::cout << "  ng2 =" << ng2  << std::endl;
    std::cout << "  nPol=" << nPol << std::endl;
    std::cout << "  needs ext.basis : " <<  (nel*ng1*ng2 < nPol) << std::endl;
#endif

    // Special case for basis functions with too many zero knot spans by using
    // the extended support
    // if(nel*ng1*ng2 < nPol)
    if(true)
    {
      supportElements = b->getExtendedSupport();
      elStart = supportElements.begin();
      elEnd   = supportElements.end();
#if SP_DEBUG > 2
      std::cout << "Extended basis: \n";
      for (el = elStart; el != elEnd; el++)
	std::cout << "  " << **el << std::endl;
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
    StdVector B(nPol*nCmp);


    // Loop over all non-zero knot-spans in the support of
    // the basis function associated with current Greville point
    for (el = elStart; el != elEnd; el++)
    {
      int ig, iel = (**el).getId()+1;
      // workingEl = iel;

      // evaluate all gauss points for this element
      RealArray gaussPt[2], unstrGauss[2];
      for (int d = 0; d < 2; d++)
	this->getGaussPointParameters(gaussPt[d],d,(d==0)?ng1:ng2,iel,(d==0)?xg:yg);

#if SP_DEBUG > 2
      std::cout << "Element " << **el << std::endl;
#endif

      // convert to unstructred mesh representation
      expandTensorGrid(gaussPt, unstrGauss);

      // Evaluate the secondary solution at all Gauss points
      Matrix sField;
      if (!this->evalSolution(sField,integrand,unstrGauss))
	return 0;

      // Loop over the Gauss points in current knot-span
      ig=1;
      for(int j=0; j<ng2; j++)
	for(int i=0; i<ng1; i++, ig++)
	{
	  lrspline->point(X,gaussPt[0][i],gaussPt[1][j]);
	  // Evaluate the polynomial expansion at current Gauss point
	  evalMonomials(p1,p2,X[0]-G[0],X[1]-G[1],P);
#if SP_DEBUG > 2
	  std::cout << "Greville     point: " << G << std::endl;
	  std::cout << "Gauss        point: " << gaussPt[0][i] << ", " << gaussPt[0][j] << std::endl;
	  std::cout << "Mapped gauss point: " << X << std::endl;
	  std::cout << "P-matrix:\n " << P << "\n--------------------" << std::endl;
#endif

	  for (k = 1; k <= nPol; k++)
	  {
	    // Accumulate the projection matrix, A += P^t * P
	    for (l = 1; l <= nPol; l++)
	      A(k,l) += P(k)*P(l);

	    // Accumulate the right-hand-side matrix, B += P^t * sigma
	    for (l = 1; l <= nCmp; l++)
	      B(k+(l-1)*nPol) += P(k)*sField(l,ig);
	  }
	}
      }
      workingEl = -1;

#if SP_DEBUG > 2
    std::cout << " ---- Matrix A -----\n";
    std::cout << A << std::endl;
    std::cout << " -------------------\n";
    std::cout << " ---- Vector B -----\n";
    std::cout << B << std::endl;
    std::cout << " -------------------\n";
#endif

      // Solve the local equation system
      if (!A.solve(B)) return false;

      // Evaluate the projected field at current Greville point
      for (k = 1; k <= nCmp; k++)
	sValues(k,ip+1) = B[(k-1)*nPol];

#if SP_DEBUG > 1
      std::cout <<"Greville point "<< ip << " (u,v) = "
		<< gpar[0][ip] <<" "<< gpar[1][ip] <<" :";
      for (k = 1; k <= nCmp; k++)
	std::cout <<" "<< sValues(k,ip);
      std::cout << std::endl;
#endif
    }

  // Project the Greville point results onto the spline basis
  // to find the control point values

  if (lrspline->rational()) {
    std::cerr << "Error: rational LR B-splines not supported yet\n";
    return NULL;
  }

  return regularInterpolation(lrspline,
			      gpar[0], gpar[1],
			      sValues,
			      nCmp);
}

LR::LRSplineSurface* ASMu2D::regularInterpolation(LR::LRSplineSurface *basis,
					 const std::vector<double>& upar,
					 const std::vector<double>& vpar,
					 const Matrix& points,
					 size_t dim) const
{
  // sanity check on input parameters
  if( (upar.size() != vpar.size())   ||
      (upar.size() != points.cols()) ||
      (dim         != points.rows()) ) {
    std::cerr << "Error in ASMu2D::regularInterpolation() - mismatching input vector/matrix sizes " << std::endl;
    return NULL;
  }

  size_t nBasis = basis->nBasisFunctions();
  Go::BasisPtsSf splineValues; 

  // copy all basis functions and mesh and later swap around the control points
  LR::LRSplineSurface *ans = basis->copy();
  ans->rebuildDimension(dim);

  DenseMatrix A(nBasis, nBasis);

  // evaluate all basis functions at all points in A matrix (same row = same evaluation point)
  for(size_t i=0; i<upar.size(); i++) {
    lrspline->computeBasis(upar[i], vpar[i], splineValues);
    // optimization note: without an element id, splineValues will be stored as a full dense vector
    std::vector<double> N = splineValues.basisValues;
    for(size_t j=0; j<nBasis; j++)
      A(i+1,j+1) = N[j];
  }

  // for all solution components, build a right-hand-side B
  for(size_t i=0; i<dim; i++) {
    // optimization note: we are basically unwrapping the matrix "points" into several column vectors.
    //                    would be better to solve for all right hand sides at once (reuse LU decomposition
    //                    and don't copy-paste so much data around)
    StdVector B(points.cols());
    for(size_t j=1; j<=points.cols(); j++)
      B(j) = points(i+1,j);

    if(!A.solve(B))
      return NULL;

    for(size_t j=0; j<nBasis; j++)
      ans->getBasisfunction(j)->controlpoint_[i] = B(j+1);
  }

  return ans;
    
}

