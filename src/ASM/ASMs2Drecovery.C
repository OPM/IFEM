// $Id$
//==============================================================================
//!
//! \file ASMs2Drecovery.C
//!
//! \date Feb 20 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Recovery of secondary solutions for structured 2D spline FE models.
//!
//==============================================================================

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SurfaceInterpolator.h"

#include "ASMs2D.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "SparseMatrix.h"
#include "DenseMatrix.h"
#include "Utilities.h"
#include "Profiler.h"


bool ASMs2D::getGrevilleParameters (RealArray& prm, int dir) const
{
  if (!surf || dir < 0 || dir > 1) return false;

  const Go::BsplineBasis& basis = surf->basis(dir);

  prm.resize(basis.numCoefs());
  for (size_t i = 0; i < prm.size(); i++)
    prm[i] = basis.grevilleParameter(i);

  return true;
}


bool ASMs2D::getQuasiInterplParameters (RealArray& prm, int dir) const
{
  if (!surf || dir < 0 || dir > 1) return false;

  const Go::BsplineBasis& basis = surf->basis(dir);

  RealArray knots_simple;
  basis.knotsSimple(knots_simple);

  prm.clear();
  for (size_t i = 0; i+1 < knots_simple.size(); i++)
  {
    prm.push_back(knots_simple[i]);
    prm.push_back(0.5*(knots_simple[i]+knots_simple[i+1]));
  }
  prm.push_back(knots_simple.back());

  return true;
}


bool ASMs2D::evaluate (const ASMbase* basis, const Vector& locVec,
		       Vector& vec) const
{
  const ASMs2D* pch = dynamic_cast<const ASMs2D*>(basis);
  if (!pch) return false;

  // Compute parameter values of the result sampling points (Greville points)
  RealArray gpar[2];
  for (int dir = 0; dir < 2; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir))
      return false;

  // Evaluate the result field at all sampling points.
  // Note: it is here assumed that *basis and *this have spline bases
  // defined over the same parameter domain.
  Matrix sValues;
  if (!pch->evalSolution(sValues,locVec,gpar))
    return false;

  // Project the results onto the spline basis to find control point
  // values based on the result values evaluated at the Greville points.
  // Note that we here implicitly assume that the number of Greville points
  // equals the number of control points such that we don't have to resize
  // the result array. Think that is always the case, but beware if trying
  // other projection schemes later.

  RealArray weights;
  if (surf->rational())
    surf->getWeights(weights);

  const Vector& vec2 = sValues;
  Go::SplineSurface* surf_new =
    Go::SurfaceInterpolator::regularInterpolation(surf->basis(0),
						  surf->basis(1),
						  gpar[0], gpar[1],
						  const_cast<Vector&>(vec2),
						  sValues.rows(),
						  surf->rational(),
						  weights);

  vec.resize(surf_new->coefs_end()-surf_new->coefs_begin());
  std::copy(surf_new->coefs_begin(),surf_new->coefs_end(),vec.begin());
  delete surf_new;

  return true;
}


Go::SplineSurface* ASMs2D::projectSolution (const IntegrandBase& integrnd) const
{
  PROFILE2("ASMs2D::projectSolution");

  // Compute parameter values of the result sampling points (Greville points)
  RealArray gpar[2];
  for (int dir = 0; dir < 2; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir))
      return 0;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrnd,gpar) || sValues.rows() == 0)
    return 0;

  // Project the results onto the spline basis to find control point
  // values based on the result values evaluated at the Greville points.
  // Note that we here implicitly assume that the number of Greville points
  // equals the number of control points such that we don't have to resize
  // the result array. Think that is always the case, but beware if trying
  // other projection schemes later.

  RealArray weights;
  if (surf->rational())
    surf->getWeights(weights);

  const Vector& vec = sValues;
  return Go::SurfaceInterpolator::regularInterpolation(surf->basis(0),
						       surf->basis(1),
						       gpar[0], gpar[1],
						       const_cast<Vector&>(vec),
						       sValues.rows(),
						       surf->rational(),
						       weights);
}


Go::GeomObject* ASMs2D::evalSolution (const IntegrandBase& integrand) const
{
  return this->projectSolution(integrand);
}


bool ASMs2D::globalL2projection (Matrix& sField,
				 const IntegrandBase& integrand,
				 bool continuous) const
{
  if (!surf) return true; // silently ignore empty patches

  PROFILE2("ASMs2D::globalL2");

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();
  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;

  // Get Gaussian quadrature point coordinates (and weights if continuous)
  const int ng1 = continuous ? nGauss : p1 - 1;
  const int ng2 = continuous ? nGauss : p2 - 1;
  const double* xg = GaussQuadrature::getCoord(ng1);
  const double* yg = GaussQuadrature::getCoord(ng2);
  const double* wg = continuous ? GaussQuadrature::getWeight(nGauss) : 0;
  if (!xg || !yg) return false;
  if (continuous && !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  Matrix gp;
  RealArray gpar[2];
  gpar[0] = this->getGaussPointParameters(gp,0,ng1,xg);
  gpar[1] = this->getGaussPointParameters(gp,1,ng2,yg);

  // Evaluate basis functions at all integration points
  std::vector<Go::BasisPtsSf>    spl0;
  std::vector<Go::BasisDerivsSf> spl1;
  if (continuous)
    surf->computeBasisGrid(gpar[0],gpar[1],spl1);
  else
    surf->computeBasisGrid(gpar[0],gpar[1],spl0);

  // Evaluate the secondary solution at all integration points
  if (!this->evalSolution(sField,integrand,gpar))
  {
    std::cerr <<" *** ASMs2D::globalL2projection: Failed for patch "<< idx+1
	      <<" nPoints="<< gpar[0].size()*gpar[1].size() << std::endl;
    return false;
  }

  // Set up the projection matrices
  const size_t nnod = this->getNoNodes(1);
  const size_t ncomp = sField.rows();
  SparseMatrix A(SparseMatrix::SUPERLU);
  StdVector B(nnod*ncomp);
  A.redim(nnod,nnod);

  double dA = 1.0;
  Vector phi(p1*p2);
  Matrix dNdu, Xnod, J;


  // === Assembly loop over all elements in the patch ==========================

  int iel = 0;
  for (int i2 = 0; i2 < nel2; i2++)
    for (int i1 = 0; i1 < nel1; i1++, iel++)
    {
      if (MLGE[iel] < 1) continue; // zero-area element

      if (continuous)
      {
	// Set up control point (nodal) coordinates for current element
	if (!this->getElementCoordinates(Xnod,1+iel))
	  return false;
	else if ((dA = 0.25*this->getParametricArea(1+iel)) < 0.0)
	  return false; // topology error (probably logic error)
      }

      // --- Integration loop over all Gauss points in each direction ----------

      int ip = (i2*ng1*nel1 + i1)*ng2;
      for (int j = 0; j < ng2; j++, ip += ng1*(nel1-1))
	for (int i = 0; i < ng1; i++, ip++)
	{
	  if (continuous)
	    extractBasis(spl1[ip],phi,dNdu);
	  else
	    phi = spl0[ip].basisValues;

	  // Compute the Jacobian inverse and derivatives
	  double dJw = 1.0;
	  if (continuous)
	  {
	    dJw = dA*wg[i]*wg[j]*utl::Jacobian(J,dNdu,Xnod,dNdu,false);
	    if (dJw == 0.0) continue; // skip singular points
	  }

	  // Integrate the linear system A*x=B
	  for (size_t ii = 0; ii < phi.size(); ii++)
	  {
	    int inod = MNPC[iel][ii]+1;
	    for (size_t jj = 0; jj < phi.size(); jj++)
	    {
	      int jnod = MNPC[iel][jj]+1;
	      A(inod,jnod) += phi[ii]*phi[jj]*dJw;
	    }
	    for (size_t r = 1; r <= ncomp; r++)
	      B(inod+(r-1)*nnod) += phi[ii]*sField(r,ip+1)*dJw;
	  }
	}
    }

#if SP_DEBUG > 1
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


Go::SplineSurface* ASMs2D::scRecovery (const IntegrandBase& integrand) const
{
  PROFILE2("ASMs2D::scRecovery");

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int nel1 = surf->numCoefs_u() - p1 + 1;

  // Get Gaussian quadrature point coordinates
  const int ng1 = p1 - 1;
  const int ng2 = p2 - 1;
  const double* xg = GaussQuadrature::getCoord(ng1);
  const double* yg = GaussQuadrature::getCoord(ng2);
  if (!xg || !yg) return 0;

  // Compute parameter values of the Gauss points over the whole patch
  Matrix gaussPt[2];
  RealArray gpar[2];
  gpar[0] = this->getGaussPointParameters(gaussPt[0],0,ng1,xg);
  gpar[1] = this->getGaussPointParameters(gaussPt[1],1,ng2,yg);

#if SP_DEBUG > 2
  for (size_t j = 0; j < gpar[1].size(); j++)
    for (size_t i = 0; i < gpar[0].size(); i++)
      std::cout <<"Gassian point "<< i <<","<< j <<" = "
		<< gpar[0][i] <<" "<< gpar[1][j] << std::endl;
#endif

  // Evaluate the secondary solution at all Gauss points
  Matrix sField;
  if (!this->evalSolution(sField,integrand,gpar))
    return 0;

  // Compute parameter values of the Greville points
  if (!this->getGrevilleParameters(gpar[0],0)) return 0;
  if (!this->getGrevilleParameters(gpar[1],1)) return 0;

  const size_t nCmp = sField.rows(); // Number of result components
  const size_t nPol = (p1+1)*(p2+1); // Number of terms in polynomial expansion

  Matrix sValues(nCmp,gpar[0].size()*gpar[1].size());
  Vector P(nPol);
  Go::Point X, G;

  // Loop over all Greville points
  int istart, jstart = 0;
  size_t ig, jg, k, l, ip = 1;
  for (jg = 0; jg < gpar[1].size(); jg++)
    for (ig = 0; ig < gpar[0].size(); ig++, ip++)
    {
      // Special case for the first and last Greville points in each direction
      if (ig == 0)
	istart = 1;
      else if (ig < gpar[0].size()-1)
	istart = ig;
      if (jg == 0)
	jstart = 1;
      else if (jg < gpar[1].size()-1)
	jstart = jg;

      // Physical coordinates of current Greville point
      surf->point(G,gpar[0][ig],gpar[1][jg]);

      // Set up the local projection matrices
      DenseMatrix A(nPol,nPol,true);
      Matrix B(nPol,nCmp);

      // Loop over all non-zero knot-spans in the support of
      // the basis function associated with current Greville point
      for (int js = jstart; js < jstart+p2; js++)
	if (js >= p2-1 && surf->knotSpan(1,js) > 0.0)
	  for (int is = istart; is < istart+p1; is++)
	    if (is >= p1-1 && surf->knotSpan(0,is) > 0.0)
	    {
	      // Loop over the Gauss points in current knot-span
	      int jp = ((js-p2+1)*ng1*nel1 + is-p1+1)*ng2 + 1;
	      for (int j = 1; j <= ng2; j++, jp += ng1*(nel1-1))
		for (int i = 1; i <= ng1; i++, jp++)
		{
		  // Fetch parameter values of current integration point
		  double u = gaussPt[0](i,is-p1+2);
		  double v = gaussPt[1](j,js-p2+2);

		  // Evaluate the polynomial expansion at current Gauss point
		  surf->point(X,u,v);
		  evalMonomials(p1+1,p2+1,X[0]-G[0],X[1]-G[1],P);

		  for (k = 1; k <= nPol; k++)
		  {
		    // Accumulate the projection matrix, A += P^t * P
		    for (l = k; l <= nPol; l++) // do the upper triangle only
		      A(k,l) += P(k)*P(l);

		    // Accumulate the right-hand-side matrix, B += P^t * sigma
		    for (l = 1; l <= nCmp; l++)
		      B(k,l) += P(k)*sField(l,jp);
		  }
		}
	    }

      // Solve the local equation system
      if (!A.solve(B)) return false;

      // Evaluate the projected field at current Greville point (first row of B)
      for (l = 1; l <= nCmp; l++)
	sValues(l,ip) = B(1,l);

#if SP_DEBUG > 1
      std::cout <<"Greville point "<< ig <<","<< jg <<" (u,v) = "
		<< gpar[0][ig] <<" "<< gpar[1][jg] <<" :";
      for (k = 1; k <= nCmp; k++)
	std::cout <<" "<< sValues(k,ip);
      std::cout << std::endl;
#endif
    }

  // Project the Greville point results onto the spline basis
  // to find the control point values

  RealArray weights;
  if (surf->rational())
    surf->getWeights(weights);

  const Vector& vec = sValues;
  return Go::SurfaceInterpolator::regularInterpolation(surf->basis(0),
						       surf->basis(1),
						       gpar[0], gpar[1],
						       const_cast<Vector&>(vec),
						       nCmp, surf->rational(),
						       weights);
}


#include "ASMs2DInterpolate.C" // TODO: inline these methods instead...

// L2-Projection: Least-square approximation; global approximation
Go::SplineSurface* ASMs2D::projectSolutionLeastSquare (const IntegrandBase& integrand) const
{
  if (!surf) return false;

  // Compute parameter values of the result sampling points (Gauss-Interpl. points)
  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  Matrix ggpar[2];
  for (int dir = 0; dir < 2; dir++)
    this->getGaussPointParameters(ggpar[dir],dir,nGauss,xg);

  std::vector<double> par_u;
  par_u = ggpar[0];
  std::vector<double> par_v;
  par_v = ggpar[1];

  // gauss weights at parameter values
  std::vector<double> wgpar_u;
  std::vector<double> wgpar_v;
  for (int dir = 0; dir < 2; dir++){
    const Go::BsplineBasis& basis = surf->basis(dir);
    RealArray::const_iterator knotit = basis.begin();
    std::vector<double> tmp;
    tmp.reserve(nGauss*(basis.numCoefs()-basis.order()));
    for (size_t i = 0; i<=(basis.numCoefs()-basis.order());i++)
    {
      double d = knotit[i+basis.order()]-knotit[i+basis.order()-1];
      for (int j = 0; j < nGauss; j++)
        tmp.push_back(d > 0.0 ? wg[j]*d*0.5 : 0.0);
    }
    if (dir == 0)
      wgpar_u = tmp;
    else if (dir == 1)
      wgpar_v = tmp;
  }

  RealArray gpar[2];
  gpar[0] = par_u;
  gpar[1] = par_v;

  // Evaluate the secondary solution at all sampling points (Gauss points)
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar))
    return 0;

  RealArray weights;
  if (surf->rational())
    surf->getWeights(weights);

  return leastsquare_approximation(surf->basis(0),
				   surf->basis(1),
				   par_u, par_v,
				   wgpar_u, wgpar_v,
				   sValues,
				   sValues.rows(),
				   surf->rational(),
				   weights);
}


// Quasi-Interpolation; local interpolation method
Go::SplineSurface* ASMs2D::projectSolutionLocal (const IntegrandBase& integrand) const
{
  // Secondary solution evaluated at Quasi-Interpolation points
  // Compute parameter values of the result sampling points (Quasi-Interpl. points)

  RealArray gpar[2];
  for (int dir = 0; dir < 2; dir++)
    if (!this->getQuasiInterplParameters(gpar[dir],dir))
      return 0;

  // Evaluate the secondary solution at all sampling points (Quasi-Interpl. points)
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar))
    return 0;

  RealArray weights;
  if (surf->rational())
    surf->getWeights(weights);

  return quasiInterpolation(surf->basis(0),
			    surf->basis(1),
			    gpar[0], gpar[1],
			    sValues,
			    sValues.rows(),
			    surf->rational(),
			    weights);
}


// Variation Diminishing Spline Approximation (VDSA); local approx. method
Go::SplineSurface* ASMs2D::projectSolutionLocalApprox(const IntegrandBase& integrand) const
{
  // Compute parameter values of the result sampling points (Greville points)
  RealArray gpar[2];
  for (int dir = 0; dir < 2; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir))
      return 0;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar))
    return 0;

  RealArray weights;
  if (surf->rational())
    surf->getWeights(weights);

  return VariationDiminishingSplineApproximation(surf->basis(0),
						 surf->basis(1),
						 gpar[0], gpar[1],
						 sValues,
						 sValues.rows(),
						 surf->rational(),
						 weights);
}
