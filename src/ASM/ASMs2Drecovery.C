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
#include "ItgPoint.h"
#include "Field.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "GlbL2projector.h"
#include "SparseMatrix.h"
#include "DenseMatrix.h"
#include "SplineUtils.h"
#include "Profiler.h"
#include <array>


bool ASMs2D::getGrevilleParameters (RealArray& prm, int dir, int basisNum) const
{
  if (dir < 0 || dir > 1) return false;

  const Go::BsplineBasis& basis = this->getBasis(basisNum)->basis(dir);

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
		       RealArray& vec, int basisNum) const
{
  const ASMs2D* pch = dynamic_cast<const ASMs2D*>(basis);
  if (!pch) return false;

  // Compute parameter values of the result sampling points (Greville points)
  std::array<RealArray,2> gpar;
  for (int dir = 0; dir < 2; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir,basisNum))
      return false;

  // Evaluate the result field at all sampling points.
  // Note: it is here assumed that *basis and *this have spline bases
  // defined over the same parameter domain.
  Matrix sValues;
  if (!pch->evalSolution(sValues,locVec,gpar.data()))
    return false;

  const Go::SplineSurface* surf = this->getBasis(basisNum);

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

  vec.assign(surf_new->coefs_begin(),surf_new->coefs_end());
  delete surf_new;

  return true;
}


Go::SplineSurface* ASMs2D::projectSolution (const IntegrandBase& integrnd) const
{
  PROFILE2("ASMs2D::projectSolution");

  const int basis = 1;

  // Compute parameter values of the result sampling points (Greville points)
  std::array<RealArray,2> gpar;
  for (int dir = 0; dir < 2; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir,basis))
      return nullptr;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrnd,gpar.data()) || sValues.rows() == 0)
    return nullptr;

  const Go::SplineSurface* psurf = this->getBasis(basis);

  // Project the results onto the spline basis to find control point
  // values based on the result values evaluated at the Greville points.
  // Note that we here implicitly assume that the number of Greville points
  // equals the number of control points such that we don't have to resize
  // the result array. Think that is always the case, but beware if trying
  // other projection schemes later.

  RealArray weights;
  if (psurf->rational())
    psurf->getWeights(weights);

  const Vector& vec = sValues;
  return Go::SurfaceInterpolator::regularInterpolation(psurf->basis(0),
                                                       psurf->basis(1),
                                                       gpar[0], gpar[1],
                                                       const_cast<Vector&>(vec),
                                                       sValues.rows(),
                                                       psurf->rational(),
                                                       weights);
}


Go::GeomObject* ASMs2D::evalSolution (const IntegrandBase& integrand) const
{
  return this->projectSolution(integrand);
}


bool ASMs2D::assembleL2matrices (SparseMatrix& A, StdVector& B,
                                 const L2Integrand& integrand,
                                 bool continuous) const
{
  const size_t nnod = this->getNoProjectionNodes();

  const Go::SplineSurface* geo = this->getBasis(ASM::GEOMETRY_BASIS);
  const Go::SplineSurface* proj = this->getBasis(ASM::PROJECTION_BASIS);
  const bool separateProjBasis = proj != geo;
  const bool singleBasis = !separateProjBasis && this->getNoBasis() == 1;

  const int p1 = proj->order_u();
  const int p2 = proj->order_v();
  const int n1 = proj->numCoefs_u();
  int nel1 = proj->numCoefs_u() - p1 + 1;
  int nel2 = proj->numCoefs_v() - p2 + 1;

  const int pmax = p1 > p2 ? p1 : p2;

  // Get Gaussian quadrature point coordinates (and weights if continuous)
  const int ng1 = continuous ? this->getNoGaussPt(pmax,true) : p1 - 1;
  const int ng2 = continuous ? ng1 : p2 - 1;
  const double* xg = GaussQuadrature::getCoord(ng1);
  const double* yg = GaussQuadrature::getCoord(ng2);
  const double* wg = continuous ? GaussQuadrature::getWeight(ng1) : nullptr;
  if (!xg || !yg) return false;
  if (continuous && !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  Matrix gp;
  std::array<RealArray,2> gpar;
  gpar[0] = this->getGaussPointParameters(gp,0,ng1,xg,proj);
  gpar[1] = this->getGaussPointParameters(gp,1,ng2,yg,proj);

  // Evaluate basis functions at all integration points
  std::vector<Go::BasisPtsSf>    spl1;
  std::vector<Go::BasisDerivsSf> spl2;
  if (continuous)
    geo->computeBasisGrid(gpar[0],gpar[1],spl2);

  if (!continuous || separateProjBasis)
    proj->computeBasisGrid(gpar[0],gpar[1],spl1);

  // Evaluate the secondary solution at all integration points
  Matrix sField;
  if (!integrand.evaluate(sField,gpar.data()))
  {
    std::cerr <<" *** ASMs2D::assembleL2matrices: Failed for patch "<< idx+1
              <<" nPoints="<< gpar[0].size()*gpar[1].size() << std::endl;
    return false;
  }

  double dA = 1.0;
  Vector phi(p1*p2);
  Matrix dNdu, Xnod, J;


  // === Assembly loop over all elements in the patch ==========================

  int iel = 0;
  for (int i2 = 0; i2 < nel2; i2++)
    for (int i1 = 0; i1 < nel1; i1++, iel++)
    {
      int ip = (i2*ng1*nel1 + i1)*ng2;
      IntVec lmnpc;
      if (!singleBasis && proj->knotSpan(0,i1+p1-1) > 0.0
                       && proj->knotSpan(1,i2+p2-1) > 0.0)
      {
        // Establish nodal point correspondance for the projection element
        int vidx;
        lmnpc.reserve(phi.size());
        if (separateProjBasis)
          vidx = (spl1[ip].left_idx[1]-p1+1)*n1 + (spl1[ip].left_idx[0]-p1+1);
        else
          vidx = (spl2[ip].left_idx[1]-p1+1)*n1 + (spl2[ip].left_idx[0]-p1+1);
        for (int j = 0; j < p2; j++, vidx += n1)
          for (int i = 0; i < p1; i++)
            lmnpc.push_back(vidx+i);
      }
      const IntVec& mnpc = singleBasis ? MNPC[iel] : lmnpc;
      if (mnpc.empty())
        continue;

      if (continuous)
      {
        // Set up control point (nodal) coordinates for current element
        if (singleBasis)
        {
          if (!this->getElementCoordinates(Xnod,1+iel))
            return false;
          else if ((dA = this->getParametricArea(1+iel)) < 0.0)
            return false; // topology error (probably logic error)
        }
        else
        {
          if (!this->getElementCoordinatesPrm(Xnod,gpar[0][i1*ng1],gpar[1][i2*ng2]))
            return false;
          else if ((dA = 0.25 * proj->knotSpan(0,mnpc.back() % n1)
                              * proj->knotSpan(1,mnpc.back() / n1)) < 0.0)
            return false; // topology error (probably logic error)
        }
      }

      // --- Integration loop over all Gauss points in each direction ----------

      Matrix eA(p1*p2, p1*p2);
      Vectors eB(sField.rows(), Vector(p1*p2));
      for (int j = 0; j < ng2; j++, ip += ng1*(nel1-1))
        for (int i = 0; i < ng1; i++, ip++)
        {
          if (continuous)
            SplineUtils::extractBasis(spl2[ip],phi,dNdu);

          if (!continuous || separateProjBasis)
            phi = spl1[ip].basisValues;

          // Compute the Jacobian inverse and derivatives
          double dJw = 1.0;
          if (continuous)
          {
            dJw = dA*wg[i]*wg[j]*utl::Jacobian(J,dNdu,Xnod,dNdu,false);
            if (dJw == 0.0) continue; // skip singular points
          }

          // Integrate the mass matrix
          eA.outer_product(phi, phi, true, dJw);

          // Integrate the rhs vector B
          for (size_t r = 1; r <= sField.rows(); r++)
            eB[r-1].add(phi,sField(r,ip+1)*dJw);
        }

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

  const int m = integrand.derivativeOrder();
  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int nel1 = surf->numCoefs_u() - p1 + 1;

  // Get Gaussian quadrature point coordinates
  const int ng1 = p1 - m;
  const int ng2 = p2 - m;
  const double* xg = GaussQuadrature::getCoord(ng1);
  const double* yg = GaussQuadrature::getCoord(ng2);
  if (!xg || !yg) return nullptr;

  // Compute parameter values of the Gauss points over the whole patch
  std::array<Matrix,2> gaussPt;
  std::array<RealArray,2> gpar;
  gpar[0] = this->getGaussPointParameters(gaussPt[0],0,ng1,xg);
  gpar[1] = this->getGaussPointParameters(gaussPt[1],1,ng2,yg);

#if SP_DEBUG > 2
  for (size_t j = 0; j < gpar[1].size(); j++)
    for (size_t i = 0; i < gpar[0].size(); i++)
      std::cout <<"Gauss point "<< i <<","<< j <<" = "
		<< gpar[0][i] <<" "<< gpar[1][j] << std::endl;
#endif

  // Evaluate the secondary solution at all Gauss points
  Matrix sField;
  if (!this->evalSolution(sField,integrand,gpar.data()))
    return nullptr;

  // Compute parameter values of the Greville points
  if (!this->getGrevilleParameters(gpar[0],0)) return nullptr;
  if (!this->getGrevilleParameters(gpar[1],1)) return nullptr;

  const int n1 = p1 - m + 1; // Patch size in first parameter direction
  const int n2 = p2 - m + 1; // Patch size in second parameter direction

  const size_t nCmp = sField.rows(); // Number of result components
  const size_t nPol = (n1+1)*(n2+1); // Number of terms in polynomial expansion

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

#if SP_DEBUG > 1
      std::cout <<"\nGreville point "<< ig <<","<< jg <<" (u,v) = "
		<< gpar[0][ig] <<" "<< gpar[1][jg] <<" X = "<< G << std::endl;
#endif

      // Set up the local projection matrices
      DenseMatrix A(nPol,nPol,true);
      Matrix B(nPol,nCmp);

      // Loop over all non-zero knot-spans in the support of
      // the basis function associated with current Greville point
      for (int js = jstart; js < jstart+n2; js++)
	if (js >= p2-1 && surf->knotSpan(1,js) > 0.0)
	  for (int is = istart; is < istart+n1; is++)
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
		  evalMonomials(n1+1,n2+1,X[0]-G[0],X[1]-G[1],P);

#if SP_DEBUG > 1
		  std::cout <<"Itg. point "<< i <<","<< j <<" (u,v) = "
			    << u <<" "<< v <<" X = "<< X <<" P-matrix:"<< P;
#endif

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

#if SP_DEBUG > 2
      std::cout <<"---- Matrix A -----"<< A
                <<"-------------------"<< std::endl;
      std::cout <<"---- Vector B -----"<< B
                <<"-------------------"<< std::endl;
#endif

      // Solve the local equation system
      if (!A.solve(B)) return nullptr;

      // Evaluate the projected field at current Greville point (first row of B)
      for (l = 1; l <= nCmp; l++)
	sValues(l,ip) = B(1,l);

#if SP_DEBUG > 1
      std::cout <<"Greville point "<< ig <<","<< jg <<" :";
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


/*!
  \note A Variation Diminishing Spline Approximation is used here as the
  regular interpolation method in GoTools only works with uniform knots.
*/

bool ASMs2D::evaluate (const Field* field, RealArray& vec, int basisNum) const
{
  // Compute parameter values of the result sampling points (Greville points)
  std::array<RealArray,2> gpar;
  for (int dir = 0; dir < 2; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir,basisNum))
      return false;

  // Evaluate the result field at all sampling points.
  // Note: it is here assumed that *basis and *this have spline bases
  // defined over the same parameter domain.
  Vector sValues;
  sValues.reserve(gpar[0].size()*gpar[1].size());
  for (double v : gpar[1])
    for (double u : gpar[0])
      sValues.push_back(field->valueFE(ItgPoint(u,v)));

  // Project the results onto the spline basis to find control point
  // values based on the result values evaluated at the Greville points.
  // Note that we here implicitly assume that the number of Greville points
  // equals the number of control points such that we don't have to resize
  // the result array. Think that is always the case, but beware if trying
  // other projection schemes later.

  Go::SplineSurface* surf_new =
    VariationDiminishingSplineApproximation(this->getBasis(basisNum),sValues,1);

  vec.assign(surf_new->coefs_begin(),surf_new->coefs_end());
  delete surf_new;

  return true;
}


// L2-Projection: Least-square approximation; global approximation
Go::SplineSurface* ASMs2D::projectSolutionLeastSquare (const IntegrandBase& integrand) const
{
  if (!surf) return nullptr;

  // Get Gaussian quadrature points and weights
  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int ng = this->getNoGaussPt(p1 > p2 ? p1 : p2, true);

  const double* xg = GaussQuadrature::getCoord(ng);
  const double* wg = GaussQuadrature::getWeight(ng);
  if (!xg || !wg) return nullptr;

  // Compute parameter values of the result sampling points (Gauss points)
  std::array<Matrix,2> ggpar;
  std::array<RealArray,2> gpar, wgpar;
  for (int dir = 0; dir < 2; dir++)
  {
    this->getGaussPointParameters(ggpar[dir],dir,ng,xg);
    gpar[dir] = ggpar[dir];

    // Gauss weights at parameter values
    const Go::BsplineBasis& basis = surf->basis(dir);
    RealArray::const_iterator knotit = basis.begin();
    RealArray& tmp = wgpar[dir];
    tmp.reserve(ng*(basis.numCoefs()-basis.order()));
    for (int i = basis.order(); i <= basis.numCoefs(); i++)
    {
      double d = knotit[i] - knotit[i-1];
      for (int j = 0; j < ng; j++)
        tmp.push_back(d > 0.0 ? wg[j]*d*0.5 : 0.0);
    }
  }

  // Evaluate the secondary solution at all sampling points (Gauss points)
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar.data()))
    return nullptr;

  RealArray weights;
  if (surf->rational())
    surf->getWeights(weights);

  return leastsquare_approximation(surf->basis(0),
				   surf->basis(1),
				   gpar[0], gpar[1],
				   wgpar[0], wgpar[1],
				   sValues,
				   sValues.rows(),
				   surf->rational(),
				   weights);
}


// Quasi-Interpolation; local interpolation method
Go::SplineSurface* ASMs2D::projectSolutionLocal (const IntegrandBase& integrand) const
{
  // Compute parameter values of the result sampling points
  // (the Quasi-Interpolation points)
  std::array<RealArray,2> gpar;
  for (int dir = 0; dir < 2; dir++)
    if (!this->getQuasiInterplParameters(gpar[dir],dir))
      return nullptr;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar.data()))
    return nullptr;

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
Go::SplineSurface* ASMs2D::projectSolutionLocalApprox (const IntegrandBase& integrand) const
{
  // Compute parameter values of the result sampling points (Greville points)
  std::array<RealArray,2> gpar;
  for (int dir = 0; dir < 2; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir))
      return nullptr;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar.data()))
    return nullptr;

  // Project onto the geometry basis
  return VariationDiminishingSplineApproximation(surf.get(),sValues,sValues.rows());
}
