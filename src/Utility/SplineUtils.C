// $Id$
//==============================================================================
//!
//! \file SplineUtils.C
//!
//! \date Mar 29 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Various utility functions on spline objects.
//!
//==============================================================================

#include "SplineUtils.h"
#include "MatVec.h"
#include "Vec3.h"

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/CurveInterpolator.h"
#include "GoTools/geometry/SurfaceInterpolator.h"


Vec3 SplineUtils::toVec3 (const Go::Point& X, int nsd)
{
  Vec3 Y;
  for (int i = 0; i < nsd && i < X.size() && i < 3; i++) Y[i] = X[i];
  return Y;
}


Vec4 SplineUtils::toVec4 (const Go::Point& X, Real time)
{
  Vec4 Y;
  for (int i = 0; i < X.size() && i < 3; i++) Y[i] = X[i];
  Y.t = time;
  return Y;
}


void SplineUtils::extractBasis (const Go::BasisDerivsSf& spline,
                                Vector& N, Matrix& dNdu)
{
  dNdu.resize(N.size(),2);

  size_t jp, n = 1;
  for (jp = 0; jp < N.size(); jp++, n++)
  {
     N  (n)   = spline.basisValues[jp];
    dNdu(n,1) = spline.basisDerivs_u[jp];
    dNdu(n,2) = spline.basisDerivs_v[jp];
  }
}


void SplineUtils::extractBasis (const Go::BasisDerivsSf2& spline,
                                Vector& N, Matrix& dNdu, Matrix3D& d2Ndu2)
{
   dNdu .resize(N.size(),2);
  d2Ndu2.resize(N.size(),2,2);

  size_t jp, n = 1;
  for (jp = 0; jp < N.size(); jp++, n++)
  {
      N   (n)     = spline.basisValues[jp];
     dNdu (n,1)   = spline.basisDerivs_u[jp];
     dNdu (n,2)   = spline.basisDerivs_v[jp];
    d2Ndu2(n,1,1) = spline.basisDerivs_uu[jp];
    d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = spline.basisDerivs_uv[jp];
    d2Ndu2(n,2,2) = spline.basisDerivs_vv[jp];
  }
}


void SplineUtils::extractBasis (const Go::BasisDerivs& spline,
                                Vector& N, Matrix& dNdu)
{
  dNdu.resize(N.size(),3);

  size_t jp, n = 1;
  for (jp = 0; jp < N.size(); jp++, n++)
  {
     N  (n)   = spline.basisValues[jp];
    dNdu(n,1) = spline.basisDerivs_u[jp];
    dNdu(n,2) = spline.basisDerivs_v[jp];
    dNdu(n,3) = spline.basisDerivs_w[jp];
  }
}


void SplineUtils::extractBasis (const Go::BasisDerivs2& spline,
                                Vector& N, Matrix& dNdu, Matrix3D& d2Ndu2)
{
   dNdu .resize(N.size(),3);
  d2Ndu2.resize(N.size(),3,3);

  size_t jp, n = 1;
  for (jp = 0; jp < N.size(); jp++, n++)
  {
      N   (n)     = spline.basisValues[jp];
     dNdu (n,1)   = spline.basisDerivs_u[jp];
     dNdu (n,2)   = spline.basisDerivs_v[jp];
     dNdu (n,3)   = spline.basisDerivs_w[jp];
    d2Ndu2(n,1,1) = spline.basisDerivs_uu[jp];
    d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = spline.basisDerivs_uv[jp];
    d2Ndu2(n,1,3) = d2Ndu2(n,3,1) = spline.basisDerivs_uw[jp];
    d2Ndu2(n,2,2) = spline.basisDerivs_vv[jp];
    d2Ndu2(n,2,3) = d2Ndu2(n,3,2) = spline.basisDerivs_vw[jp];
    d2Ndu2(n,3,3) = spline.basisDerivs_ww[jp];
  }
}


Go::SplineCurve* SplineUtils::project (const Go::SplineCurve* curve,
				       const RealFunc& f, Real time)
{
  if (!curve) return NULL;

  const Go::BsplineBasis& basis = curve->basis();
  const int nPoints = basis.numCoefs();

  RealArray gpar(nPoints), fval(nPoints);
  Go::Point X;

  // Compute parameter values of the function sampling points (Greville points)
  // and evaluate the function at these points
  for (int i = 0; i < nPoints; i++)
  {
    gpar[i] = basis.grevilleParameter(i);
    curve->point(X,gpar[i]);
    fval[i] = f(toVec4(X,time));
  }

  // Get weights for rational spline curves (NURBS)
  RealArray weights;
  if (curve->rational())
    curve->getWeights(weights);

  // Project the function onto the spline basis to find control point values
  return Go::CurveInterpolator::regularInterpolation(basis,gpar,fval,1,
						     curve->rational(),weights);
}


Go::SplineCurve* SplineUtils::project (const Go::SplineCurve* curve,
				       const VecFunc& f, int nComp, Real time)
{
  if (!curve || nComp < 1) return NULL;
  if (nComp > 3) nComp = 3;

  const Go::BsplineBasis& basis = curve->basis();
  const int nPoints = basis.numCoefs();

  RealArray gpar(nPoints), fval(nComp*nPoints);
  Go::Point X;
  Vec3 fOfX;

  // Compute parameter values of the function sampling points (Greville points)
  // and evaluate the function at these points
  int i, j, k;
  for (i = k = 0; i < nPoints; i++)
  {
    gpar[i] = basis.grevilleParameter(i);
    curve->point(X,gpar[i]);
    fOfX = f(toVec4(X,time));
    for (j = 0; j < nComp; j++, k++)
      fval[k] = fOfX[j];
  }

  // Get weights for rational spline curves (NURBS)
  RealArray weights;
  if (curve->rational())
    curve->getWeights(weights);

  // Project the function onto the spline basis to find control point values
  return Go::CurveInterpolator::regularInterpolation(basis,gpar,fval,nComp,
						     curve->rational(),weights);
}


Go::SplineSurface* SplineUtils::project (const Go::SplineSurface* surface,
					 const RealFunc& f, Real time)
{
  if (!surface) return NULL;

  const Go::BsplineBasis& ubas = surface->basis(0);
  const Go::BsplineBasis& vbas = surface->basis(1);
  const int nu = ubas.numCoefs();
  const int nv = vbas.numCoefs();

  RealArray upar(nu), vpar(nv);

  // Compute parameter values of the function sampling points (Greville points)
  int i, j;
  for (i = 0; i < nu; i++)
    upar[i] = ubas.grevilleParameter(i);
  for (j = 0; j < nv; j++)
    vpar[j] = vbas.grevilleParameter(j);

  // Evaluate the function at the sampling points
  Go::Point X;
  RealArray fval(nu*nv);
  size_t k = 0;
  for (j = 0; j < nv; j++)
    for (i = 0; i < nu; i++, k++)
    {
      surface->point(X,upar[i],vpar[j]);
      fval[k] = f(toVec4(X,time));
    }

  // Get weights for rational spline curves (NURBS)
  RealArray weights;
  if (surface->rational())
    surface->getWeights(weights);

  // Project the function onto the spline basis to find control point values
  return Go::SurfaceInterpolator::regularInterpolation(ubas,vbas,
						       upar,vpar,fval,1,
						       surface->rational(),
						       weights);
}


Go::SplineSurface* SplineUtils::project (const Go::SplineSurface* surface,
					 const VecFunc& f, int nComp, Real time)
{
  if (!surface || nComp < 1) return NULL;
  if (nComp > 3) nComp = 3;

  const Go::BsplineBasis& ubas = surface->basis(0);
  const Go::BsplineBasis& vbas = surface->basis(1);
  const int nu = ubas.numCoefs();
  const int nv = vbas.numCoefs();

  RealArray upar(nu), vpar(nv);

  // Compute parameter values of the function sampling points (Greville points)
  int i, j;
  for (i = 0; i < nu; i++)
    upar[i] = ubas.grevilleParameter(i);
  for (j = 0; j < nv; j++)
    vpar[j] = vbas.grevilleParameter(j);

  // Evaluate the function at the sampling points
  Go::Point X;
  Vec3 fOfX;
  RealArray fval(nComp*nu*nv);
  size_t k = 0;
  for (j = 0; j < nv; j++)
    for (i = 0; i < nu; i++)
    {
      surface->point(X,upar[i],vpar[j]);
      fOfX = f(toVec4(X,time));
      for (int l = 0; l < nComp; l++, k++)
	fval[k] = fOfX[l];
    }

  // Get weights for rational spline curves (NURBS)
  RealArray weights;
  if (surface->rational())
    surface->getWeights(weights);

  // Project the function onto the spline basis to find control point values
  return Go::SurfaceInterpolator::regularInterpolation(ubas,vbas,
						       upar,vpar,fval,nComp,
						       surface->rational(),
						       weights);
}
