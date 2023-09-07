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
#include "Function.h"
#include "Vec3.h"

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveInterpolator.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SurfaceInterpolator.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/VolumeInterpolator.h"


Vec3 SplineUtils::toVec3 (const Go::Point& X, int nsd)
{
  Vec3 Y;
  for (int i = 0; i < nsd && i < X.size() && i < 3; i++) Y[i] = X[i];
  return Y;
}


Vec4 SplineUtils::toVec4 (const Go::Point& X, Real time,
                          const double* u)
{
  Vec4 Y;
  for (int i = 0; i < X.size() && i < 3; i++) Y[i] = X[i];
  Y.t = time;
  Y.u = u;
  return Y;
}


void SplineUtils::point (Vec3& X, double u, const Go::SplineCurve* curve)
{
  Go::Point Y;
#pragma omp critical
  curve->point(Y,u);
  for (int i = 0; i < Y.size() && i < 3; i++) X[i] = Y[i];  
}


void SplineUtils::point (Vec3& X, double u, double v, const Go::SplineSurface* surf)
{
  Go::Point Y;
#pragma omp critical
  surf->point(Y,u,v);
  for (int i = 0; i < Y.size() && i < 3; i++) X[i] = Y[i];  
}


void SplineUtils::point (Vec3& X, double u, double v, double w,
                         const Go::SplineVolume* vol)
{
  Go::Point Y;
#pragma omp critical
  vol->point(Y,u,v,w);
  for (int i = 0; i < Y.size() && i < 3; i++) X[i] = Y[i];  
}


void SplineUtils::extractBasis (const Go::BasisDerivsSf& spline,
                                Vector& N, Matrix& dNdu)
{
   N  .resize(spline.basisValues.size());
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
    N   .resize(spline.basisValues.size());
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


void SplineUtils::extractBasis (const Go::BasisDerivsSf3& spline,
                                Vector& N, Matrix& dNdu,
                                Matrix3D& d2Ndu2, Matrix4D& d3Ndu3)
{
    N   .resize(spline.basisValues.size());
   dNdu .resize(N.size(),2);
  d2Ndu2.resize(N.size(),2,2);
  d3Ndu3.resize(N.size(),2,2,2);

  size_t jp, n = 1;
  for (jp = 0; jp < N.size(); jp++, n++)
  {
      N   (n)     = spline.basisValues[jp];
     dNdu (n,1)   = spline.basisDerivs_u[jp];
     dNdu (n,2)   = spline.basisDerivs_v[jp];
    d2Ndu2(n,1,1) = spline.basisDerivs_uu[jp];
    d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = spline.basisDerivs_uv[jp];
    d2Ndu2(n,2,2) = spline.basisDerivs_vv[jp];
    d3Ndu3(n,1,1,1) = spline.basisDerivs_uuu[jp];
    d3Ndu3(n,2,1,1) = spline.basisDerivs_uuv[jp];
    d3Ndu3(n,1,2,1) = spline.basisDerivs_uuv[jp];
    d3Ndu3(n,1,1,2) = spline.basisDerivs_uuv[jp];
    d3Ndu3(n,2,1,1) = spline.basisDerivs_uuv[jp];
    d3Ndu3(n,2,2,1) = spline.basisDerivs_uvv[jp];
    d3Ndu3(n,1,2,2) = spline.basisDerivs_uvv[jp];
    d3Ndu3(n,2,1,2) = spline.basisDerivs_uvv[jp];
    d3Ndu3(n,2,2,2) = spline.basisDerivs_vvv[jp];
  }
}


void SplineUtils::extractBasis (const Go::BasisDerivs& spline,
                                Vector& N, Matrix& dNdu)
{
   N  .resize(spline.basisValues.size());
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
    N   .resize(spline.basisValues.size());
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
                                       const FunctionBase& f,
                                       int nComp, Real time)
{
  if (!curve || nComp < 1) return nullptr;

  const Go::BsplineBasis& basis = curve->basis();
  const int nPoints = basis.numCoefs();

  Go::Point X;
  RealArray gpar(nPoints), fval, fOfX;
  fval.reserve(nComp*nPoints);

  // Compute parameter values of the function sampling points (Greville points)
  // and evaluate the function at these points
  for (int i = 0; i < nPoints; i++)
  {
    gpar[i] = basis.grevilleParameter(i);
    curve->point(X,gpar[i]);
    fOfX = f.getValue(toVec4(X,time,&gpar[i]));
    fval.insert(fval.end(),fOfX.begin(),fOfX.begin()+nComp);
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
                                         const FunctionBase& f,
                                         int nComp, Real time)
{
  if (!surface || nComp < 1) return nullptr;

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
  RealArray fval, fOfX;
  fval.reserve(nComp*nu*nv);
  for (j = 0; j < nv; j++)
    for (i = 0; i < nu; i++)
    {
      surface->point(X,upar[i],vpar[j]);
      double u[2] = {upar[i], vpar[j]};
      fOfX = f.getValue(toVec4(X,time,u));
      fval.insert(fval.end(),fOfX.begin(),fOfX.begin()+nComp);
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


Go::SplineVolume* SplineUtils::project (const Go::SplineVolume* volume,
                                        const FunctionBase& f,
                                        int nComp, Real time)
{
  if (!volume || nComp < 1) return nullptr;

  const Go::BsplineBasis& ubas = volume->basis(0);
  const Go::BsplineBasis& vbas = volume->basis(1);
  const Go::BsplineBasis& wbas = volume->basis(2);
  const int nu = ubas.numCoefs();
  const int nv = vbas.numCoefs();
  const int nw = wbas.numCoefs();

  RealArray upar(nu), vpar(nv), wpar(nw);

  // Compute parameter values of the function sampling points (Greville points)
  int i, j, k;
  for (i = 0; i < nu; i++)
    upar[i] = ubas.grevilleParameter(i);
  for (j = 0; j < nv; j++)
    vpar[j] = vbas.grevilleParameter(j);
  for (k = 0; k < nw; k++)
    wpar[k] = wbas.grevilleParameter(k);

  // Evaluate the function at the sampling points
  Go::Point X;
  RealArray fval, fOfX;
  fval.reserve(nComp*nu*nv*nw);
  for (k = 0; k < nw; k++)
    for (j = 0; j < nv; j++)
      for (i = 0; i < nu; i++)
      {
        volume->point(X,upar[i],vpar[j],wpar[k]);
        double u[3] = {upar[i], vpar[j], wpar[k]};
        fOfX = f.getValue(toVec4(X,time,u));
        fval.insert(fval.end(),fOfX.begin(),fOfX.begin()+nComp);
      }

  // Get weights for rational spline curves (NURBS)
  RealArray weights;
  if (volume->rational())
    volume->getWeights(weights);

  // Project the function onto the spline basis to find control point values
  return Go::VolumeInterpolator::regularInterpolation(ubas,vbas,wbas,
                                                      upar,vpar,wpar,fval,nComp,
                                                      volume->rational(),
                                                      weights);
}


std::vector<double> SplineUtils::buildKnotVector(int p,
                                                 const std::vector<double>& knots,
                                                 const std::vector<int>& cont)
{
  std::vector<double> result;
  result.reserve(p+knots.size()+1);
  for (size_t i = 0; i < knots.size(); ++i)
    for (int k = 0; k < p-cont[i]; ++k)
      result.push_back(knots[i]);

  return result;
}
