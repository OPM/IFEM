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
#include "GoTools/geometry/CurveInterpolator.h"


Vec3 SplineUtils::toVec3 (const Go::Point& X, int nsd)
{
  Vec3 Y;
  for (int i = 0; i < nsd && i < X.size() && i < 3; i++) Y[i] = X[i];
  return Y;
}


Vec4 SplineUtils::toVec4 (const Go::Point& X, real time)
{
  Vec4 Y;
  for (int i = 0; i < X.size() && i < 3; i++) Y[i] = X[i];
  Y.t = time;
  return Y;
}


Go::SplineCurve* SplineUtils::project (const Go::SplineCurve* curve,
				       const RealFunc& f, real time)
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


Go::SplineCurve* SplineUtils::project(const Go::SplineCurve* curve,
				      const VecFunc& f, int nComp, real time)
{
  if (!curve || nComp < 1 || nComp > 3) return NULL;

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
