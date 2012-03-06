// $Id$
//==============================================================================
//!
//! \file ASMs3Drecovery.C
//!
//! \date March 06 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Recovery of secondary solutions for structured 3D spline FE models.
//!
//==============================================================================

#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/VolumeInterpolator.h"

#include "ASMs3D.h"
#include "GaussQuadrature.h"
#include "Profiler.h"


bool ASMs3D::getGrevilleParameters (RealArray& prm, int dir) const
{
  if (!svol) return false;

  const Go::BsplineBasis& basis = svol->basis(dir);

  prm.resize(basis.numCoefs());
  for (size_t i = 0; i < prm.size(); i++)
    prm[i] = basis.grevilleParameter(i);

  return true;
}


bool ASMs3D::getQuasiInterplParameters (RealArray& prm, int dir) const
{
  if (!svol) return false;
  const Go::BsplineBasis& basis = svol->basis(dir);
  
  std::vector< double > knots_simple;
  basis.knotsSimple(knots_simple);
  
  prm.clear();
  for (size_t i = 0; i < knots_simple.size(); i++){
    prm.push_back(knots_simple[i]);
    prm.push_back(0.5*(knots_simple[i]+knots_simple[i+1]));}
  prm.pop_back();
  
  return true;
}


Go::GeomObject* ASMs3D::evalSolution (const IntegrandBase& integrand) const
{
  return this->projectSolution(integrand);
}


Go::SplineVolume* ASMs3D::projectSolution (const IntegrandBase& integrand) const
{
  PROFILE1("test Global projection");
  // Compute parameter values of the result sampling points (Greville points)
  RealArray gpar[3];
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir))
      return 0;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar) || sValues.rows() == 0)
    return 0;

  // Project the results onto the spline basis to find control point
  // values based on the result values evaluated at the Greville points.
  // Note that we here implicitly assume that the number of Greville points
  // equals the number of control points such that we don't have to resize
  // the result array. Think that is always the case, but beware if trying
  // other projection schemes later.

  RealArray weights;
  if (svol->rational())
    svol->getWeights(weights);

  const Vector& vec = sValues;
  return Go::VolumeInterpolator::regularInterpolation(svol->basis(0),
						      svol->basis(1),
						      svol->basis(2),
						      gpar[0], gpar[1], gpar[2],
						      const_cast<Vector&>(vec),
						      sValues.rows(),
						      svol->rational(),
						      weights);
}


Go::SplineVolume* ASMs3D::projectSolutionLeastSquare (const IntegrandBase& integrand) const
{
  PROFILE1("test L2- projection");
  // Compute parameter values of the result sampling points (Gauss-Interpl. points)
  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  Matrix ggpar[3];
  for (int dir = 0; dir < 3; dir++)
    this->getGaussPointParameters(ggpar[dir],dir,nGauss,xg);

  std::vector<double> par_u;
  par_u = ggpar[0];
  std::vector<double> par_v;
  par_v = ggpar[1];
  std::vector<double> par_w;
  par_w = ggpar[2];

  // gauss weights at parameter values
  std::vector<double> wgpar_u;
  std::vector<double> wgpar_v;
  std::vector<double> wgpar_w;
  double d;
  for (int dir = 0; dir < 3; dir++){
    if (!svol) return false;
    const Go::BsplineBasis& basis = svol->basis(dir);
    RealArray::const_iterator knotit = svol->basis(dir).begin();
    std::vector<double> tmp;
    tmp.reserve(nGauss*(basis.numCoefs()-basis.order()));
    for (size_t i = 0; i<=(basis.numCoefs()-basis.order());i++)
    {
      d = knotit[i+basis.order()]-knotit[i+basis.order()-1];
      if (d > 0)
      {for (int j = 0;j<nGauss;j++)
        tmp.push_back(wg[j]/(2/d));
      }
      else if (d == 0)
      {for (int j = 0;j<nGauss;j++)
        tmp.push_back(0);}
    }  
    if (dir == 0) 
      wgpar_u = tmp; 
    else if (dir == 1)
      wgpar_v = tmp;
    else if (dir == 2)
      wgpar_w = tmp;
  }

  RealArray gpar[3];
  gpar[0] = par_u;
  gpar[1] = par_v;
  gpar[2] = par_w;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar))
    return 0;

  RealArray weights;
  if (svol->rational())
    svol->getWeights(weights);

  const Vector& vec = sValues;
  return leastsquare_approximation(svol->basis(0),
                                   svol->basis(1),
                                   svol->basis(2),
                                   par_u, par_v, par_w,
                                   wgpar_u, wgpar_v, wgpar_w,
                                   const_cast<Vector&>(vec),
                                   sValues.rows(),
                                   svol->rational(),
                                   weights);
}


Go::SplineVolume* ASMs3D::projectSolutionLocal (const IntegrandBase& integrand) const
{
  PROFILE1("test Quasi projection");
  // Compute parameter values of the result sampling points (Greville points)
  RealArray gpar[3];
  for (int dir = 0; dir < 2; dir++)
    if (!this->getQuasiInterplParameters(gpar[dir],dir))
      return 0;
  
  for (int dir = 2; dir < 3; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir))
      return 0;
  
  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar))
    return 0;
  
  RealArray weights;
  if (svol->rational())
    svol->getWeights(weights);
  
  const Vector& vec = sValues;
  return quasiInterpolation(svol->basis(0),
			    svol->basis(1),
			    svol->basis(2),
			    gpar[0], gpar[1], gpar[2],
			    const_cast<Vector&>(vec),
			    sValues.rows(),
			    svol->rational(),
			    weights);
}


Go::SplineVolume* ASMs3D::projectSolutionLocalApprox(const IntegrandBase& integrand) const
{
  PROFILE1("test VDSA projection");      
  // Compute parameter values of the result sampling points (Greville points)
  RealArray gpar[3];
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir))
      return 0;
  
  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar))
    return 0;

  RealArray weights;
  if (svol->rational())
    svol->getWeights(weights);
  
  const Vector& vec = sValues;
  return VariationDiminishingSplineApproximation(svol->basis(0),
						 svol->basis(1),
						 svol->basis(2),
						 gpar[0], gpar[1], gpar[2],
						 const_cast<Vector&>(vec),
						 sValues.rows(),
						 svol->rational(),
						 weights);
}
