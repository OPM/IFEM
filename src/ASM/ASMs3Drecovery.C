// $Id$
//==============================================================================
//!
//! \file ASMs3Drecovery.C
//!
//! \date Mar 06 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Recovery of secondary solutions for structured 3D spline FE models.
//!
//==============================================================================

#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/VolumeInterpolator.h"

#include "ASMs3D.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "SparseMatrix.h"
#include "Utilities.h"
#include "Profiler.h"


bool ASMs3D::getGrevilleParameters (RealArray& prm, int dir) const
{
  if (!svol || dir < 0 || dir > 2) return false;

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


Go::SplineVolume* ASMs3D::projectSolution (const IntegrandBase& integrand) const
{
  PROFILE1("ASMs3D::projectSolution");

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


Go::GeomObject* ASMs3D::evalSolution (const IntegrandBase& integrand) const
{
  return this->projectSolution(integrand);
}


bool ASMs3D::globalL2projection (Matrix& sField,
				 const IntegrandBase& integrand,
				 bool continuous) const
{
  if (!svol) return true; // silently ignore empty patches

  PROFILE2("ASMs3D::globalL2");

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);
  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);
  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;
  const int nel3 = n3 - p3 + 1;

  // Get Gaussian quadrature point coordinates (and weights if continuous)
  const int ng1 = continuous ? nGauss : p1 - 1;
  const int ng2 = continuous ? nGauss : p2 - 1;
  const int ng3 = continuous ? nGauss : p3 - 1;
  const double* xg = GaussQuadrature::getCoord(ng1);
  const double* yg = GaussQuadrature::getCoord(ng2);
  const double* zg = GaussQuadrature::getCoord(ng3);
  const double* wg = continuous ? GaussQuadrature::getWeight(nGauss) : 0;
  if (!xg || !yg || !zg) return false;
  if (continuous && !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  Matrix gp;
  RealArray gpar[3];
  gpar[0] = this->getGaussPointParameters(gp,0,ng1,xg);
  gpar[1] = this->getGaussPointParameters(gp,1,ng2,yg);
  gpar[2] = this->getGaussPointParameters(gp,2,ng3,zg);

  // Evaluate basis functions at all integration points
  std::vector<Go::BasisPts>    spl0;
  std::vector<Go::BasisDerivs> spl1;
  if (continuous)
    svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],spl1);
  else
    svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],spl0);

  // Evaluate the secondary solution at all integration points
  if (!this->evalSolution(sField,integrand,gpar))
    return false;

  // Set up the projection matrices
  const size_t nnod = this->getNoNodes(1);
  const size_t ncomp = sField.rows();
  SparseMatrix A(SparseMatrix::SUPERLU);
  StdVector B(nnod*ncomp);
  A.redim(nnod,nnod);

  double dV = 1.0;
  Vector phi(p1*p2*p3);
  Matrix dNdu, Xnod, J;


  // === Assembly loop over all elements in the patch ==========================

  int iel = 0;
  for (int i3 = 0; i3 < nel3; i3++)
    for (int i2 = 0; i2 < nel2; i2++)
      for (int i1 = 0; i1 < nel1; i1++, iel++)
      {
	if (MLGE[iel] < 1) continue; // zero-volume element

	if (continuous)
	{
	  // Set up control point (nodal) coordinates for current element
	  if (!this->getElementCoordinates(Xnod,1+iel))
	    return false;
	  else if ((dV = 0.125*this->getParametricVolume(1+iel)) < 0.0)
	    return false; // topology error (probably logic error)
	}

	// --- Integration loop over all Gauss points in each direction --------

        int ip = ((i3*ng2*nel2 + i2)*ng1*nel1 + i1)*ng3;
        for (int k = 0; k < ng3; k++, ip += ng2*(nel2-1)*ng1*nel1)
          for (int j = 0; j < ng2; j++, ip += ng1*(nel1-1))
            for (int i = 0; i < ng1; i++, ip++)
	    {
	      if (continuous)
		extractBasis(spl1[ip],phi,dNdu);
	      else
		phi = spl0[ip].basisValues;

	      // Compute the Jacobian inverse and derivatives
	      double dJw = dV;
	      if (continuous)
	      {
		dJw *= wg[i]*wg[j]*wg[k]*utl::Jacobian(J,dNdu,Xnod,dNdu,false);
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

  // Solve the patch-global equation system
  if (!A.solve(B)) return false;

  // Store the control-point values of the projected field
  sField.resize(ncomp,nnod);
  for (size_t i = 1; i <= nnod; i++)
    for (size_t j = 1; j <= ncomp; j++)
      sField(j,i) = B(i+(j-1)*nnod);

  return true;
}

#include "ASMs3DInterpolate.C" // TODO: inline these methods instead...

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

  return leastsquare_approximation(svol->basis(0),
                                   svol->basis(1),
                                   svol->basis(2),
                                   par_u, par_v, par_w,
                                   wgpar_u, wgpar_v, wgpar_w,
                                   sValues,
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
  
  return quasiInterpolation(svol->basis(0),
			    svol->basis(1),
			    svol->basis(2),
			    gpar[0], gpar[1], gpar[2],
			    sValues,
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

  return VariationDiminishingSplineApproximation(svol->basis(0),
						 svol->basis(1),
						 svol->basis(2),
						 gpar[0], gpar[1], gpar[2],
						 sValues,
						 sValues.rows(),
						 svol->rational(),
						 weights);
}
