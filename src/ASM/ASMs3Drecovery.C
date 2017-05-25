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
#include "FiniteElement.h"
#include "Field.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "SparseMatrix.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Profiler.h"
#include <array>


bool ASMs3D::getGrevilleParameters (RealArray& prm, int dir, int basisNum) const
{
  if (dir < 0 || dir > 2) return false;

  const Go::BsplineBasis& basis = this->getBasis(basisNum)->basis(dir);

  prm.resize(basis.numCoefs());
  for (size_t i = 0; i < prm.size(); i++)
    prm[i] = basis.grevilleParameter(i);

  return true;
}


bool ASMs3D::getQuasiInterplParameters (RealArray& prm, int dir) const
{
  if (!svol || dir < 0 || dir > 2) return false;

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


bool ASMs3D::evaluate (const ASMbase* basis, const Vector& locVec,
		       RealArray& vec, int basisNum) const
{
  const ASMs3D* pch = dynamic_cast<const ASMs3D*>(basis);
  if (!pch) return false;

  // Compute parameter values of the result sampling points (Greville points)
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir,basisNum))
      return false;

  // Evaluate the result field at all sampling points.
  // Note: it is here assumed that *basis and *this have spline bases
  // defined over the same parameter domain.
  Matrix sValues;
  if (!pch->evalSolution(sValues,locVec,gpar.data()))
    return false;

  Go::SplineVolume* svol = this->getBasis(basisNum);

  // Project the results onto the spline basis to find control point
  // values based on the result values evaluated at the Greville points.
  // Note that we here implicitly assume that the number of Greville points
  // equals the number of control points such that we don't have to resize
  // the result array. Think that is always the case, but beware if trying
  // other projection schemes later.

  RealArray weights;
  if (svol->rational())
    svol->getWeights(weights);

  const Vector& vec2 = sValues;
  Go::SplineVolume* vol_new =
    Go::VolumeInterpolator::regularInterpolation(svol->basis(0),
						 svol->basis(1),
						 svol->basis(2),
						 gpar[0], gpar[1], gpar[2],
						 const_cast<Vector&>(vec2),
						 sValues.rows(),
						 svol->rational(),
						 weights);

  vec.assign(vol_new->coefs_begin(),vol_new->coefs_end());
  delete vol_new;

  return true;
}


Go::SplineVolume* ASMs3D::projectSolution (const IntegrandBase& integrand) const
{
  PROFILE2("ASMs3D::projectSolution");

  const int basis = 1;

  // Compute parameter values of the result sampling points (Greville points)
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir,basis))
      return nullptr;

  const Go::SplineVolume* pvol = this->getBasis(basis);

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar.data()) || sValues.rows() == 0)
    return nullptr;

  // Project the results onto the spline basis to find control point
  // values based on the result values evaluated at the Greville points.
  // Note that we here implicitly assume that the number of Greville points
  // equals the number of control points such that we don't have to resize
  // the result array. Think that is always the case, but beware if trying
  // other projection schemes later.

  RealArray weights;
  if (pvol->rational())
    pvol->getWeights(weights);

  const Vector& vec = sValues;
  return Go::VolumeInterpolator::regularInterpolation(pvol->basis(0),
                                                      pvol->basis(1),
                                                      pvol->basis(2),
                                                      gpar[0], gpar[1], gpar[2],
                                                      const_cast<Vector&>(vec),
                                                      sValues.rows(),
                                                      pvol->rational(),
                                                      weights);
}


Go::GeomObject* ASMs3D::evalSolution (const IntegrandBase& integrand) const
{
  return this->projectSolution(integrand);
}


bool ASMs3D::assembleL2matrices (SparseMatrix& A, StdVector& B,
                                 const IntegrandBase& integrand,
                                 bool continuous) const
{
  const size_t nnod = this->getNoNodes(1);

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
  std::array<RealArray,3> gpar;
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
  Matrix sField;
  if (!this->evalSolution(sField,integrand,gpar.data()))
  {
    std::cerr <<" *** ASMs3D::assembleL2matrices: Failed for patch "<< idx+1
	      <<" nPoints="<< gpar[0].size()*gpar[1].size()*gpar[2].size()
	      << std::endl;
    return false;
  }

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
		SplineUtils::extractBasis(spl1[ip],phi,dNdu);
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
		for (size_t r = 1; r <= sField.rows(); r++)
		  B(inod+(r-1)*nnod) += phi[ii]*sField(r,ip+1)*dJw;
	      }
	    }
      }

  return true;
}


#include "ASMs3DInterpolate.C" // TODO: inline these methods instead...


/*!
  \note A Variation Diminishing Spline Approximation is used here as the
  regular interpolation method in GoTools only works with uniform knots.
*/

bool ASMs3D::evaluate (const Field* field, RealArray& vec, int basisNum) const
{
  // Compute parameter values of the result sampling points (Greville points)
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir,basisNum))
      return false;

  // Evaluate the result field at all sampling points.
  // Note: it is here assumed that *basis and *this have spline bases
  // defined over the same parameter domain.
  Vector sValues(gpar[0].size()*gpar[1].size()*gpar[2].size());
  Vector::iterator it=sValues.begin();
  for (size_t l=0;l<gpar[2].size();++l) {
    FiniteElement fe;
    fe.w = gpar[2][l];
    for (size_t j=0;j<gpar[1].size();++j) {
      fe.v = gpar[1][j];
      for (size_t i=0;i<gpar[0].size();++i) {
        fe.u = gpar[0][i];
        *it++ = field->valueFE(fe);
      }
    }
  }

  Go::SplineVolume* svol = this->getBasis(basisNum);

  // Project the results onto the spline basis to find control point
  // values based on the result values evaluated at the Greville points.
  // Note that we here implicitly assume that the number of Greville points
  // equals the number of control points such that we don't have to resize
  // the result array. Think that is always the case, but beware if trying
  // other projection schemes later.

  RealArray weights;
  if (svol->rational())
    svol->getWeights(weights);

  Go::SplineVolume* vol_new =
    VariationDiminishingSplineApproximation(svol->basis(0),
                                            svol->basis(1),
                                            svol->basis(2),
                                            gpar[0], gpar[1], gpar[2],
                                            sValues, 1, svol->rational(),
                                            weights);

  vec.assign(vol_new->coefs_begin(),vol_new->coefs_end());
  delete vol_new;

  return true;
}


Go::SplineVolume* ASMs3D::projectSolutionLeastSquare (const IntegrandBase& integrand) const
{
  if (!svol) return nullptr;

  PROFILE1("test L2- projection");
  // Compute parameter values of the result sampling points (Gauss-Interpl. points)
  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return nullptr;

  std::array<Matrix,3> ggpar;
  std::array<RealArray,3> gpar, wgpar;
  for (int dir = 0; dir < 3; dir++)
  {
    this->getGaussPointParameters(ggpar[dir],dir,nGauss,xg);
    gpar[dir] = ggpar[dir];

    // Gauss weights at parameter values
    const Go::BsplineBasis& basis = svol->basis(dir);
    RealArray::const_iterator knotit = basis.begin();
    RealArray& tmp = wgpar[dir];
    tmp.reserve(nGauss*(basis.numCoefs()-basis.order()));
    for (int i = basis.order(); i <= basis.numCoefs(); i++)
    {
      double d = knotit[i]-knotit[i-1];
      for (int j = 0; j < nGauss; j++)
        tmp.push_back(d > 0.0 ? wg[j]*d*0.5 : 0.0);
    }
  }

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar.data()))
    return nullptr;

  RealArray weights;
  if (svol->rational())
    svol->getWeights(weights);

  return leastsquare_approximation(svol->basis(0),
                                   svol->basis(1),
                                   svol->basis(2),
                                   gpar[0], gpar[1], gpar[2],
                                   wgpar[0], wgpar[1], wgpar[2],
                                   sValues,
                                   sValues.rows(),
                                   svol->rational(),
                                   weights);
}


Go::SplineVolume* ASMs3D::projectSolutionLocal (const IntegrandBase& integrand) const
{
  PROFILE1("test Quasi projection");
  // Compute parameter values of the result sampling points (Greville points)
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 2; dir++)
    if (!this->getQuasiInterplParameters(gpar[dir],dir))
      return nullptr;

  for (int dir = 2; dir < 3; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir))
      return nullptr;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar.data()))
    return nullptr;

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
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir))
      return nullptr;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar.data()))
    return nullptr;

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
