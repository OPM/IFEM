// $Id$
//==============================================================================
//!
//! \file SplineField2D.C
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for spline-based finite element scalar field in 2D.
//!
//==============================================================================

#include "GoTools/geometry/SplineSurface.h"

#include "SplineField2D.h"
#include "ASMs2D.h"
#include "FiniteElement.h"
#include "CoordinateMapping.h"
#include "Utilities.h"
#include "Vec3.h"
#include <array>


SplineField2D::SplineField2D (const ASMs2D* patch,
                              const RealArray& v, char nbasis,
                              const char* name)
  : FieldBase(name), basis(patch->getBasis(nbasis)), surf(patch->getSurface())
{
  const int n1 = basis->numCoefs_u();
  const int n2 = basis->numCoefs_v();
  nno = n1*n2;

  const int p1 = basis->order_u();
  const int p2 = basis->order_v();
  nelm = (n1-p1+1)*(n2-p2+1);

  // Ensure the values array has compatible length, pad with zeros if necessary
  values.resize(nno);
  RealArray::const_iterator end = v.size() > nno ? v.begin()+nno : v.end();
  std::copy(v.begin(),end,values.begin());
}


double SplineField2D::valueNode (size_t node) const
{
  return node > 0 && node <= nno ? values(node) : 0.0;
}


double SplineField2D::valueFE (const FiniteElement& fe) const
{
  if (!basis) return false;

  // Evaluate the basis functions at the given point
  Go::BasisPtsSf spline;
#pragma omp critical
  basis->computeBasis(fe.u,fe.v,spline);

  // Evaluate the solution field at the given point
  IntVec ip;
  ASMs2D::scatterInd(basis->numCoefs_u(),basis->numCoefs_v(),
		     basis->order_u(),basis->order_v(),
		     spline.left_idx,ip);

  Vector Vnod;
  utl::gather(ip,1,values,Vnod);
  return Vnod.dot(spline.basisValues);
}


double SplineField2D::valueCoor (const Vec3& x) const
{
  // use with caution
  Go::Point pt(3), clopt(3);
  pt[0] = x[0];
  pt[1] = x[1];
  pt[2] = x[2];
  double clo_u, clo_v, dist;
#pragma omp critical
  surf->closestPoint(pt, clo_u, clo_v, clopt, dist, 1e-5);

  FiniteElement fe;
  fe.u = clopt[0];
  fe.v = clopt[1];

  return valueFE(fe);
}


bool SplineField2D::valueGrid (RealArray& val, const int* npe) const
{
  val.clear();
  if (!basis) return false;

  // Compute parameter values of the visualization points
  std::array<RealArray,2> gpar;
  for (int dir = 0; dir < 2; dir++)
  {
    int nSegPerSpan = npe[dir] - 1;
    if (nSegPerSpan < 1) return false;

    RealArray::const_iterator uit = basis->basis(dir).begin();
    double ucurr = 0.0, uprev = *(uit++);
    while (uit != basis->basis(dir).end())
    {
      ucurr = *(uit++);
      if (ucurr > uprev)
        if (nSegPerSpan == 1)
          gpar[dir].push_back(uprev);
        else for (int i = 0; i < nSegPerSpan; i++)
        {
          double xg = (double)(2*i-nSegPerSpan)/(double)nSegPerSpan;
          gpar[dir].push_back(0.5*(ucurr*(1.0+xg) + uprev*(1.0-xg)));
        }
      uprev = ucurr;
    }

    if (ucurr > gpar[dir].back())
      gpar[dir].push_back(ucurr);
  }

  // Evaluate the field in the visualization points
  val.reserve(gpar[0].size()*gpar[1].size());
  for (size_t j = 0; j < gpar[1].size(); j++)
    for (size_t i = 0; i < gpar[0].size(); i++)
    {
      Go::BasisPtsSf spline;
#pragma omp critical
      basis->computeBasis(gpar[0][i],gpar[1][j],spline);

      IntVec ip;
      ASMs2D::scatterInd(basis->numCoefs_u(),basis->numCoefs_v(),
			 basis->order_u(),basis->order_v(),
			 spline.left_idx,ip);

      Vector Vnod;
      utl::gather(ip,1,values,Vnod);
      val.push_back(Vnod.dot(spline.basisValues));
    }

  return true;
}


bool SplineField2D::gradFE (const FiniteElement& fe, Vector& grad) const
{
  if (!basis) return false;
  if (!surf)  return false;

  // Evaluate the basis functions at the given point
  Go::BasisDerivsSf spline;
#pragma omp critical
  surf->computeBasis(fe.u,fe.v,spline);

  const int uorder = surf->order_u();
  const int vorder = surf->order_v();
  const size_t nen = uorder*vorder;

  Matrix dNdu(nen,2), dNdX;
  for (size_t n = 1; n <= nen; n++)
  {
    dNdu(n,1) = spline.basisDerivs_u[n-1];
    dNdu(n,2) = spline.basisDerivs_v[n-1];
  }

  IntVec ip;
  ASMs2D::scatterInd(surf->numCoefs_u(),surf->numCoefs_v(),
		     uorder,vorder,spline.left_idx,ip);

  // Evaluate the Jacobian inverse
  Matrix Xnod, Jac;
  Vector Xctrl(&(*surf->coefs_begin()),surf->coefs_end()-surf->coefs_begin());
  utl::gather(ip,surf->dimension(),Xctrl,Xnod);
  utl::Jacobian(Jac,dNdX,Xnod,dNdu);

  // Evaluate the gradient of the solution field at the given point
  if (basis != surf)
  {
    // Mixed formulation, the solution uses a different basis than the geometry
#pragma omp critical
    basis->computeBasis(fe.u,fe.v,spline);

    const size_t nbf = basis->order_u()*basis->order_v();
    dNdu.resize(nbf,2);
    for (size_t n = 1; n <= nbf; n++)
    {
      dNdu(n,1) = spline.basisDerivs_u[n-1];
      dNdu(n,2) = spline.basisDerivs_v[n-1];
    }
    dNdX.multiply(dNdu,Jac); // dNdX = dNdu * Jac

    ip.clear();
    ASMs2D::scatterInd(basis->numCoefs_u(),basis->numCoefs_v(),
		       basis->order_u(),basis->order_v(),
		       spline.left_idx,ip);
  }

  Vector Vnod;
  utl::gather(ip,1,values,Vnod);
  return dNdX.multiply(Vnod,grad,true); // grad = dNdX * Vnod^t
}


bool SplineField2D::hessianFE(const FiniteElement& fe, Matrix& H) const
{
  if (!basis) return false;
  if (!surf)  return false;

  // Order of basis
  const int uorder = surf->order_u();
  const int vorder = surf->order_v();
  const size_t nen = uorder*vorder;

  // Evaluate the basis functions at the given point
  Go::BasisDerivsSf  spline;
  Go::BasisDerivsSf2 spline2;
  Matrix3D d2Ndu2;
  Matrix dNdu, dNdX;
  IntVec ip;
  if (surf == basis) {
#pragma omp critical
    surf->computeBasis(fe.u,fe.v,spline2);

    dNdu.resize(nen,2);
    d2Ndu2.resize(nen,2,2);
    for (size_t n = 1; n <= nen; n++) {
      dNdu(n,1) = spline2.basisDerivs_u[n-1];
      dNdu(n,2) = spline2.basisDerivs_v[n-1];
      d2Ndu2(n,1,1) = spline2.basisDerivs_uu[n-1];
      d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = spline2.basisDerivs_uv[n-1];
      d2Ndu2(n,2,2) = spline2.basisDerivs_vv[n-1];
    }

    ASMs2D::scatterInd(surf->numCoefs_u(),surf->numCoefs_v(),
		       uorder,vorder,spline2.left_idx,ip);
  }
  else {
#pragma omp critical
    surf->computeBasis(fe.u,fe.v,spline);

    dNdu.resize(nen,2);
    for (size_t n = 1; n <= nen; n++) {
      dNdu(n,1) = spline.basisDerivs_u[n-1];
      dNdu(n,2) = spline.basisDerivs_v[n-1];
    }

    ASMs2D::scatterInd(surf->numCoefs_u(),surf->numCoefs_v(),
		       uorder,vorder,spline.left_idx,ip);
  }

  // Evaluate the Jacobian inverse
  Matrix Xnod, Jac;
  Vector Xctrl(&(*surf->coefs_begin()),surf->coefs_end()-surf->coefs_begin());
  utl::gather(ip,surf->dimension(),Xctrl,Xnod);
  utl::Jacobian(Jac,dNdX,Xnod,dNdu);

  // Evaluate the gradient of the solution field at the given point
  if (basis != surf)
  {
    // Mixed formulation, the solution uses a different basis than the geometry
#pragma omp critical
    basis->computeBasis(fe.u,fe.v,spline2);

    const size_t nbf = basis->order_u()*basis->order_v();
    dNdu.resize(nbf,2);
    d2Ndu2.resize(nbf,2,2);
    for (size_t n = 1; n <= nbf; n++) {
      dNdu(n,1) = spline2.basisDerivs_u[n-1];
      dNdu(n,2) = spline2.basisDerivs_v[n-1];
      d2Ndu2(n,1,1) = spline2.basisDerivs_uu[n-1];
      d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = spline2.basisDerivs_uv[n-1];
      d2Ndu2(n,2,2) = spline2.basisDerivs_vv[n-1];
    }

    ip.clear();
    ASMs2D::scatterInd(basis->numCoefs_u(),basis->numCoefs_v(),
		       basis->order_u(),basis->order_v(),
		       spline2.left_idx,ip);
  }

  Vector Vnod;
  utl::gather(ip,1,values,Vnod);
  return H.multiply(d2Ndu2,Vnod);
}


bool SplineField2D::gradCoor (const Vec3& x, Vector& grad) const
{
  // Not implemented yet
  return false;
}
