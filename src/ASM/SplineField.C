// $Id$
//==============================================================================
//!
//! \file SplineField.C
//!
//! \date Mar 15 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Utility class for spline-based finite element fields.
//!
//==============================================================================

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"

#include "SplineField.h"

#include "ASMs2D.h"
#include "ASMs3D.h"
#include "CoordinateMapping.h"
#include "ItgPoint.h"
#include "SplineUtils.h"


bool SplineField::evalMapping (const Go::SplineSurface& surf,
                               int nsd,
                               const ItgPoint& x,
                               std::vector<int>& ip,
                               Matrix& Xnod,
                               Matrix& Jac,
                               Matrix& dNdX,
                               Matrix3D* d2NdX2,
                               Matrix3D* Hess)
{
  Vector N;
  Matrix dNdu;
  Matrix3D d2Ndu2;

  if (Hess) {
    Go::BasisDerivsSf2 spline2;
#pragma omp critical
    surf.computeBasis(x.u,x.v,spline2);
    SplineUtils::extractBasis(spline2, N, dNdu, d2Ndu2);
    ASMs2D::scatterInd(surf.numCoefs_u(),surf.numCoefs_v(),
                       surf.order_u(),surf.order_v(),
                       spline2.left_idx,ip);

  } else {
    Go::BasisDerivsSf spline;
#pragma omp critical
    surf.computeBasis(x.u,x.v,spline);
    SplineUtils::extractBasis(spline, N, dNdu);
    ASMs2D::scatterInd(surf.numCoefs_u(),surf.numCoefs_v(),
                       surf.order_u(),surf.order_v(),
                       spline.left_idx,ip);
  }

  Xnod.resize(nsd,ip.size());
  for (size_t i = 0; i < ip.size(); i++)
    Xnod.fillColumn(1+i,&(*surf.coefs_begin())+surf.dimension()*ip[i]);

  if (!utl::Jacobian(Jac,dNdX,Xnod,dNdu))
    return false;

  return Hess ? utl::Hessian(*Hess,*d2NdX2,Jac,Xnod,d2Ndu2,dNdX,true) : true;
}


bool SplineField::evalBasis (const Go::SplineSurface& surf,
                             const ItgPoint& x,
                             std::vector<int>& ip,
                             const Matrix& Xnod,
                             const Matrix& Jac,
                             Matrix& dNdX,
                             Matrix3D* d2NdX2,
                             Matrix3D* Hess)
{
  Vector N;
  Matrix dNdu;
  Matrix3D d2Ndu2;

  ip.clear();
  if (Hess) {
    Go::BasisDerivsSf2 spline2;
#pragma omp critical
    surf.computeBasis(x.u,x.v,spline2);
    SplineUtils::extractBasis(spline2, N, dNdu, d2Ndu2);
    ASMs2D::scatterInd(surf.numCoefs_u(),surf.numCoefs_v(),
                       surf.order_u(),surf.order_v(),
                       spline2.left_idx,ip);
  } else {
    Go::BasisDerivsSf spline;
#pragma omp critical
    surf.computeBasis(x.u,x.v,spline);
    SplineUtils::extractBasis(spline, N, dNdu);
    ASMs2D::scatterInd(surf.numCoefs_u(),surf.numCoefs_v(),
                       surf.order_u(),surf.order_v(),
                       spline.left_idx,ip);
  }

  dNdX.multiply(dNdu,Jac); // dNdX = dNdu * Jac

  return Hess ? utl::Hessian(*Hess,*d2NdX2,Jac,Xnod,d2Ndu2,dNdX,false) : true;
}


bool SplineField::evalMapping (const Go::SplineVolume& vol,
                               const ItgPoint& x,
                               std::vector<int>& ip,
                               Matrix& Xnod,
                               Matrix& Jac,
                               Matrix& dNdX,
                               Matrix3D* d2NdX2,
                               Matrix3D* Hess)
{
  Vector N;
  Matrix dNdu;
  Matrix3D d2Ndu2;

  if (Hess) {
    Go::BasisDerivs2 spline2;
#pragma omp critical
    vol.computeBasis(x.u,x.v,x.w,spline2);
    SplineUtils::extractBasis(spline2, N, dNdu, d2Ndu2);
    ASMs3D::scatterInd(vol.numCoefs(0),vol.numCoefs(1),vol.numCoefs(2),
                       vol.order(0),vol.order(1),vol.order(2),
                       spline2.left_idx,ip);

  } else {
    Go::BasisDerivs spline;
#pragma omp critical
    vol.computeBasis(x.u,x.v,x.w,spline);
    SplineUtils::extractBasis(spline, N, dNdu);
    ASMs3D::scatterInd(vol.numCoefs(0),vol.numCoefs(1),vol.numCoefs(2),
                       vol.order(0),vol.order(1),vol.order(2),
                       spline.left_idx,ip);
  }

  Xnod.resize(3,ip.size());
  for (size_t i = 0; i < ip.size(); i++)
    Xnod.fillColumn(1+i,&(*vol.coefs_begin())+vol.dimension()*ip[i]);

  if (!utl::Jacobian(Jac,dNdX,Xnod,dNdu))
    return false;

  return Hess ? utl::Hessian(*Hess,*d2NdX2,Jac,Xnod,d2Ndu2,dNdX,true) : true;
}


bool SplineField::evalBasis (const Go::SplineVolume& vol,
                             const ItgPoint& x,
                             std::vector<int>& ip,
                             const Matrix& Xnod,
                             const Matrix& Jac,
                             Matrix& dNdX,
                             Matrix3D* d2NdX2,
                             Matrix3D* Hess)
{
  Vector N;
  Matrix dNdu;
  Matrix3D d2Ndu2;

  ip.clear();
  if (Hess) {
    Go::BasisDerivs2 spline2;
#pragma omp critical
    vol.computeBasis(x.u,x.v,x.w,spline2);
    SplineUtils::extractBasis(spline2, N, dNdu, d2Ndu2);
    ASMs3D::scatterInd(vol.numCoefs(0),vol.numCoefs(1),vol.numCoefs(2),
                       vol.order(0),vol.order(1),vol.order(2),
                       spline2.left_idx,ip);
  } else {
    Go::BasisDerivs spline;
#pragma omp critical
    vol.computeBasis(x.u,x.v,x.w,spline);
    SplineUtils::extractBasis(spline, N, dNdu);
    ASMs3D::scatterInd(vol.numCoefs(0),vol.numCoefs(1),vol.numCoefs(2),
                       vol.order(0),vol.order(1),vol.order(2),
                       spline.left_idx,ip);
  }

  dNdX.multiply(dNdu,Jac); // dNdX = dNdu * Jac

  return Hess ? utl::Hessian(*Hess,*d2NdX2,Jac,Xnod,d2Ndu2,dNdX,false) : true;
}
