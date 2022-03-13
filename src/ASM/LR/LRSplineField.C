// $Id$
//==============================================================================
//!
//! \file LRSplineField.C
//!
//! \date Mar 22 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Utility class for LR spline-based finite element fields.
//!
//==============================================================================

#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/LRSplineVolume.h"

#include "LRSplineField.h"
#include "CoordinateMapping.h"
#include "ItgPoint.h"
#include "SplineUtils.h"


bool LRSplineField::evalMapping (const LR::LRSplineSurface& surf,
                                 const ItgPoint& x,
                                 const LR::Element*& elm,
                                 Matrix& Xnod,
                                 Matrix& Jac,
                                 Matrix& dNdX,
                                 Matrix3D* d2NdX2,
                                 Matrix3D* Hess)
{
  int iel = surf.getElementContaining(x.u,x.v);
  elm = surf.getElement(iel);
  Vector N;
  Matrix dNdu;
  Matrix3D d2Ndu2;
  if (Hess) {
    Go::BasisDerivsSf2 spline2;
    surf.computeBasis(x.u,x.v,spline2,iel);
    SplineUtils::extractBasis(spline2, N, dNdu, d2Ndu2);
  } else {
    Go::BasisDerivsSf spline;
    surf.computeBasis(x.u,x.v,spline,iel);
    SplineUtils::extractBasis(spline, N, dNdu);
  }

  Xnod.resize(2,elm->nBasisFunctions());
  size_t i = 1;
  for (const LR::Basisfunction* f : elm->support()) {
    for (size_t j = 1; j <= 2; ++j)
      Xnod(j,i) = f->cp(j-1);
    ++i;
  }

  // Evaluate the Jacobian inverse
  if (!utl::Jacobian(Jac,dNdX,Xnod,dNdu))
    return false;

  return Hess ? utl::Hessian(*Hess,*d2NdX2,Jac,Xnod,d2Ndu2,dNdX,true) : true;
}


bool LRSplineField::evalBasis (const LR::LRSplineSurface& surf,
                               const ItgPoint& x,
                               const LR::Element*& elm,
                               const Matrix& Xnod,
                               const Matrix& Jac,
                               Matrix& dNdX,
                               Matrix3D* d2NdX2,
                               Matrix3D* Hess)
{
  int iel = surf.getElementContaining(x.u,x.v);
  elm = surf.getElement(iel);
  Vector N;
  Matrix dNdu;
  Matrix3D d2Ndu2;
  if (Hess) {
    Go::BasisDerivsSf2 spline2;
    surf.computeBasis(x.u,x.v,spline2,iel);
    SplineUtils::extractBasis(spline2, N, dNdu, d2Ndu2);
  } else {
    Go::BasisDerivsSf spline;
    surf.computeBasis(x.u,x.v,spline,iel);
    SplineUtils::extractBasis(spline, N, dNdu);
  }

  dNdX.multiply(dNdu,Jac); // dNdX = dNdu * Jac

  return Hess ? utl::Hessian(*Hess,*d2NdX2,Jac,Xnod,d2Ndu2,dNdX,false) : true;
}


bool LRSplineField::evalMapping (const LR::LRSplineVolume& vol,
                                 const ItgPoint& x,
                                 const LR::Element*& elm,
                                 Matrix& Xnod,
                                 Matrix& Jac,
                                 Matrix& dNdX,
                                 Matrix3D* d2NdX2,
                                 Matrix3D* Hess)
{
  int iel = vol.getElementContaining(x.u,x.v,x.w);
  elm = vol.getElement(iel);
  Vector N;
  Matrix dNdu;
  Matrix3D d2Ndu2;
  if (Hess) {
    Go::BasisDerivs2 spline2;
    vol.computeBasis(x.u,x.v,x.w,spline2,iel);
    SplineUtils::extractBasis(spline2, N, dNdu, d2Ndu2);
  } else {
    Go::BasisDerivs spline;
    vol.computeBasis(x.u,x.v,x.w,spline,iel);
    SplineUtils::extractBasis(spline, N, dNdu);
  }

  Xnod.resize(3,elm->nBasisFunctions());
  size_t i = 1;
  for (const LR::Basisfunction* f : elm->support()) {
    for (size_t j = 1; j <= 3; ++j)
      Xnod(j,i) = f->cp(j-1);
    ++i;
  }

  // Evaluate the Jacobian inverse
  if (!utl::Jacobian(Jac,dNdX,Xnod,dNdu))
    return false;

  return Hess ? utl::Hessian(*Hess,*d2NdX2,Jac,Xnod,d2Ndu2,dNdX,true) : true;
}


bool LRSplineField::evalBasis (const LR::LRSplineVolume& vol,
                               const ItgPoint& x,
                               const LR::Element*& elm,
                               const Matrix& Xnod,
                               const Matrix& Jac,
                               Matrix& dNdX,
                               Matrix3D* d2NdX2,
                               Matrix3D* Hess)
{
  int iel = vol.getElementContaining(x.u,x.v,x.w);
  elm = vol.getElement(iel);
  Vector N;
  Matrix dNdu;
  Matrix3D d2Ndu2;
  if (Hess) {
    Go::BasisDerivs2 spline2;
    vol.computeBasis(x.u,x.v,x.w,spline2,iel);
    SplineUtils::extractBasis(spline2, N, dNdu, d2Ndu2);
  } else {
    Go::BasisDerivs spline;
    vol.computeBasis(x.u,x.v,x.w,spline,iel);
    SplineUtils::extractBasis(spline, N, dNdu);
  }

  dNdX.multiply(dNdu,Jac); // dNdX = dNdu * Jac

  return Hess ? utl::Hessian(*Hess,*d2NdX2,Jac,Xnod,d2Ndu2,dNdX,false) : true;
}
