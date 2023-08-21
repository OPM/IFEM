// $Id$
//==============================================================================
//!
//! \file SplineFields2Dmx.C
//!
//! \date Oct 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for mixed spline-based finite element vector fields in 2D.
//!
//==============================================================================

#include "GoTools/geometry/SplineSurface.h"

#include "SplineFields2Dmx.h"
#include "SplineField.h"

#include "ASMs2Dmx.h"
#include "ItgPoint.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Vec3.h"


SplineFields2Dmx::SplineFields2Dmx (const ASMs2Dmx* patch,
                                    const RealArray& v, char basis,
                                    const char* name)
  : Fields(name), surf(patch), bases(utl::getDigits(basis))
{
  nf = 2;
  auto vit = v.begin();
  size_t ofs = 0;
  for (int i = 1; i < *bases.begin(); ++i)
    ofs += patch->getNoNodes(i)*patch->getNoFields(i);
  vit += ofs;
  for (int b : bases) {
    size_t nno = patch->getNoNodes(b)*patch->getNoFields(b);
    RealArray::const_iterator end = v.size() > nno+ofs ? vit+nno : v.end();
    std::copy(vit,end,std::back_inserter(values));
    vit += nno;
    ofs += nno;
    values.resize(ofs);
  }
}


bool SplineFields2Dmx::valueNode (size_t node, Vector& vals) const
{
  return false;
}


bool SplineFields2Dmx::valueCoor (const Vec4& x, Vector& vals) const
{
  if (x.u)
    return this->valueFE(ItgPoint(x.u[0],x.u[1]),vals);

  // Use with caution, very slow!
  Go::Point pt(x.x,x.y,x.z), clopt(3);
  double clo_u, clo_v, dist;
  const Go::SplineSurface* geo = surf->getBasis(ASMmxBase::elmBasis);
  if (!geo)
    geo = surf->getBasis(1);
#pragma omp critical
  geo->closestPoint(pt, clo_u, clo_v, clopt, dist, 1.0e-5);

  return this->valueFE(ItgPoint(clo_u,clo_v),vals);
}


bool SplineFields2Dmx::valueFE (const ItgPoint& x, Vector& vals) const
{
  if (!surf) return false;

  vals.resize(2);

  // Evaluate the basis functions at the given point
  auto vit = values.begin();
  auto rit = vals.begin();
  for (int b : bases) {
    Go::SplineSurface* basis = surf->getBasis(b);
    Go::BasisPtsSf spline;
#pragma omp critical
    basis->computeBasis(x.u,x.v,spline);

    // Evaluate the solution field at the given point
    std::vector<int> ip;
    ASMs2D::scatterInd(basis->numCoefs_u(),basis->numCoefs_v(),
                       basis->order_u(),basis->order_v(),
                       spline.left_idx,ip);

    Matrix Vnod;
    utl::gather(ip,1,Vector(&*vit,surf->getNoNodes(b)),Vnod);
    Vector val2;
    Vnod.multiply(spline.basisValues,val2); // vals = Vnod * basisValues
    *rit++ = val2.front();
    vit += surf->getNoNodes(b);
  }

  return true;
}


bool SplineFields2Dmx::gradFE (const ItgPoint& x, Matrix& grad) const
{
  if (!surf) return false;

  // Evaluate the basis functions at the given point
  Matrix Xnod, Jac, dNdX;
  std::vector<int> ip;
  const Go::SplineSurface* geo = surf->getBasis(ASMmxBase::elmBasis);
  if (!geo)
    geo = surf->getBasis(1);
  if (!SplineField::evalMapping(*geo,surf->getNoSpaceDim(),x,ip,Xnod,Jac,dNdX))
    return false;

  // Evaluate the gradient of the solution field at the given point
  auto vit = values.begin();
  size_t row = 1;
  grad.resize(2,2);
  for (int b : bases) {
    const Go::SplineSurface* basis = surf->getBasis(b);
    if (!SplineField::evalBasis(*basis,x,ip,Xnod,Jac,dNdX))
      return false;

    const size_t nval = surf->getNoNodes(b)*surf->getNoFields(b);
    Matrix Vnod;
    utl::gather(ip,1,Vector(&*vit,nval),Vnod);
    Matrix grad2;
    grad2.multiply(Vnod,dNdX); // grad = Vnod * dNdX
    grad(row,1) = grad2(1,1);
    grad(row++,2) = grad2(1,2);
    vit += nval;
  }

  return true;
}


bool SplineFields2Dmx::hessianFE (const ItgPoint& x, Matrix3D& H) const
{
  if (!surf)  return false;

  // Evaluate the basis functions at the given point
  Matrix Xnod, Jac, dNdX;
  Matrix3D d2NdX2, Hess;
  std::vector<int> ip;
  const Go::SplineSurface* geo = surf->getBasis(ASMmxBase::elmBasis);
  if (!geo)
    geo = surf->getBasis(1);
  if (!SplineField::evalMapping(*geo,surf->getNoSpaceDim(),x,ip,Xnod,Jac,dNdX,&d2NdX2,&Hess))
    return false;

  H.resize(2,2,2);
  auto vit = values.begin();
  for (int b : bases) {
    const Go::SplineSurface* basis = surf->getBasis(b);
    if (!SplineField::evalBasis(*basis,x,ip,Xnod,Jac,dNdX,&d2NdX2,&Hess))
      return false;

    const size_t nval = surf->getNoNodes(b)*surf->getNoFields(b);
    Matrix Vnod;
    utl::gather(ip,1,Vector(&*vit,nval),Vnod);

    Matrix3D hess(1,2,2);
    hess.multiply(Vnod,d2NdX2);
    for (size_t i = 1; i <= 2; ++i)
      for (size_t j = 1; j <= 2; ++j)
        H(b,i,j) = hess(1,i,j);

    vit += nval;
  }

  return true;
}
