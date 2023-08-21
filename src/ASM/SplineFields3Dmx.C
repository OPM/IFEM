// $Id$
//==============================================================================
//!
//! \file SplineFields3Dmx.C
//!
//! \date Oct 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for mixed spline-based finite element vector fields in 3D.
//!
//==============================================================================

#include "GoTools/trivariate/SplineVolume.h"

#include "SplineFields3Dmx.h"
#include "SplineField.h"

#include "ASMs3Dmx.h"
#include "ItgPoint.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Vec3.h"


SplineFields3Dmx::SplineFields3Dmx (const ASMs3Dmx* patch,
                                    const RealArray& v, char basis,
                                    const char* name)
  : Fields(name), svol(patch), bases(utl::getDigits(basis))
{
  nf = 3;
  auto vit = v.begin();
  size_t ofs = 0;
  for (int i = 1; i < *bases.begin(); ++i)
    ofs += patch->getNoNodes(i)*patch->getNoFields(i);
  vit += ofs;
  for (int b : bases) {
    size_t nno = patch->getNoFields(b)*patch->getNoNodes(b);
    RealArray::const_iterator end = v.size() > nno+ofs ? vit+nno : v.end();
    std::copy(vit,end,std::back_inserter(values));
    vit += nno;
    ofs += nno;
    values.resize(ofs);
  }
}


bool SplineFields3Dmx::valueNode (size_t node, Vector& vals) const
{
  return false;
}


bool SplineFields3Dmx::valueCoor (const Vec4& x, Vector& vals) const
{
  if (x.u)
    return this->valueFE(ItgPoint(x.u[0],x.u[1],x.u[2]),vals);

  // Use with caution, very slow!
  Go::Point pt(x.x,x.y,x.z), clopt(3);
  double clo_u, clo_v, clo_w, dist;
  const Go::SplineVolume* geo = svol->getBasis(ASMmxBase::elmBasis);
  if (!geo)
    geo = svol->getBasis(1);
#pragma omp critical
  geo->closestPoint(pt, clo_u, clo_v, clo_w, clopt, dist, 1.0e-5);

  return this->valueFE(ItgPoint(clo_u,clo_v,clo_w),vals);
}


bool SplineFields3Dmx::valueFE (const ItgPoint& x, Vector& vals) const
{
  if (!svol) return false;

  vals.resize(3);

  // Evaluate the basis functions at the given point
  auto vit = values.begin();
  auto rit = vals.begin();
  for (int b : bases) {
    Go::SplineVolume* basis = svol->getBasis(b);
    Go::BasisPts spline;
#pragma omp critical
    basis->computeBasis(x.u,x.v,x.w,spline);

    // Evaluate the solution field at the given point
    IntVec ip;
    ASMs3D::scatterInd(basis->numCoefs(0),basis->numCoefs(1),
                       basis->numCoefs(2),
                       basis->order(0),basis->order(1),basis->order(2),
                       spline.left_idx,ip);

    Matrix Vnod;
    utl::gather(ip,1,Vector(&*vit,svol->getNoNodes(b)),Vnod);
    Vector val2;
    Vnod.multiply(spline.basisValues,val2); // vals = Vnod * basisValues
    *rit++ = val2.front();
    vit += svol->getNoNodes(b);
  }

  return true;
}


bool SplineFields3Dmx::gradFE (const ItgPoint& x, Matrix& grad) const
{
  if (!svol)  return false;

  // Evaluate the basis functions at the given point
  Matrix Xnod, Jac, dNdX;
  std::vector<int> ip;
  const Go::SplineVolume* geo = svol->getBasis(ASMmxBase::elmBasis);
  if (!geo)
    geo = svol->getBasis(1);
  if (!SplineField::evalMapping(*geo,x,ip,Xnod,Jac,dNdX))
    return false;

  // Evaluate the gradient of the solution field at the given point
  auto vit = values.begin();
  size_t row = 1;
  grad.resize(3,3);
  for (int b : bases) {
    const Go::SplineVolume* basis = svol->getBasis(b);
    if (!SplineField::evalBasis(*basis,x,ip,Xnod,Jac,dNdX))
      return false;

    const size_t nval = svol->getNoNodes(b)*svol->getNoFields(b);
    utl::gather(ip,1,Vector(&*vit,nval),Xnod);
    Matrix grad2;
    grad2.multiply(Xnod,dNdX); // grad = Xnod * dNdX
    grad(row,1) = grad2(1,1);
    grad(row,2) = grad2(1,2);
    grad(row++,3) = grad2(1,3);
    vit += nval;
  }

  return true;
}


bool SplineFields3Dmx::hessianFE (const ItgPoint& x, Matrix3D& H) const
{
  if (!svol)  return false;

  Matrix Xnod, Jac, dNdX;
  Matrix3D d2NdX2, Hess;
  std::vector<int> ip;
  const Go::SplineVolume* geo = svol->getBasis(ASMmxBase::elmBasis);
  if (!geo)
    geo = svol->getBasis(1);
  if (!SplineField::evalMapping(*geo,x,ip,Xnod,Jac,dNdX,&d2NdX2,&Hess))
    return false;

  H.resize(3,3,3);
  auto vit = values.begin();
  for (int b : bases) {
    const Go::SplineVolume* basis = svol->getBasis(b);
    if (!SplineField::evalBasis(*basis,x,ip,Xnod,Jac,dNdX,&d2NdX2,&Hess))
      return false;

    const size_t nval = svol->getNoNodes(b)*svol->getNoFields(b);
    Matrix Vnod;
    utl::gather(ip,1,Vector(&*vit,nval),Vnod);

    Matrix3D hess(1,3,3);
    hess.multiply(Vnod,d2NdX2);
    for (size_t i = 1; i <= 3; ++i)
      for (size_t j = 1; j <= 3; ++j)
        H(b,i,j) = hess(1,i,j);

    vit += nval;
  }

  return true;
}
