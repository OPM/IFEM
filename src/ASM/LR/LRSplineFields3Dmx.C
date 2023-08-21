// $Id$
//==============================================================================
//!
//! \file LRSplineFields3Dmx.C
//!
//! \date Mar 9 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for LR spline-based finite element mixed vector fields in 3D.
//!
//==============================================================================

#include "LRSpline/LRSplineVolume.h"

#include "LRSplineFields3Dmx.h"
#include "LRSplineField.h"

#include "ASMu3Dmx.h"
#include "ItgPoint.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Vec3.h"


LRSplineFields3Dmx::LRSplineFields3Dmx (const ASMu3Dmx* patch,
                                        const RealArray& v, char basis,
                                        const char* name)
  : Fields(name), vol(patch), bases(utl::getDigits(basis))
{
  nf = 3;
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


bool LRSplineFields3Dmx::valueNode (size_t node, Vector& vals) const
{
  return false;
}


bool LRSplineFields3Dmx::valueCoor (const Vec4& x, Vector& vals) const
{
  if (x.u)
    return this->valueFE(ItgPoint(x.u[0],x.u[1],x.u[2]),vals);

  std::cerr << "** LRSplineFields3Dmx::valueCoor: "
            << "not implemented without parameters\n";

  return false;
}


bool LRSplineFields3Dmx::valueFE (const ItgPoint& x, Vector& vals) const
{
  if (!vol) return false;

  vals.resize(3);

  // Evaluate the basis functions at the given point
  size_t ofs = 0;
  auto rit = vals.begin();
  for (int b : bases) {
    const LR::LRSplineVolume* basis = vol->getBasis(b);

    int iel = basis->getElementContaining(x.u,x.v,x.w);
    const LR::Element* elm = basis->getElement(iel);
    Go::BasisPts spline;
    basis->computeBasis(x.u,x.v,x.w,spline,iel);

    // Evaluate the solution field at the given point

    size_t i = 1;
    Vector Vnod(elm->nBasisFunctions());
    for (const LR::Basisfunction* f : elm->support())
      Vnod(i++) = values(f->getId()+ofs+1);

    *rit++ = Vnod.dot(spline.basisValues);
    ofs += vol->getNoNodes(b);
  }

  return true;
}


bool LRSplineFields3Dmx::gradFE (const ItgPoint& x, Matrix& grad) const
{
  if (!vol)  return false;

  // Evaluate the basis functions at the given point
  Matrix Xnod, Jac, dNdX;
  const LR::LRSplineVolume* geo = vol->getBasis(ASMmxBase::itgBasis);
  if (!geo)
    geo = vol->getBasis(1);
  const LR::Element* elm;
  if (!LRSplineField::evalMapping(*geo,x,elm,Xnod,Jac,dNdX))
    return false;

  // Evaluate the gradient of the solution field at the given point
  size_t ofs = 0;
  size_t row = 1;
  grad.resize(3,3);
  for (int b : bases) {
    const LR::LRSplineVolume* basis = vol->getBasis(b);
    if (!LRSplineField::evalBasis(*basis,x,elm,Xnod,Jac,dNdX))
      return false;

    size_t i = 1;
    Matrix Vnod(1,elm->nBasisFunctions());
    for (const LR::Basisfunction* f : elm->support())
      Vnod(1,i++) = values(f->getId()+ofs+1);

    Matrix grad2;
    grad2.multiply(Vnod,dNdX); // grad = Xnod * dNdX
    grad(row,1) = grad2(1,1);
    grad(row,2) = grad2(1,2);
    grad(row++,3) = grad2(1,3);
    ofs += vol->getNoNodes(b)*vol->getNoFields(b);
  }

  return true;
}


bool LRSplineFields3Dmx::hessianFE (const ItgPoint& x, Matrix3D& H) const
{
  if (!vol)  return false;

  // Evaluate the basis functions at the given point
  Matrix Xnod, Jac, dNdX;
  Matrix3D d2NdX2, Hess;
  const LR::LRSplineVolume* geo = vol->getBasis(ASMmxBase::itgBasis);
  if (!geo)
    geo = vol->getBasis(1);
  const LR::Element* elm;
  if (!LRSplineField::evalMapping(*geo,x,elm,Xnod,Jac,dNdX,&d2NdX2,&Hess))
    return false;

  H.resize(3,3,3);
  size_t ofs = 0;
  for (int b : bases) {
    const LR::LRSplineVolume* basis = vol->getBasis(b);
    if (!LRSplineField::evalBasis(*basis,x,elm,Xnod,Jac,dNdX,&d2NdX2,&Hess))
      return false;

    size_t i = 1;
    Matrix Vnod(1, elm->nBasisFunctions());
    for (const LR::Basisfunction* f : elm->support())
      Vnod(1,i++) = values(f->getId()+ofs+1);

    Matrix3D hess(1,3,3);
    hess.multiply(Vnod,d2NdX2);
    for (size_t i = 1; i <= 3; ++i)
      for (size_t j = 1; j <= 3; ++j)
        H(b,i,j) = hess(1,i,j);
    ofs += vol->getNoNodes(b)*vol->getNoFields(b);
  }

  return true;
}
