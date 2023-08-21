// $Id$
//==============================================================================
//!
//! \file LRSplineFields2Dmx.C
//!
//! \date Mar 9 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for mixed LR spline-based finite element vector fields in 2D.
//!
//==============================================================================

#include "LRSpline/LRSplineSurface.h"

#include "LRSplineFields2Dmx.h"
#include "LRSplineField.h"

#include "ASMu2Dmx.h"
#include "ItgPoint.h"
#include "Utilities.h"
#include "Vec3.h"


LRSplineFields2Dmx::LRSplineFields2Dmx (const ASMu2Dmx* patch,
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


bool LRSplineFields2Dmx::valueNode (size_t node, Vector& vals) const
{
  return false;
}


bool LRSplineFields2Dmx::valueCoor (const Vec4& x, Vector& vals) const
{
  if (x.u)
    return this->valueFE(ItgPoint(x.u[0],x.u[1]),vals);

  std::cerr << "** LRSplineFields2Dmx::valueCoor: "
            << "not implemented without parameters\n";

  return false;
}


bool LRSplineFields2Dmx::valueFE (const ItgPoint& x, Vector& vals) const
{
  if (!surf) return false;

  vals.resize(2);

  // Evaluate the basis functions at the given point
  size_t ofs = 0;
  auto rit = vals.begin();
  for (int b : bases) {
    const LR::LRSplineSurface* basis = surf->getBasis(b);

    int iel = basis->getElementContaining(x.u,x.v);
    const LR::Element* elm = basis->getElement(iel);
    Go::BasisPtsSf spline;
    if (surf->rational())
      ASMu2D::computeBasisNurbs(x.u,x.v,spline,iel,*basis);
    else
      basis->computeBasis(x.u,x.v,spline,iel);

    // Evaluate the solution field at the given point
    size_t i = 1;
    Vector Vnod(elm->nBasisFunctions());
    for (const LR::Basisfunction* f : elm->support())
      Vnod(i++) = values(f->getId()+ofs+1);

    *rit++ = Vnod.dot(spline.basisValues);
    ofs += surf->getNoNodes(b);
  }

  return true;
}


bool LRSplineFields2Dmx::gradFE (const ItgPoint& x, Matrix& grad) const
{
  if (!surf)  return false;

  // Evaluate the basis functions at the given point
  Matrix Xnod, Jac, dNdX;
  const LR::Element* elm;
  const LR::LRSplineSurface* geo = surf->getBasis(ASMmxBase::elmBasis);
  if (!geo)
    geo = surf->getBasis(1);
  if (!LRSplineField::evalMapping(*geo,x,elm,Xnod,Jac,dNdX,surf->rational()))
    return false;

  // Evaluate the gradient of the solution field at the given point
  size_t ofs = 0;
  size_t row = 1;
  grad.resize(2,2);
  for (int b : bases) {
    const LR::LRSplineSurface* basis = surf->getBasis(b);
    if (!LRSplineField::evalBasis(*basis,x,elm,Xnod,Jac,dNdX,surf->rational()))
      return false;

    size_t i = 1;
    Matrix Vnod(1,elm->nBasisFunctions());
    for (const LR::Basisfunction* f : elm->support())
      Vnod(1,i++) = values(f->getId()+ofs+1);

    Matrix grad2;
    grad2.multiply(Vnod,dNdX); // grad = Xnod * dNdX
    grad(row,1) = grad2(1,1);
    grad(row++,2) = grad2(1,2);
    ofs += surf->getNoNodes(b)*surf->getNoFields(b);
  }

  return true;
}


bool LRSplineFields2Dmx::hessianFE (const ItgPoint& x, Matrix3D& H) const
{
  if (!surf)  return false;

  // Evaluate the basis functions at the given point
  Matrix Xnod, Jac, dNdX;
  Matrix3D d2NdX2, Hess;
  const LR::Element* elm;
  const LR::LRSplineSurface* geo = surf->getBasis(ASMmxBase::elmBasis);
  if (!geo)
    geo = surf->getBasis(1);
  if (!LRSplineField::evalMapping(*geo,x,elm,Xnod,Jac,
                                  dNdX,surf->rational(),&d2NdX2,&Hess))
    return false;

  H.resize(2,2,2);
  size_t ofs = 0;
  for (int b : bases) {
    const LR::LRSplineSurface* basis = surf->getBasis(b);
    if (!LRSplineField::evalBasis(*basis,x,elm,Xnod,Jac,
                                  dNdX,surf->rational(),&d2NdX2,&Hess))
      return false;

    size_t i = 1;
    Matrix Vnod(1, elm->nBasisFunctions());
    for (const LR::Basisfunction* f : elm->support())
      Vnod(1,i++) = values(f->getId()+ofs+1);

    Matrix3D hess(1,2,2);
    hess.multiply(Vnod,d2NdX2);
    for (size_t i = 1; i <= 2; ++i)
      for (size_t j = 1; j <= 2; ++j)
        H(b,i,j) = hess(1,i,j);
    ofs += surf->getNoNodes(b)*surf->getNoFields(b);
  }

  return true;
}
