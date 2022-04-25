// $Id$
//==============================================================================
//!
//! \file LRSplineField2D.C
//!
//! \date Mar 9 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for LR spline-based finite element scalar field in 2D.
//!
//==============================================================================

#include "LRSpline/LRSplineSurface.h"

#include "LRSplineField2D.h"
#include "LRSplineField.h"

#include "ASMu2D.h"
#include "ItgPoint.h"
#include "SplineUtils.h"
#include "Vec3.h"


LRSplineField2D::LRSplineField2D (const ASMu2D* patch,
                                  const RealArray& v, char nbasis,
                                  char cmp, const char* name)
  : FieldBase(name), basis(patch->getBasis(nbasis)), surf(patch->getSurface())
{
  nno = basis->nBasisFunctions();
  nelm = basis->nElements();

  size_t ofs = 0;
  for (char i = 1; i < nbasis; ++i)
    ofs += patch->getNoNodes(i)*patch->getNoFields(i);
  auto vit = v.begin()+ofs;

  values.resize(nno);
  int nf = patch->getNoFields(nbasis);
  int ndof = nf > 1 && cmp > 0 ? nf*nno : nno;
  auto end = v.size() > ofs+ndof ? vit+ndof : v.end();
  if (nf == 1 || cmp == 0)
    std::copy(vit,end,values.begin());
  else
    for (size_t i = 0; i < nno && vit != end; ++i, vit += nf)
      values[i] = *(vit+cmp-1);

  basis->generateIDs();
}


LRSplineField2D::LRSplineField2D (const LR::LRSplineSurface* srf,
                                  const RealArray& v, const char* name)
  : FieldBase(name), basis(srf), surf(srf)
{
  values = v;
}


double LRSplineField2D::valueNode (size_t node) const
{
  return node > 0 && node <= nno ? values(node) : 0.0;
}


double LRSplineField2D::valueFE (const ItgPoint& x) const
{
  if (!basis) return 0.0;

  // Evaluate the basis functions at the given point
  int iel = basis->getElementContaining(x.u,x.v);
  auto elm = basis->getElement(iel);
  Go::BasisPtsSf spline;
  basis->computeBasis(x.u,x.v,spline,iel);

  Vector Vnod;
  Vnod.reserve(elm->nBasisFunctions());
  for (const LR::Basisfunction* f : elm->support())
    Vnod.push_back(values(f->getId()+1));

  return Vnod.dot(spline.basisValues);
}


double LRSplineField2D::valueCoor (const Vec4& x) const
{
  // Just produce a segmentation fault, if invoked without the parameters
  if (x.u)
    return this->valueFE(ItgPoint(x.u[0],x.u[1]));

  std::cerr << "** LRSplineField2D::valueCoor: "
            << "not implemented without parameters\n";

  return false;
}


bool LRSplineField2D::gradFE (const ItgPoint& x, Vector& grad) const
{
  if (!basis) return false;
  if (!surf)  return false;

  // Evaluate the basis functions at the given point
  Matrix Xnod, Jac, dNdX;
  const LR::Element* elm;
  if (!LRSplineField::evalMapping(*surf,x,elm,Xnod,Jac,dNdX))
    return false;

  // Evaluate the gradient of the solution field at the given point
  Vector Vnod;
  if (basis != surf)
    if (!LRSplineField::evalBasis(*basis,x,elm,Xnod,Jac,dNdX))
      return false;

  Vnod.reserve(elm->nBasisFunctions());
  for (const LR::Basisfunction* f : elm->support())
    Vnod.push_back(values(f->getId()+1));

  return dNdX.multiply(Vnod,grad,true); // grad = dNdX * Vnod^t
}


bool LRSplineField2D::hessianFE (const ItgPoint& x, Matrix& H) const
{
  if (!basis) return false;
  if (!surf)  return false;

  // Evaluate the basis functions at the given point
  Matrix Xnod, Jac, dNdX;
  const LR::Element* elm;
  Matrix3D d2NdX2, Hess;
  if (!LRSplineField::evalMapping(*surf,x,elm,Xnod,Jac,dNdX,&d2NdX2,&Hess))
    return false;

  if (surf != basis)
    if (!LRSplineField::evalBasis(*basis,x,elm,Xnod,Jac,dNdX,&d2NdX2,&Hess))
      return false;

  Matrix Vnod(1,elm->nBasisFunctions());
  size_t i = 1;
  for (const LR::Basisfunction* f : elm->support())
    Vnod(1,i++) = values(f->getId()+1);

  Matrix3D hess(1,2,2);
  hess.multiply(Vnod,d2NdX2);

  H.resize(2,2);
  for (size_t i = 1; i <= 2; ++i)
    for (size_t j = 1; j <= 2; ++j)
      H(i,j) = hess(1,i,j);

  return true;
}
