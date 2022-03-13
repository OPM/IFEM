// $Id$
//==============================================================================
//!
//! \file LRSplineField3D.C
//!
//! \date Mar 9 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for LR spline-based finite element scalar field in 3D.
//!
//==============================================================================

#include "LRSpline/LRSplineVolume.h"

#include "LRSplineField3D.h"
#include "LRSplineField.h"

#include "ASMu3D.h"
#include "ItgPoint.h"
#include "SplineUtils.h"
#include "Vec3.h"


LRSplineField3D::LRSplineField3D (const ASMu3D* patch,
                                  const RealArray& v, char nbasis,
                                  char cmp, const char* name)
  : FieldBase(name), basis(patch->getBasis(nbasis)), vol(patch->getVolume())
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


double LRSplineField3D::valueNode (size_t node) const
{
  return node > 0 && node <= nno ? values(node) : 0.0;
}


double LRSplineField3D::valueFE (const ItgPoint& x) const
{
  if (!basis) return 0.0;

  // Evaluate the basis functions at the given point
  int iel = basis->getElementContaining(x.u,x.v,x.w);
  auto elm = basis->getElement(iel);
  Go::BasisPts spline;
  basis->computeBasis(x.u,x.v,x.w,spline,iel);

  Vector Vnod;
  Vnod.reserve(elm->nBasisFunctions());
  for (const LR::Basisfunction* f : elm->support())
    Vnod.push_back(values(f->getId()+1));

  return Vnod.dot(spline.basisValues);
}


double LRSplineField3D::valueCoor (const Vec4& x) const
{
  if (x.u)
    return this->valueFE(ItgPoint(x.u[0],x.u[1],x.u[2]));

  std::cerr << "** LRSplineField3D::valueCoor: "
            << "not implemented without parameters\n";

  return false;
}


bool LRSplineField3D::gradFE (const ItgPoint& x, Vector& grad) const
{
  if (!basis) return false;
  if (!vol)   return false;

  // Evaluate the basis functions at the given point
  Matrix Xnod, Jac, dNdX;
  const LR::Element* elm;
  if (!LRSplineField::evalMapping(*vol,x,elm,Xnod,Jac,dNdX))
    return false;

  // Evaluate the gradient of the solution field at the given point
  Vector Vnod;
  if (basis != vol)
    if (!LRSplineField::evalBasis(*basis,x,elm,Xnod,Jac,dNdX))
      return false;

  Vnod.reserve(elm->nBasisFunctions());
  for (const LR::Basisfunction* f : elm->support())
    Vnod.push_back(values(f->getId()+1));

  return dNdX.multiply(Vnod,grad,true); // grad = dNdX * Vnod^t
}


bool LRSplineField3D::hessianFE (const ItgPoint& x, Matrix& H) const
{
  if (!basis) return false;
  if (!vol)  return false;

  // Evaluate the basis functions at the given point
  Matrix Xnod, Jac, dNdX;
  const LR::Element* elm;
  Matrix3D d2NdX2, Hess;
  if (!LRSplineField::evalMapping(*vol,x,elm,Xnod,Jac,dNdX,&d2NdX2,&Hess))
    return false;

  if (vol != basis)
    if (!LRSplineField::evalBasis(*vol,x,elm,Xnod,Jac,dNdX,&d2NdX2,&Hess))
      return false;

  Matrix Vnod(1,elm->nBasisFunctions());
  size_t i = 1;
  for (const LR::Basisfunction* f : elm->support())
    Vnod(1,i++) = values(f->getId()+1);

  Matrix3D hess(1,3,3);
  hess.multiply(Vnod,d2NdX2);

  H.resize(3,3);
  for (size_t i = 1; i <= 3; ++i)
    for (size_t j = 1; j <= 3; ++j)
      H(i,j) = hess(1,i,j);

  return true;
}
