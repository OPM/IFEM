// $Id$
//==============================================================================
//!
//! \file LRSplineFields2D.C
//!
//! \date Mar 9 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for LR spline-based finite element vector fields in 2D.
//!
//==============================================================================

#include "LRSpline/LRSplineSurface.h"

#include "LRSplineFields2D.h"
#include "LRSplineField.h"

#include "ASMu2D.h"
#include "ItgPoint.h"
#include "SplineUtils.h"
#include "Vec3.h"


LRSplineFields2D::LRSplineFields2D (const ASMu2D* patch,
                                    const RealArray& v, char nbasis,
                                    int nnf, const char* name)
  : Fields(name), basis(patch->getBasis(nbasis)), surf(patch->getSurface())
{
  nno = basis->nBasisFunctions();
  nelm = basis->nElements();

  size_t ofs = 0;
  for (char i = 1; i < nbasis; ++i)
    ofs += patch->getNoNodes(i)*patch->getNoFields(i);
  auto vit = v.begin()+ofs;

  if (nnf == 0)
    nnf = 2;
  nf = nnf;
  int nfc = patch->getNoFields(nbasis);
  values.resize(nno*nf);
  int ndof = nfc*nno;
  auto end = v.size() > ofs+ndof ? vit+ndof : v.end();
  if (nfc == nf)
    std::copy(vit,end,values.begin());
  else
    for (size_t i = 0; i < nno && vit != end; ++i, vit += nfc)
      for (size_t j = 0; j < nf; ++j)
        values[nf*i+j] = *(vit+j);

  basis->generateIDs();
}


LRSplineFields2D::LRSplineFields2D (const LR::LRSplineSurface* srf,
                                    const RealArray& v, int cmp, const char* name)
  : Fields(name), basis(srf), surf(srf)
{
  values = v;
  nf = cmp;
}


bool LRSplineFields2D::valueNode (size_t node, Vector& vals) const
{
  if (node < 1 || node > nno) return false;

  vals.resize(nf);
  vals.fill(values.ptr()+(node-1)*nf);
  return true;
}


bool LRSplineFields2D::valueFE (const ItgPoint& x, Vector& vals) const
{
  if (!basis) return false;

  // Evaluate the basis functions at the given point
  int iel = basis->getElementContaining(x.u,x.v);
  auto elm = basis->getElement(iel);

  Go::BasisPtsSf spline;
  basis->computeBasis(x.u,x.v,spline,iel);

  // Evaluate the solution field at the given point
  Matrix Vnod(nf, elm->nBasisFunctions());
  size_t i = 1;
  for (const LR::Basisfunction* f : elm->support()) {
    for (size_t j = 1; j <= nf; ++j)
      Vnod(j,i) = values(f->getId()*nf+j);
    ++i;
  }

  Vnod.multiply(spline.basisValues,vals); // vals = Vnod * basisValues

  return true;
}


bool LRSplineFields2D::valueCoor (const Vec4& x, Vector& vals) const
{
  if (x.u)
    return this->valueFE(ItgPoint(x.u[0],x.u[1]),vals);

  return false;
}


bool LRSplineFields2D::gradFE (const ItgPoint& x, Matrix& grad) const
{
  if (!basis) return false;
  if (!surf)  return false;

  // Evaluate the basis functions at the given point
  Matrix Xnod, Jac, dNdX;
  const LR::Element* elm;
  if (!LRSplineField::evalMapping(*surf,x,elm,Xnod,Jac,dNdX))
    return false;

  // Evaluate the gradient of the solution field at the given point
  Matrix Vnod;
  if (basis != surf)
    if (!LRSplineField::evalBasis(*surf,x,elm,Xnod,Jac,dNdX))
      return false;

  Vnod.resize(nf, elm->nBasisFunctions());
  size_t i = 1;
  for (const LR::Basisfunction* f : elm->support()) {
    for (size_t j = 1; j <= nf; ++j)
      Vnod(j,i) = values(f->getId()*nf+j);
    ++i;
  }

  return !grad.multiply(Vnod,dNdX).empty(); // grad = Vnod * dNdX
}


bool LRSplineFields2D::hessianFE (const ItgPoint& x, Matrix3D& H) const
{
  if (!basis) return false;
  if (!surf)  return false;

  // Evaluate the basis functions at the given point
  Matrix Xnod, Jac, dNdX;
  const LR::Element* elm;
  Matrix3D d2NdX2, Hess;
  if (!LRSplineField::evalMapping(*surf,x,elm,Xnod,Jac,dNdX,&d2NdX2,&Hess))
    return false;

  // Evaluate the gradient of the solution field at the given point
  if (surf != basis)
    if (!LRSplineField::evalBasis(*surf,x,elm,Xnod,Jac,dNdX,&d2NdX2,&Hess))
      return false;

  Matrix Vnod(nf, elm->nBasisFunctions());
  size_t i = 1;
  for (const LR::Basisfunction* f : elm->support()) {
    for (size_t j = 1; j <= nf; ++j)
      Vnod(j,i) = values(f->getId()*nf+j);
    ++i;
  }

  return H.multiply(Vnod,d2NdX2);
}
