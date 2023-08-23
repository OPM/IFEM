// $Id$
//==============================================================================
//!
//! \file LRSplineFields3D.C
//!
//! \date Mar 9 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for LR spline-based finite element vector fields in 3D.
//!
//==============================================================================

#include "LRSpline/LRSplineVolume.h"

#include "LRSplineFields3D.h"
#include "LRSplineField.h"

#include "ASMu3D.h"
#include "ItgPoint.h"
#include "Vec3.h"


LRSplineFields3D::LRSplineFields3D (const ASMu3D* patch,
                                    const RealArray& v, char nbasis,
                                    int nnf, const char* name)
  : Fields(name),
    basis(patch->getBasis(nbasis)),
    vol(patch->getBasis(ASM::GEOMETRY_BASIS))
{
  nno = basis->nBasisFunctions();
  nelm = basis->nElements();

  size_t ofs = 0;
  for (char i = 1; i < nbasis; ++i)
    ofs += patch->getNoNodes(i)*patch->getNoFields(i);
  auto vit = v.begin()+ofs;

  if (nnf == 0)
    nnf = 3;
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


LRSplineFields3D::LRSplineFields3D (const LR::LRSplineVolume* svol,
                                    const RealArray& v, int cmp, const char* name)
  : Fields(name), basis(svol), vol(svol)
{
  values = v;
  nf = cmp;
}


bool LRSplineFields3D::valueNode (size_t node, Vector& vals) const
{
  if (node < 1 || node > nno) return false;

  vals.resize(nf);
  vals.fill(values.ptr()+(node-1)*nf);
  return true;
}


bool LRSplineFields3D::valueFE (const ItgPoint& x, Vector& vals) const
{
  if (!basis) return false;

  // Evaluate the basis functions at the given point
  int iel = basis->getElementContaining(x.u,x.v,x.w);
  auto elm = basis->getElement(iel);

  Go::BasisPts spline;
  basis->computeBasis(x.u,x.v,x.w,spline,iel);

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


bool LRSplineFields3D::valueCoor (const Vec4& x, Vector& vals) const
{
  if (x.u)
    return this->valueFE(ItgPoint(x.u[0],x.u[1],x.u[2]),vals);

  return false;
}


bool LRSplineFields3D::gradFE (const ItgPoint& x, Matrix& grad) const
{
  if (!basis) return false;
  if (!vol)   return false;

  // Evaluate the basis functions at the given point
  Matrix Xnod, Jac, dNdX;
  const LR::Element* elm;
  if (!LRSplineField::evalMapping(*vol,x,elm,Xnod,Jac,dNdX))
    return false;

  // Evaluate the gradient of the solution field at the given point
  Matrix Vnod;
  if (basis != vol)
    if (!LRSplineField::evalBasis(*basis,x,elm,Xnod,Jac,dNdX))
      return false;

  Vnod.resize(nf, elm->nBasisFunctions());
  size_t i = 1;
  for (const LR::Basisfunction* f : elm->support()) {
    for (size_t j = 1; j <= nf; ++j)
      Vnod(j,i) = values(f->getId()*nf+j);
    ++i;
  }

  grad.multiply(Vnod,dNdX); // grad = Vnod * dNdX

  return true;
}


bool LRSplineFields3D::hessianFE (const ItgPoint& x, Matrix3D& H) const
{
  if (!basis) return false;
  if (!vol)  return false;

  // Evaluate the basis functions at the given point
  Matrix Xnod, Jac, dNdX;
  const LR::Element* elm;
  Matrix3D Hess, d2NdX2;
  if (!LRSplineField::evalMapping(*vol,x,elm,Xnod,Jac,dNdX,&d2NdX2,&Hess))
    return false;

  if (vol != basis)
    if (!LRSplineField::evalBasis(*basis,x,elm,Xnod,Jac,dNdX,&d2NdX2,&Hess))
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
