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

#include "LRSplineField2D.h"
#include "ASMu2D.h"
#include "FiniteElement.h"
#include "CoordinateMapping.h"
#include "Vec3.h"

#include "LRSpline/LRSplineSurface.h"

#include <cassert>


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


double LRSplineField2D::valueNode (size_t node) const
{
  return node > 0 && node <= nno ? values(node) : 0.0;
}


double LRSplineField2D::valueFE (const FiniteElement& fe) const
{
  if (!basis) return 0.0;

  // Evaluate the basis functions at the given point
  int iel = basis->getElementContaining(fe.u,fe.v);
  auto elm = basis->getElement(iel);
  Go::BasisPtsSf spline;
  basis->computeBasis(fe.u,fe.v,spline,iel);

  Vector Vnod;
  Vnod.reserve(elm->nBasisFunctions());
  for (auto it  = elm->constSupportBegin();
            it != elm->constSupportEnd(); ++it)
    Vnod.push_back(values((*it)->getId()+1));

  return Vnod.dot(spline.basisValues);
}


double LRSplineField2D::valueCoor (const Vec4& x) const
{
  if (x.u) {
    FiniteElement fe;
    fe.u = x.u[0];
    fe.v = x.u[1];

    return this->valueFE(fe);
  }

  assert(0);
  return 0.0;
}


bool LRSplineField2D::gradFE (const FiniteElement& fe, Vector& grad) const
{
  if (!basis) return false;
  if (!surf)  return false;

  // Evaluate the basis functions at the given point
  int iel = surf->getElementContaining(fe.u,fe.v);
  auto elm = surf->getElement(iel);
  Go::BasisDerivsSf spline;
  basis->computeBasis(fe.u,fe.v,spline,iel);

  const size_t nen = elm->nBasisFunctions();

  Matrix dNdu(nen,2), dNdX;
  for (size_t n = 1; n <= nen; n++)
  {
    dNdu(n,1) = spline.basisDerivs_u[n-1];
    dNdu(n,2) = spline.basisDerivs_v[n-1];
  }

  Matrix Xnod(2,nen), Jac;
  Vector Vnod;
  size_t i = 1;
  Vnod.reserve(nen);
  for (auto it  = elm->constSupportBegin();
            it != elm->constSupportEnd(); ++it, ++i) {
    for (size_t j = 1; j <= 2; ++j)
      Xnod(j, i) = (*it)->cp(j-1);

    if (surf == basis)
      Vnod.push_back(values((*it)->getId()+1));
  }

  // Evaluate the Jacobian inverse
  utl::Jacobian(Jac,dNdX,Xnod,dNdu);

  // Evaluate the gradient of the solution field at the given point
  if (basis != surf)
  {
    // Mixed formulation, the solution uses a different basis than the geometry
    int iel = basis->getElementContaining(fe.u,fe.v);
    auto belm = basis->getElement(iel);
    basis->computeBasis(fe.u,fe.v,spline,iel);

    const size_t nbf = belm->nBasisFunctions();
    dNdu.resize(nbf,2);
    for (size_t n = 1; n <= nbf; n++)
    {
      dNdu(n,1) = spline.basisDerivs_u[n-1];
      dNdu(n,2) = spline.basisDerivs_v[n-1];
    }
    dNdX.multiply(dNdu,Jac); // dNdX = dNdu * Jac

    Vnod.reserve(nbf);
    for (auto it  = belm->constSupportBegin();
              it != belm->constSupportEnd(); ++it)
      Vnod.push_back(values((*it)->getId()+1));
  }

  return dNdX.multiply(Vnod,grad,true); // grad = dNdX * Vnod^t
}


bool LRSplineField2D::hessianFE (const FiniteElement& fe, Matrix& H) const
{
  if (!basis) return false;
  if (!surf)  return false;

  int iel = surf->getElementContaining(fe.u,fe.v);
  auto elm = surf->getElement(iel);
  const size_t nen = elm->nBasisFunctions();

  // Evaluate the basis functions at the given point
  Go::BasisDerivsSf  spline;
  Go::BasisDerivsSf2 spline2;
  Matrix3D d2Ndu2;
  Matrix dNdu(nen,2), dNdX;
  Matrix Xnod(nen,2), Jac;
  Vector Vnod;
  size_t i = 1;
  for (auto it  = elm->constSupportBegin();
            it != elm->constSupportEnd(); ++it, ++i)
    for (size_t j = 1; j <= 2; ++j)
      Xnod(i, j) = (*it)->cp(j-1);

  if (surf == basis) {
    surf->computeBasis(fe.u,fe.v,spline2,iel);
    d2Ndu2.resize(nen,2,2);
    for (size_t n = 1; n <= nen; n++) {
      dNdu(n,1) = spline2.basisDerivs_u[n-1];
      dNdu(n,2) = spline2.basisDerivs_v[n-1];
      d2Ndu2(n,1,1) = spline2.basisDerivs_uu[n-1];
      d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = spline2.basisDerivs_uv[n-1];
      d2Ndu2(n,2,2) = spline2.basisDerivs_vv[n-1];
    }

    Vnod.reserve(nen);
    for (auto it  = elm->constSupportBegin();
              it != elm->constSupportEnd(); ++it)
      Vnod.push_back(values((*it)->getId()+1));
  }
  else {
    surf->computeBasis(fe.u,fe.v,spline,iel);
    for (size_t n = 1; n <= nen; n++) {
      dNdu(n,1) = spline.basisDerivs_u[n-1];
      dNdu(n,2) = spline.basisDerivs_v[n-1];
    }
  }

  // Evaluate the Jacobian inverse
  utl::Jacobian(Jac,dNdX,Xnod,dNdu);

  // Evaluate the gradient of the solution field at the given point
  if (basis != surf) {
    // Mixed formulation, the solution uses a different basis than the geometry
    int iel = basis->getElementContaining(fe.u,fe.v);
    auto belm = basis->getElement(iel);
    basis->computeBasis(fe.u,fe.v,spline2,iel);

    const size_t nbf = belm->nBasisFunctions();
    dNdu.resize(nbf,2);
    d2Ndu2.resize(nbf,2,2);
    for (size_t n = 1; n <= nbf; n++) {
      dNdu(n,1) = spline2.basisDerivs_u[n-1];
      dNdu(n,2) = spline2.basisDerivs_v[n-1];
      d2Ndu2(n,1,1) = spline2.basisDerivs_uu[n-1];
      d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = spline2.basisDerivs_uv[n-1];
      d2Ndu2(n,2,2) = spline2.basisDerivs_vv[n-1];
    }

    Vnod.reserve(nbf);
    for (auto it  = belm->constSupportBegin();
              it != belm->constSupportEnd(); ++it)
      Vnod.push_back(values((*it)->getId()+1));
  }

  return H.multiply(d2Ndu2,Vnod);
}
