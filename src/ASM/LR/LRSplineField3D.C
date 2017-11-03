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

#include "LRSplineField3D.h"
#include "ASMu3D.h"
#include "FiniteElement.h"
#include "CoordinateMapping.h"
#include "Vec3.h"

#include "LRSpline/LRSplineVolume.h"


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
}


double LRSplineField3D::valueNode (size_t node) const
{
  return node > 0 && node <= nno ? values(node) : 0.0;
}


double LRSplineField3D::valueFE (const FiniteElement& fe) const
{
  if (!basis) return 0.0;

  // Evaluate the basis functions at the given point
  int iel = basis->getElementContaining(fe.u,fe.v,fe.w);
  auto elm = basis->getElement(iel);
  Go::BasisPts spline;
  basis->computeBasis(fe.u,fe.v,fe.w,spline,iel);

  Vector Vnod;
  Vnod.reserve(elm->nBasisFunctions());
  for (auto it  = elm->constSupportBegin();
            it != elm->constSupportEnd(); ++it)
    Vnod.push_back(values((*it)->getId()+1));

  return Vnod.dot(spline.basisValues);
}


double LRSplineField3D::valueCoor (const Vec4& x) const
{
  assert(0);
  return 0.0;
}


bool LRSplineField3D::gradFE (const FiniteElement& fe, Vector& grad) const
{
  if (!basis) return false;
  if (!vol)   return false;

  // Evaluate the basis functions at the given point
  int iel = vol->getElementContaining(fe.u,fe.v,fe.w);
  auto elm = vol->getElement(iel);
  Go::BasisDerivs spline;
  vol->computeBasis(fe.u,fe.v,fe.w,spline,iel);

  const size_t nen = elm->nBasisFunctions();

  Matrix dNdu(nen,3), dNdX;
  for (size_t n = 1; n <= nen; n++)
  {
    dNdu(n,1) = spline.basisDerivs_u[n-1];
    dNdu(n,2) = spline.basisDerivs_v[n-1];
    dNdu(n,3) = spline.basisDerivs_w[n-1];
  }

  Matrix Xnod(3,nen), Jac;
  Vector Vnod;
  size_t i = 1;
  Vnod.reserve(nen);
  for (auto it  = elm->constSupportBegin();
            it != elm->constSupportEnd(); ++it, ++i) {
    for (size_t j = 1; j <= 3; ++j)
      Xnod(j, i) = (*it)->cp(j-1);

    if (vol == basis)
      Vnod.push_back(values((*it)->getId()+1));
  }

  // Evaluate the Jacobian inverse
  utl::Jacobian(Jac,dNdX,Xnod,dNdu);

  // Evaluate the gradient of the solution field at the given point
  if (basis != vol)
  {
    // Mixed formulation, the solution uses a different basis than the geometry
    int iel = basis->getElementContaining(fe.u,fe.v,fe.w);
    auto belm = basis->getElement(iel);
    basis->computeBasis(fe.u,fe.v,fe.w,spline,iel);

    const size_t nbf = belm->nBasisFunctions();
    dNdu.resize(nbf,3);
    for (size_t n = 1; n <= nbf; n++)
    {
      dNdu(n,1) = spline.basisDerivs_u[n-1];
      dNdu(n,2) = spline.basisDerivs_v[n-1];
      dNdu(n,3) = spline.basisDerivs_w[n-1];
    }
    dNdX.multiply(dNdu,Jac); // dNdX = dNdu * Jac

    Vnod.reserve(nbf);
    for (auto it  = belm->constSupportBegin();
              it != belm->constSupportEnd(); ++it)
      Vnod.push_back(values((*it)->getId()+1));
  }

  return dNdX.multiply(Vnod,grad,true); // grad = dNdX * Vnod^t
}


bool LRSplineField3D::hessianFE (const FiniteElement& fe, Matrix& H) const
{
  if (!basis) return false;
  if (!vol)  return false;

  int iel = vol->getElementContaining(fe.u,fe.v,fe.w);
  auto elm = vol->getElement(iel);
  const size_t nen = elm->nBasisFunctions();

  // Evaluate the basis functions at the given point
  Go::BasisDerivs  spline;
  Go::BasisDerivs2 spline2;
  Matrix3D d2Ndu2;
  Matrix dNdu(nen,3), dNdX;
  Matrix Xnod(nen,3), Jac;
  Vector Vnod;
  size_t i = 1;
  for (auto it  = elm->constSupportBegin();
            it != elm->constSupportEnd(); ++it, ++i)
    for (size_t j = 1; j <= 3; ++j)
      Xnod(i, j) = (*it)->cp(j-1);

  if (vol == basis) {
    vol->computeBasis(fe.u,fe.v,fe.w,spline2,iel);
    d2Ndu2.resize(nen,3,3);
    for (size_t n = 1; n <= nen; n++) {
      dNdu(n,1) = spline2.basisDerivs_u[n-1];
      dNdu(n,2) = spline2.basisDerivs_v[n-1];
      dNdu(n,3) = spline2.basisDerivs_w[n-1];
      d2Ndu2(n,1,1) = spline2.basisDerivs_uu[n-1];
      d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = spline2.basisDerivs_uv[n-1];
      d2Ndu2(n,1,3) = d2Ndu2(n,3,1) = spline2.basisDerivs_uw[n-1];
      d2Ndu2(n,2,2) = spline2.basisDerivs_vv[n-1];
      d2Ndu2(n,2,3) = d2Ndu2(n,3,2) = spline2.basisDerivs_vw[n-1];
      d2Ndu2(n,3,3) = spline2.basisDerivs_ww[n-1];
    }

    Vnod.reserve(nen);
    for (auto it  = elm->constSupportBegin();
              it != elm->constSupportEnd(); ++it)
      Vnod.push_back(values((*it)->getId()+1));
  }
  else {
    vol->computeBasis(fe.u,fe.v,fe.w,spline,iel);
    for (size_t n = 1; n <= nen; n++) {
      dNdu(n,1) = spline.basisDerivs_u[n-1];
      dNdu(n,2) = spline.basisDerivs_v[n-1];
      dNdu(n,3) = spline.basisDerivs_w[n-1];
    }
  }

  // Evaluate the Jacobian inverse
  utl::Jacobian(Jac,dNdX,Xnod,dNdu);

  // Evaluate the gradient of the solution field at the given point
  if (basis != vol) {
    // Mixed formulation, the solution uses a different basis than the geometry
    int iel = basis->getElementContaining(fe.u,fe.v,fe.w);
    auto belm = basis->getElement(iel);
    basis->computeBasis(fe.u,fe.v,fe.w,spline2,iel);

    const size_t nbf = belm->nBasisFunctions();
    dNdu.resize(nbf,3);
    d2Ndu2.resize(nbf,3,3);
    for (size_t n = 1; n <= nbf; n++) {
      dNdu(n,1) = spline2.basisDerivs_u[n-1];
      dNdu(n,2) = spline2.basisDerivs_v[n-1];
      dNdu(n,3) = spline2.basisDerivs_w[n-1];
      d2Ndu2(n,1,1) = spline2.basisDerivs_uu[n-1];
      d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = spline2.basisDerivs_uv[n-1];
      d2Ndu2(n,1,3) = d2Ndu2(n,3,1) = spline2.basisDerivs_uw[n-1];
      d2Ndu2(n,2,2) = spline2.basisDerivs_vv[n-1];
      d2Ndu2(n,2,3) = d2Ndu2(n,3,2) = spline2.basisDerivs_vw[n-1];
      d2Ndu2(n,3,3) = spline2.basisDerivs_ww[n-1];
    }

    Vnod.reserve(nbf);
    for (auto it  = belm->constSupportBegin();
              it != belm->constSupportEnd(); ++it)
      Vnod.push_back(values((*it)->getId()+1));
  }

  return H.multiply(d2Ndu2,Vnod);
}
