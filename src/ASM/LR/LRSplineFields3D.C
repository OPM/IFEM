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

#include "LRSplineFields3D.h"
#include "ASMu3D.h"
#include "FiniteElement.h"
#include "CoordinateMapping.h"
#include "Vec3.h"

#include "LRSpline/LRSplineVolume.h"

#include <cassert>


LRSplineFields3D::LRSplineFields3D (const ASMu3D* patch,
                                    const RealArray& v, char nbasis,
                                    int nnf, const char* name)
  : Fields(name), basis(patch->getBasis(nbasis)), vol(patch->getVolume())
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


bool LRSplineFields3D::valueFE (const FiniteElement& fe, Vector& vals) const
{
  if (!basis) return false;

  // Evaluate the basis functions at the given point
  int iel = basis->getElementContaining(fe.u,fe.v,fe.w);
  auto elm = basis->getElement(iel);

  Go::BasisPts spline;
  basis->computeBasis(fe.u,fe.v,fe.w,spline,iel);

  // Evaluate the solution field at the given point
  Matrix Vnod(nf, elm->nBasisFunctions());
  size_t i = 1;
  for (auto it  = elm->constSupportBegin();
            it != elm->constSupportEnd(); ++it, ++i)
    for (size_t j = 1; j <= nf; ++j)
      Vnod(j,i) = values((*it)->getId()*nf+j);

  Vnod.multiply(spline.basisValues,vals); // vals = Vnod * basisValues

  return true;
}


bool LRSplineFields3D::valueCoor (const Vec4& x, Vector& vals) const
{
  if (x.u) {
    FiniteElement fe;
    fe.u = x.u[0];
    fe.v = x.u[1];
    fe.w = x.u[2];

    return this->valueFE(fe, vals);
  }

  assert(0);
  return false;
}


bool LRSplineFields3D::gradFE (const FiniteElement& fe, Matrix& grad) const
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

  Matrix Xnod(3, nen), Jac;
  size_t i = 1;
  for (auto it  = elm->constSupportBegin();
            it != elm->constSupportEnd(); ++it, ++i)
    for (size_t j = 1; j <= 3; ++j)
      Xnod(j,i) = (*it)->cp(j-1);

  // Evaluate the Jacobian inverse
  utl::Jacobian(Jac,dNdX,Xnod,dNdu);

  // Evaluate the gradient of the solution field at the given point
  Matrix Vnod;
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

    Vnod.resize(nf, nbf);
    i = 1;
    for (auto it  = belm->constSupportBegin();
              it != belm->constSupportEnd(); ++it, ++i)
      for (size_t j = 1; j <= nf; ++j)
        Vnod(j,i) = values((*it)->getId()*nf+j);
  }
  else {
    Vnod.resize(nf, nen);
    i = 1;
    for (auto it  = elm->constSupportBegin();
              it != elm->constSupportEnd(); ++it, ++i)
      for (size_t j = 1; j <= nf; ++j)
        Vnod(j,i) = values((*it)->getId()*nf+j);
  }

  grad.multiply(Vnod,dNdX); // grad = Vnod * dNdX

  return true;
}


bool LRSplineFields3D::hessianFE (const FiniteElement& fe, Matrix3D& H) const
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
  Matrix dNdu, dNdX;
  Matrix Xnod(3, nen), Jac, Vnod;
  size_t i = 1;
  for (auto it  = elm->constSupportBegin();
            it != elm->constSupportEnd(); ++it, ++i)
    for (size_t j = 1; j <= 3; ++j)
      Xnod(j,i) = (*it)->cp(j-1);

  if (vol == basis) {
    vol->computeBasis(fe.u,fe.v,fe.w,spline2,iel);

    dNdu.resize(nen,3);
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

    Vnod.resize(nf, nen);
    i = 1;
    for (auto it  = elm->constSupportBegin();
              it != elm->constSupportEnd(); ++it, ++i)
      for (size_t j = 1; j <= nf; ++j)
        Vnod(j,i) = values((*it)->getId()*nf+j);
  }
  else {
    vol->computeBasis(fe.u,fe.v,fe.w,spline,iel);

    dNdu.resize(nen,3);
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

    Vnod.resize(nf, nbf);
    i = 1;
    for (auto it  = belm->constSupportBegin();
              it != belm->constSupportEnd(); ++it, ++i)
      for (size_t j = 1; j <= nf; ++j)
        Vnod(j,i) = values((*it)->getId()*nf+j);
  }

  return H.multiply(Vnod,d2Ndu2);
}
