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

#include "LRSplineFields2Dmx.h"
#include "ASMu2Dmx.h"
#include "FiniteElement.h"
#include "CoordinateMapping.h"
#include "Utilities.h"
#include "Vec3.h"

#include "LRSpline/LRSplineSurface.h"


LRSplineFields2Dmx::LRSplineFields2Dmx (const ASMu2Dmx* patch,
                                        const RealArray& v, char basis,
                                        const char* name)
  : Fields(name), surf(patch)
{
  bases = utl::getDigits(basis);
  nf = 2;
  auto vit = v.begin();
  size_t ofs = 0;
  for (int i = 1; i < *bases.begin(); ++i)
    ofs += patch->getNoNodes(i)*patch->getNoFields(i);
  vit += ofs;
  for (const auto& it : bases) {
    size_t nno = patch->getNoNodes(it)*patch->getNoFields(it);
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


bool LRSplineFields2Dmx::valueFE (const FiniteElement& fe, Vector& vals) const
{
  if (!surf) return false;

  vals.resize(2);

  // Evaluate the basis functions at the given point
  size_t ofs = 0;
  auto rit = vals.begin();
  for (const auto& it : bases) {
    const LR::LRSplineSurface* basis = surf->getBasis(it);

    int iel = basis->getElementContaining(fe.u,fe.v);
    auto elm = basis->getElement(iel);
    Go::BasisPtsSf spline;
    basis->computeBasis(fe.u,fe.v,spline,iel);

    // Evaluate the solution field at the given point

    size_t i = 1;
    Vector Vnod(elm->nBasisFunctions());
    for (auto it  = elm->constSupportBegin();
        it != elm->constSupportEnd(); ++it, ++i)
      Vnod(i) = values((*it)->getId()+ofs+1);

    *rit++ = Vnod.dot(spline.basisValues);
    ofs += surf->getNoNodes(it);
  }

  return true;
}


bool LRSplineFields2Dmx::valueCoor (const Vec3& x, Vector& vals) const
{
  // Not implemented yet
  return false;
}


bool LRSplineFields2Dmx::gradFE (const FiniteElement& fe, Matrix& grad) const
{
  if (!surf)  return false;

  // Evaluate the basis functions at the given point
  Go::BasisDerivsSf spline;
  const LR::LRSplineSurface* gsurf = surf->getBasis(ASMmxBase::geoBasis);
  int iel = gsurf->getElementContaining(fe.u,fe.v);
  auto elm = gsurf->getElement(iel);
  gsurf->computeBasis(fe.u,fe.v,spline,iel);

  const size_t nen = elm->nBasisFunctions();

  Matrix dNdu(nen,2), dNdX;
  for (size_t n = 1; n <= nen; n++)
  {
    dNdu(n,1) = spline.basisDerivs_u[n-1];
    dNdu(n,2) = spline.basisDerivs_v[n-1];
  }

  Matrix Xnod(nen,2), Jac;
  size_t i = 1;
  for (auto it  = elm->constSupportBegin();
            it != elm->constSupportEnd(); ++it, ++i)
    for (size_t j = 1; j <= 2; ++j)
      Xnod(i, j) = (*it)->cp(j-1);

  // Evaluate the Jacobian inverse
  utl::Jacobian(Jac,dNdX,Xnod,dNdu);

  // Evaluate the gradient of the solution field at the given point
  size_t ofs = 0;
  size_t row = 1;
  for (const auto& it : bases) {
    const LR::LRSplineSurface* basis = surf->getBasis(it);
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

    size_t i = 1;
    Matrix Vnod(nbf,1);
    for (auto it  = elm->constSupportBegin();
        it != elm->constSupportEnd(); ++it, ++i)
      Vnod(i,1) = values((*it)->getId()+ofs+1);

    Matrix grad2;
    grad2.multiply(Vnod,dNdX); // grad = Xnod * dNdX
    grad(row,1) = grad2(1,1);
    grad(row++,2) = grad2(1,2);
    ofs += surf->getNoNodes(it)*surf->getNoFields(it);
  }

  return true;
}


bool LRSplineFields2Dmx::hessianFE(const FiniteElement& fe, Matrix3D& H) const
{
  return false;
}


bool LRSplineFields2Dmx::gradCoor (const Vec3& x, Matrix& grad) const
{
  // Not implemented yet
  return false;
}
