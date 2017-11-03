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

#include "LRSplineFields3Dmx.h"
#include "ASMu3Dmx.h"
#include "FiniteElement.h"
#include "CoordinateMapping.h"
#include "Utilities.h"

#include "LRSpline/LRSplineVolume.h"


LRSplineFields3Dmx::LRSplineFields3Dmx (const ASMu3Dmx* patch,
                                        const RealArray& v, char basis,
                                        const char* name)
  : Fields(name), vol(patch)
{
  bases = utl::getDigits(basis);
  nf = 3;
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


bool LRSplineFields3Dmx::valueNode (size_t node, Vector& vals) const
{
  return false;
}


bool LRSplineFields3Dmx::valueFE (const FiniteElement& fe, Vector& vals) const
{
  if (!vol) return false;

  vals.resize(3);

  // Evaluate the basis functions at the given point
  size_t ofs = 0;
  auto rit = vals.begin();
  for (const auto& it : bases) {
    const LR::LRSplineVolume* basis = vol->getBasis(it);

    int iel = basis->getElementContaining(fe.u,fe.v,fe.w);
    auto elm = basis->getElement(iel);
    Go::BasisPts spline;
    basis->computeBasis(fe.u,fe.v,fe.w,spline,iel);

    // Evaluate the solution field at the given point

    size_t i = 1;
    Vector Vnod(elm->nBasisFunctions());
    for (auto eit  = elm->constSupportBegin();
              eit != elm->constSupportEnd(); ++eit, ++i)
      Vnod(i) = values((*eit)->getId()+ofs+1);

    *rit++ = Vnod.dot(spline.basisValues);
    ofs += vol->getNoNodes(it);
  }

  return true;
}


bool LRSplineFields3Dmx::gradFE (const FiniteElement& fe, Matrix& grad) const
{
  if (!vol)  return false;

  // Evaluate the basis functions at the given point
  Go::BasisDerivs spline;
  const LR::LRSplineVolume* gvol = vol->getBasis(ASMmxBase::geoBasis);
  int iel = gvol->getElementContaining(fe.u,fe.v,fe.w);
  auto elm = gvol->getElement(iel);
  gvol->computeBasis(fe.u,fe.v,fe.w,spline,iel);

  const size_t nen = elm->nBasisFunctions();
  Matrix dNdu(nen,3), dNdX;
  for (size_t n = 1; n <= nen; n++)
  {
    dNdu(n,1) = spline.basisDerivs_u[n-1];
    dNdu(n,2) = spline.basisDerivs_v[n-1];
    dNdu(n,3) = spline.basisDerivs_w[n-1];
  }

  Matrix Xnod(nen,3), Jac;
  size_t i = 1;
  for (auto it  = elm->constSupportBegin();
            it != elm->constSupportEnd(); ++it, ++i)
    for (size_t j = 1; j <= 3; ++j)
      Xnod(i, j) = (*it)->cp(j-1);

  // Evaluate the Jacobian inverse
  utl::Jacobian(Jac,dNdX,Xnod,dNdu);

  // Evaluate the gradient of the solution field at the given point
  size_t ofs = 0;
  size_t row = 1;
  for (const auto& it : bases) {
    const LR::LRSplineVolume* basis = vol->getBasis(it);
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

    size_t i = 1;
    Matrix Vnod(nbf,1);
    for (auto eit  = elm->constSupportBegin();
              eit != elm->constSupportEnd(); ++eit, ++i)
      Vnod(i,1) = values((*eit)->getId()+ofs+1);

    Matrix grad2;
    grad2.multiply(Vnod,dNdX); // grad = Xnod * dNdX
    grad(row,1) = grad2(1,1);
    grad(row,2) = grad2(1,2);
    grad(row++,3) = grad2(1,3);
    ofs += vol->getNoNodes(it)*vol->getNoFields(it);
  }

  return true;
}
