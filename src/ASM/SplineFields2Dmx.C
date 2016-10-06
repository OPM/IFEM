// $Id$
//==============================================================================
//!
//! \file SplineFields2Dmx.C
//!
//! \date Oct 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for mixed spline-based finite element vector fields in 2D.
//!
//==============================================================================

#include "SplineFields2Dmx.h"
#include "ASMs2Dmx.h"
#include "FiniteElement.h"
#include "CoordinateMapping.h"
#include "Utilities.h"
#include "Vec3.h"

#include "GoTools/geometry/SplineSurface.h"


SplineFields2Dmx::SplineFields2Dmx (const ASMs2Dmx* patch,
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


bool SplineFields2Dmx::valueNode (size_t node, Vector& vals) const
{
  return false;
}


bool SplineFields2Dmx::valueFE (const FiniteElement& fe, Vector& vals) const
{
  if (!surf) return false;

  vals.resize(2);

  // Evaluate the basis functions at the given point
  auto vit = values.begin();
  auto rit = vals.begin();
  for (const auto& it : bases) {
    Go::SplineSurface* basis = surf->getBasis(it);
    Go::BasisPtsSf spline;
#pragma omp critical
    basis->computeBasis(fe.u,fe.v,spline);

    // Evaluate the solution field at the given point
    std::vector<int> ip;
    ASMs2D::scatterInd(basis->numCoefs_u(),basis->numCoefs_v(),
                       basis->order_u(),basis->order_v(),
                       spline.left_idx,ip);

    Matrix Vnod;
    utl::gather(ip,1,Vector(&*vit,surf->getNoNodes(it)),Vnod);
    Vector val2;
    Vnod.multiply(spline.basisValues,val2); // vals = Vnod * basisValues
    *rit++ = val2.front();
    vit += surf->getNoNodes(it);
  }

  return true;
}


bool SplineFields2Dmx::valueCoor (const Vec3& x, Vector& vals) const
{
  // Not implemented yet
  return false;
}


bool SplineFields2Dmx::gradFE (const FiniteElement& fe, Matrix& grad) const
{
  if (!surf)  return false;

  // Evaluate the basis functions at the given point
  Go::BasisDerivsSf spline;
  const Go::SplineSurface* gsurf = surf->getBasis(ASMmxBase::geoBasis);
#pragma omp critical
  gsurf->computeBasis(fe.u,fe.v,spline);

  const int uorder = gsurf->order_u();
  const int vorder = gsurf->order_v();
  const size_t nen = uorder*vorder;

  Matrix dNdu(nen,2), dNdX;
  for (size_t n = 1; n <= nen; n++)
  {
    dNdu(n,1) = spline.basisDerivs_u[n-1];
    dNdu(n,2) = spline.basisDerivs_v[n-1];
  }

  std::vector<int> ip;
  ASMs2D::scatterInd(gsurf->numCoefs_u(),gsurf->numCoefs_v(),
		     uorder,vorder,spline.left_idx,ip);

  // Evaluate the Jacobian inverse
  Matrix Xnod, Jac;
  Vector Xctrl(&(*gsurf->coefs_begin()),gsurf->coefs_end()-gsurf->coefs_begin());
  utl::gather(ip,gsurf->dimension(),Xctrl,Xnod);
  utl::Jacobian(Jac,dNdX,Xnod,dNdu);

  // Evaluate the gradient of the solution field at the given point
  auto vit = values.begin();
  size_t row = 1;
  for (const auto& it : bases) {
    const Go::SplineSurface* basis = surf->getBasis(it);
#pragma omp critical
    basis->computeBasis(fe.u,fe.v,spline);

    const size_t nbf = basis->order_u()*basis->order_v();
    dNdu.resize(nbf,2);
    for (size_t n = 1; n <= nbf; n++)
    {
      dNdu(n,1) = spline.basisDerivs_u[n-1];
      dNdu(n,2) = spline.basisDerivs_v[n-1];
    }
    dNdX.multiply(dNdu,Jac); // dNdX = dNdu * Jac

    ip.clear();
    ASMs2D::scatterInd(basis->numCoefs_u(),basis->numCoefs_v(),
		       basis->order_u(),basis->order_v(),
		       spline.left_idx,ip);

    utl::gather(ip,1,Vector(&*vit, surf->getNoNodes(it)*
                                   surf->getNoFields(it)),Xnod);
    Matrix grad2;
    grad2.multiply(Xnod,dNdX); // grad = Xnod * dNdX
    grad(row,1) = grad2(1,1);
    grad(row++,2) = grad2(1,2);
    vit += surf->getNoNodes(it)*surf->getNoFields(it);
  }

  return true;
}


bool SplineFields2Dmx::hessianFE(const FiniteElement& fe, Matrix3D& H) const
{
  return false;
}


bool SplineFields2Dmx::gradCoor (const Vec3& x, Matrix& grad) const
{
  // Not implemented yet
  return false;
}
