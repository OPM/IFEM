// $Id$
//==============================================================================
//!
//! \file SplineFields3Dmx.C
//!
//! \date Oct 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for mixed spline-based finite element vector fields in 3D.
//!
//==============================================================================

#include "SplineFields3Dmx.h"
#include "ASMs3Dmx.h"
#include "FiniteElement.h"
#include "CoordinateMapping.h"
#include "Utilities.h"
#include "Vec3.h"

#include "GoTools/trivariate/SplineVolume.h"


SplineFields3Dmx::SplineFields3Dmx (const ASMs3Dmx* patch,
                                    const RealArray& v, char basis,
                                    const char* name)
  : Fields(name), svol(patch)
{
  bases = utl::getDigits(basis);
  nf = 3;
  auto vit = v.begin();
  size_t ofs = 0;
  for (int i = 1; i < *bases.begin(); ++i)
    ofs += patch->getNoNodes(i)*patch->getNoFields(i);
  vit += ofs;
  for (const auto& it : bases) {
    size_t nno = patch->getNoFields(it)*patch->getNoNodes(it);
    RealArray::const_iterator end = v.size() > nno+ofs ? vit+nno : v.end();
    std::copy(vit,end,std::back_inserter(values));
    vit += nno;
    ofs += nno;
    values.resize(ofs);
  }
}


bool SplineFields3Dmx::valueNode (size_t node, Vector& vals) const
{
  return false;
}


bool SplineFields3Dmx::valueFE (const FiniteElement& fe, Vector& vals) const
{
  if (!svol) return false;

  vals.resize(3);

  // Evaluate the basis functions at the given point
  auto vit = values.begin();
  auto rit = vals.begin();
  for (const auto& it : bases) {
    Go::SplineVolume* basis = svol->getBasis(it);
    Go::BasisPts spline;
#pragma omp critical
    basis->computeBasis(fe.u,fe.v,fe.w,spline);

    // Evaluate the solution field at the given point
    std::vector<int> ip;
    ASMs3D::scatterInd(basis->numCoefs(0),basis->numCoefs(1),
                       basis->numCoefs(2),
                       basis->order(0),basis->order(1),basis->order(2),
                       spline.left_idx,ip);

    Matrix Vnod;
    utl::gather(ip,1,Vector(&*vit,svol->getNoNodes(it)),Vnod);
    Vector val2;
    Vnod.multiply(spline.basisValues,val2); // vals = Vnod * basisValues
    *rit++ = val2.front();
    vit += svol->getNoNodes(it);
  }

  return true;
}


bool SplineFields3Dmx::valueCoor (const Vec3& x, Vector& vals) const
{
  // Not implemented yet
  return false;
}


bool SplineFields3Dmx::gradFE (const FiniteElement& fe, Matrix& grad) const
{
  if (!svol)  return false;

  // Evaluate the basis functions at the given point
  Go::BasisDerivs spline;
  const Go::SplineVolume* gvol = svol->getBasis(ASMmxBase::geoBasis);
#pragma omp critical
  gvol->computeBasis(fe.u,fe.v,fe.w,spline);

  const int uorder = gvol->order(0);
  const int vorder = gvol->order(1);
  const int worder = gvol->order(2);
  const size_t nen = uorder*vorder*worder;

  Matrix dNdu(nen,3), dNdX;
  for (size_t n = 1; n <= nen; n++)
  {
    dNdu(n,1) = spline.basisDerivs_u[n-1];
    dNdu(n,2) = spline.basisDerivs_v[n-1];
    dNdu(n,3) = spline.basisDerivs_w[n-1];
  }

  std::vector<int> ip;
  ASMs3D::scatterInd(gvol->numCoefs(0),gvol->numCoefs(1),gvol->numCoefs(2),
		     uorder,vorder,worder,spline.left_idx,ip);

  // Evaluate the Jacobian inverse
  Matrix Xnod, Jac;
  Vector Xctrl(&(*gvol->coefs_begin()),gvol->coefs_end()-gvol->coefs_begin());
  utl::gather(ip,gvol->dimension(),Xctrl,Xnod);
  utl::Jacobian(Jac,dNdX,Xnod,dNdu);

  // Evaluate the gradient of the solution field at the given point
  auto vit = values.begin();
  size_t row = 1;
  for (const auto& it : bases) {
    const Go::SplineVolume* basis = svol->getBasis(it);
#pragma omp critical
    basis->computeBasis(fe.u,fe.v,fe.w,spline);

    const size_t nbf = basis->order(0)*basis->order(1)*basis->order(2);
    dNdu.resize(nbf,3);
    for (size_t n = 1; n <= nbf; n++)
    {
      dNdu(n,1) = spline.basisDerivs_u[n-1];
      dNdu(n,2) = spline.basisDerivs_v[n-1];
      dNdu(n,3) = spline.basisDerivs_w[n-1];
    }
    dNdX.multiply(dNdu,Jac); // dNdX = dNdu * Jac

    ip.clear();
    ASMs3D::scatterInd(basis->numCoefs(0),basis->numCoefs(1),basis->numCoefs(2),
		       basis->order(0),basis->order(1),basis->order(2),
		       spline.left_idx,ip);

    utl::gather(ip,1,Vector(&*vit, svol->getNoNodes(it)*
                                   svol->getNoFields(it)),Xnod);
    Matrix grad2;
    grad2.multiply(Xnod,dNdX); // grad = Xnod * dNdX
    grad(row,1) = grad2(1,1);
    grad(row,2) = grad(1,2);
    grad(row++,3) = grad(1,3);
    vit += svol->getNoNodes(it)*svol->getNoFields(it);
  }

  return true;
}


bool SplineFields3Dmx::hessianFE(const FiniteElement& fe, Matrix3D& H) const
{
  return false;
}


bool SplineFields3Dmx::gradCoor (const Vec3& x, Matrix& grad) const
{
  // Not implemented yet
  return false;
}
