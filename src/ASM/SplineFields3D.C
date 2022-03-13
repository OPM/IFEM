// $Id$
//==============================================================================
//!
//! \file SplineFields3D.C
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for spline-based finite element vector fields in 3D.
//!
//==============================================================================

#include "GoTools/trivariate/SplineVolume.h"

#include "SplineFields3D.h"
#include "SplineField.h"

#include "ASMs3D.h"
#include "ItgPoint.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Vec3.h"


SplineFields3D::SplineFields3D (const ASMs3D* patch,
                                const RealArray& v, char nbasis,
                                int nnf, const char* name)
  : Fields(name), basis(patch->getBasis(nbasis)), vol(patch->getVolume())
{
  const int n1 = basis->numCoefs(0);
  const int n2 = basis->numCoefs(1);
  const int n3 = basis->numCoefs(2);
  nno = n1*n2*n3;

  const int p1 = basis->order(0);
  const int p2 = basis->order(1);
  const int p3 = basis->order(2);
  nelm = (n1-p1+1)*(n2-p2+1)*(n3-p3+1);

  nsd = patch->getNoSpaceDim();

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


SplineFields3D::SplineFields3D (const Go::SplineVolume* svol,
                                const RealArray& v, int cmp, const char* name)
  : Fields(name), basis(svol), vol(svol)
{
  nsd = 3;
  nf = cmp;
  values = v;
}


bool SplineFields3D::valueNode (size_t node, Vector& vals) const
{
  if (node < 1 || node > nno) return false;

  vals.resize(nf);
  vals.fill(values.ptr()+(node-1)*nf);
  return true;
}


bool SplineFields3D::valueFE (const ItgPoint& x, Vector& vals) const
{
  if (!basis) return false;

  // Evaluate the basis functions at the given point
  Go::BasisPts spline;
#pragma omp critical
  basis->computeBasis(x.u,x.v,x.w,spline);

  // Evaluate the solution field at the given point
  std::vector<int> ip;
  ASMs3D::scatterInd(basis->numCoefs(0),basis->numCoefs(1),basis->numCoefs(2),
                     basis->order(0),basis->order(1),basis->order(2),
                     spline.left_idx,ip);

  Matrix Vnod;
  utl::gather(ip,nf,values,Vnod);
  Vnod.multiply(spline.basisValues,vals); // vals = Vnod * basisValues

  return true;
}


bool SplineFields3D::valueCoor (const Vec4& x, Vector& vals) const
{
  if (x.u)
    return this->valueFE(ItgPoint(x.u[0],x.u[1],x.u[2]),vals);

  // Use with caution, very slow!
  Go::Point pt(x.x,x.y,x.z), clopt(3);
  double clo_u, clo_v, clo_w, dist;
#pragma omp critical
  vol->closestPoint(pt, clo_u, clo_v, clo_w, clopt, dist, 1.0e-5);

  return this->valueFE(ItgPoint(clo_u,clo_v,clo_w),vals);
}


bool SplineFields3D::gradFE (const ItgPoint& x, Matrix& grad) const
{
  if (!basis) return false;
  if (!vol)   return false;

  // Evaluate the basis functions at the given point
  IntVec ip;
  Matrix Xnod, Jac, dNdX;
  if (!SplineField::evalMapping(*vol,x,ip,Xnod,Jac,dNdX))
    return false;

  // Evaluate the gradient of the solution field at the given point
  if (basis != vol)
    if (!SplineField::evalBasis(*basis,x,ip,Xnod,Jac,dNdX))
      return false;

  Matrix Vnod;
  utl::gather(ip,nf,values,Vnod);
  return !grad.multiply(Vnod,dNdX).empty(); // grad = Vnod * dNdX
}


bool SplineFields3D::hessianFE (const ItgPoint& x, Matrix3D& H) const
{
  if (!basis) return false;
  if (!vol)  return false;

  // Evaluate the basis functions at the given point
  IntVec ip;
  Matrix Xnod, Jac, dNdX;
  Matrix3D d2NdX2, Hess;
  if (!SplineField::evalMapping(*vol,x,ip,Xnod,Jac,dNdX,&d2NdX2,&Hess))
    return false;

  if (vol != basis)
    if (!SplineField::evalBasis(*basis,x,ip,Xnod,Jac,dNdX,&d2NdX2,&Hess))
      return false;

  Matrix Vnod;
  utl::gather(ip,nf,values,Vnod);

  return H.multiply(Vnod,d2NdX2);
}
