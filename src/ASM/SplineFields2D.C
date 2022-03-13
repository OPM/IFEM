// $Id$
//==============================================================================
//!
//! \file SplineFields2D.C
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for spline-based finite element vector fields in 2D.
//!
//==============================================================================

#include "GoTools/geometry/SplineSurface.h"

#include "SplineFields2D.h"
#include "SplineField.h"

#include "ASMs2D.h"
#include "ItgPoint.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Vec3.h"


SplineFields2D::SplineFields2D (const ASMs2D* patch,
                                const RealArray& v, char nbasis,
                                int nnf, const char* name)
  : Fields(name), basis(patch->getBasis(nbasis)), surf(patch->getSurface())
{
  const int n1 = basis->numCoefs_u();
  const int n2 = basis->numCoefs_v();
  nno = n1*n2;

  const int p1 = basis->order_u();
  const int p2 = basis->order_v();
  nelm = (n1-p1+1)*(n2-p2+1);

  nsd = patch->getNoSpaceDim();

  size_t ofs = 0;
  for (char i = 1; i < nbasis; ++i)
    ofs += patch->getNoNodes(i)*patch->getNoFields(i);
  auto vit = v.begin()+ofs;

  // Ensure the values array has compatible length, pad with zeros if necessary
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
}


SplineFields2D::SplineFields2D (const Go::SplineSurface* srf,
                                const RealArray& v, int cmp, const char* name)
  : Fields(name), basis(srf), surf(srf)
{
  nsd = 2; // That is, this constructor can not be used for shells (nsd=3)
  nf = cmp;
  values = v;
}


bool SplineFields2D::valueNode (size_t node, Vector& vals) const
{
  if (node < 1 || node > nno) return false;

  vals.resize(nf);
  vals.fill(values.ptr()+(node-1)*nf);
  return true;
}


bool SplineFields2D::valueFE (const ItgPoint& x, Vector& vals) const
{
  if (!basis) return false;

  // Evaluate the basis functions at the given point
  Go::BasisPtsSf spline;
#pragma omp critical
  basis->computeBasis(x.u,x.v,spline);

  // Evaluate the solution field at the given point
  std::vector<int> ip;
  ASMs2D::scatterInd(basis->numCoefs_u(),basis->numCoefs_v(),
		     basis->order_u(),basis->order_v(),
		     spline.left_idx,ip);

  Matrix Vnod;
  utl::gather(ip,nf,values,Vnod);
  Vnod.multiply(spline.basisValues,vals); // vals = Vnod * basisValues

  return true;
}


bool SplineFields2D::valueCoor (const Vec4& x, Vector& vals) const
{
  if (x.u)
    return this->valueFE(ItgPoint(x.u[0],x.u[1]),vals);

  // Use with caution, very slow!
  Go::Point pt(x.x,x.y,x.z), clopt(3);
  double clo_u, clo_v, dist;
#pragma omp critical
  surf->closestPoint(pt, clo_u, clo_v, clopt, dist, 1.0e-5);

  return this->valueFE(ItgPoint(clo_u,clo_v),vals);
}


bool SplineFields2D::gradFE (const ItgPoint& x, Matrix& grad) const
{
  if (!basis) return false;
  if (!surf)  return false;

  // Evaluate the basis functions at the given point
  std::vector<int> ip;
  Matrix Xnod, Jac, dNdX;
  if (!SplineField::evalMapping(*surf,nsd,x,ip,Xnod,Jac,dNdX))
    return false;

  if (surf != basis)
    if (!SplineField::evalBasis(*basis,x,ip,Xnod,Jac,dNdX))
      return false;

  Matrix Vnod;
  utl::gather(ip,nf,values,Vnod);
  return !grad.multiply(Vnod,dNdX).empty(); // grad = Vnod * dNdX
}


bool SplineFields2D::hessianFE (const ItgPoint& x, Matrix3D& H) const
{
  if (!basis) return false;
  if (!surf)  return false;

  std::vector<int> ip;
  Matrix Xnod, Jac, dNdX;
  Matrix3D d2NdX2, Hess;
  if (!SplineField::evalMapping(*surf,nsd,x,ip,Xnod,Jac,dNdX,&d2NdX2,&Hess))
    return false;

  if (surf != basis)
    if (!SplineField::evalBasis(*basis,x,ip,Xnod,Jac,dNdX,&d2NdX2,&Hess))
      return false;

  Matrix Vnod;
  utl::gather(ip,nf,values,Vnod);

  return H.multiply(Vnod,d2NdX2);
}
