// $Id$
//==============================================================================
//!
//! \file SplineFields1D.C
//!
//! \date Nov 23 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for spline-based finite element vector fields in 1D.
//!
//==============================================================================

#include "GoTools/geometry/SplineCurve.h"

#include "SplineFields1D.h"
#include "FiniteElement.h"
#include "CoordinateMapping.h"
#include "Utilities.h"
#include "Vec3.h"
#include <cassert>
#include <numeric>


SplineFields1D::SplineFields1D (const Go::SplineCurve* crv,
                                const RealArray& v, int cmp, const char* name)
  : Fields(name), curv(crv)
{
  nf = cmp;
  values = v;
}


bool SplineFields1D::valueNode (size_t node, Vector& vals) const
{
  if (node < 1 || node > nno) return false;

  vals.resize(nf);
  vals.fill(values.ptr()+(node-1)*nf);
  return true;
}


bool SplineFields1D::valueFE (const FiniteElement& fe, Vector& vals) const
{
  if (!curv) return false;

  // Evaluate the basis functions at the given point
  RealArray basisVal, basisDerivs;
#pragma omp critical
  const_cast<Go::SplineCurve*>(curv)->computeBasis(fe.u,basisVal,basisDerivs);

  // Evaluate the solution field at the given point
  int first, last;
  curv->basis().coefsAffectingParam(fe.u, first, last);
  std::vector<int> ip((last-first+1)*nf);;
  std::iota(ip.begin(),ip.end(),first*nf);
  Matrix Vnod;
  utl::gather(ip,1,values,Vnod);
  Vnod.multiply(basisVal,vals); // vals = Vnod * basisValues

  return true;
}


bool SplineFields1D::valueCoor (const Vec4& x, Vector& vals) const
{
  assert(0);
  return false;
}


bool SplineFields1D::gradFE (const FiniteElement& fe, Matrix& grad) const
{
  if (!curv)  return false;

  // Evaluate the basis functions at the given point
  RealArray basisVal, basisDerivs;
#pragma omp critical
  curv->computeBasis(fe.u,basisVal,basisDerivs);

  // Evaluate the solution field at the given point
  int first, last;
  curv->basis().coefsAffectingParam(fe.u, first, last);
  std::vector<int> ip((last-first+1)*nf);;
  std::iota(ip.begin(),ip.end(),first*nf);
  Matrix Vnod(last-first+1,nf);
  utl::gather(ip,1,values,Vnod);
  for (size_t i = 1; i <= nf; ++i)
    grad(i,1) = Vnod.getColumn(i).dot(basisDerivs);
  return true;
}


bool SplineFields1D::hessianFE (const FiniteElement& fe, Matrix3D& H) const
{
  if (!curv)  return false;

  // Evaluate the basis functions at the given point
  RealArray basisVal, basisDerivs;
#pragma omp critical
  curv->computeBasis(fe.u,basisVal,basisDerivs);

  // Evaluate the solution field at the given point
  int first, last;
  curv->basis().coefsAffectingParam(fe.u, first, last);
  std::vector<int> ip((last-first+1)*nf);;
  std::iota(ip.begin(),ip.end(),first*nf);
  Matrix Vnod(last-first+1,nf);
  utl::gather(ip,1,values,Vnod);
  H.resize(nf,1,1);
  for (size_t i = 1; i <= nf; ++i)
    H(i,1,1) = Vnod.getColumn(i).dot(basisDerivs);

  return true;
}
