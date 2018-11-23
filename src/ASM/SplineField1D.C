// $Id$
//==============================================================================
//!
//! \file SplineField1D.C
//!
//! \date Nov 23 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for spline-based finite element scalar field in 1D.
//!
//==============================================================================

#include "GoTools/geometry/SplineCurve.h"

#include "SplineField1D.h"
#include "ASMs1D.h"
#include "FiniteElement.h"
#include "CoordinateMapping.h"
#include "Utilities.h"
#include "Vec3.h"
#include <numeric>
#include <cassert>


SplineField1D::SplineField1D (const ASMs1D* patch,
                              const RealArray& v, const char* name)
  : FieldBase(name), curv(patch->getCurve())
{
  values = v;
  nno = values.size();
}


double SplineField1D::valueNode (size_t node) const
{
  return node > 0 && node <= nno ? values(node) : 0.0;
}


double SplineField1D::valueFE (const FiniteElement& fe) const
{
  if (!curv) return false;

  // Evaluate the basis functions at the given point
  RealArray basisVal, basisDerivs;
#pragma omp critical
  const_cast<Go::SplineCurve*>(curv)->computeBasis(fe.u,basisVal,basisDerivs);

  // Evaluate the solution field at the given point
  int first, last;
  curv->basis().coefsAffectingParam(fe.u, first, last);
  IntVec ip(last-first+1);;
  std::iota(ip.begin(),ip.end(),first);
  Vector Vnod;
  utl::gather(ip,1,values,Vnod);
  return Vnod.dot(basisVal);
}


double SplineField1D::valueCoor (const Vec4& x) const
{
  assert(0);
  return 0.0;
}


bool SplineField1D::valueGrid (RealArray& val, const int* npe) const
{
  return false;
}


bool SplineField1D::gradFE (const FiniteElement& fe, Vector& grad) const
{
  if (!curv)  return false;

  // Evaluate the basis functions at the given point
  RealArray basisVal, basisDerivs;
#pragma omp critical
  curv->computeBasis(fe.u,basisVal,basisDerivs);

  // Evaluate the solution field at the given point
  int first, last;
  curv->basis().coefsAffectingParam(fe.u, first, last);
  IntVec ip(last-first+1);;
  std::iota(ip.begin(),ip.end(),first);
  Vector Vnod;
  utl::gather(ip,1,values,Vnod);
  grad(1) = Vnod.dot(basisDerivs);
  return true;
}


bool SplineField1D::hessianFE (const FiniteElement& fe, Matrix& H) const
{
  if (!curv)  return false;

  // Evaluate the basis functions at the given point
  RealArray basisVal, basisDerivs, basisDerivs2;
#pragma omp critical
  const_cast<Go::SplineCurve*>(curv)->computeBasis(fe.u,basisVal,basisDerivs, basisDerivs2);

  // Evaluate the solution field at the given point
  int first, last;
  curv->basis().coefsAffectingParam(fe.u, first, last);
  IntVec ip(last-first+1);;
  std::iota(ip.begin(),ip.end(),first);
  Vector Vnod;
  utl::gather(ip,1,values,Vnod);
  H.resize(1,1,1);
  H(1,1) = Vnod.dot(basisDerivs2);
  return true;
}
