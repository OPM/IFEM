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
#include "ItgPoint.h"
#include "CoordinateMapping.h"
#include "Utilities.h"
#include <numeric>


SplineFields1D::SplineFields1D (const Go::SplineCurve* crv,
                                const RealArray& v, int ncmp, const char* name)
  : Fields(name), curv(crv)
{
  nsd = 1; // That is, this constructor can not be used for curved beams (nsd>1)
  nf = ncmp;
  values = v;
}


bool SplineFields1D::valueFE (const ItgPoint& x, Vector& vals) const
{
  if (!curv) return false;

  // Find nodal indices of element containing given point
  int first, last;
  curv->basis().coefsAffectingParam(x.u,first,last);
  if (first == last)
  {
    // We are at an interpolatory point (with C^0 continuity)
    vals.resize(nf);
    vals.fill(values.ptr()+nf*first,nf);
    return true;
  }

  // Find nodal field values around this point
  std::vector<int> ip(last-first+1);
  Matrix Vnod;
  std::iota(ip.begin(),ip.end(),first);
  utl::gather(ip,nf,values,Vnod);

  // Evaluate the basis functions at the given point
  RealArray basisVal, basisDerivs;
#pragma omp critical
  curv->computeBasis(x.u,basisVal,basisDerivs);

  // Evaluate the field at the given point.
  // Notice we don't just do a matrix-vector multiplication here,
  // to account for that the basisVal array may be longer that ip.
  // This seems to happen when we have repeated knots in the affected element.
  vals.resize(nf,true);
  size_t inod = 0;
  for (double N : basisVal)
    if (N > 0.0 && inod < Vnod.cols())
      vals.add(Vnod.getColumn(++inod),N);

  return true;
}


bool SplineFields1D::gradFE (const ItgPoint& x, Matrix& grad) const
{
  if (!curv) return false;

  // Evaluate the basis functions at the given point
  RealArray basisVal, basisDerivs;
#pragma omp critical
  curv->computeBasis(x.u,basisVal,basisDerivs);
  Matrix Jac, dNdX, dNdu(basisDerivs.size(),1);
  dNdu.fillColumn(1,basisDerivs);

  // Find nodal indices of element containing the given point
  int first, last;
  curv->basis().coefsAffectingParam(x.u,first,last);
  std::vector<int> ip(last-first+1);
  std::iota(ip.begin(),ip.end(),first);

  // Find nodal coordinates of element containing the given point
  Matrix Velm, Xelm(nsd,ip.size());
  RealArray::const_iterator cit = curv->coefs_begin();
  for (int inod = first; inod <= last; inod++)
  {
    int iip = inod*curv->dimension();
    for (size_t i = 0; i < nsd; i++)
      Xelm(i+1,inod-first+1) = *(cit+(iip+i));
  }

  // Evaluate the field gradient at the given point
  utl::gather(ip,nf,values,Velm);
  if (!utl::Jacobian(Jac,dNdX,Xelm,dNdu))
    return false; // Singular Jacobian

  return !grad.multiply(Velm,dNdX).empty(); // grad = Xelm * dNdX
}


bool SplineFields1D::hessianFE (const ItgPoint& x, Matrix3D& H) const
{
  if (!curv)  return false;

  // Evaluate the basis functions at the given point
  RealArray basisVal, basisDerivs1, basisDerivs2;
  // TODO: It does not make sense that we here need to const_cast
  // when it is not needed in the first derivatives version used above.
  // Can this be fixed in GoTools?
  Go::SplineCurve* ccrv = const_cast<Go::SplineCurve*>(curv);
#pragma omp critical
  ccrv->computeBasis(x.u,basisVal,basisDerivs1,basisDerivs2);
  Matrix Jac, dNdX, dNdu(basisDerivs1.size(),1);
  dNdu.fillColumn(1,basisDerivs1);
  Matrix3D Hess, d2NdX2, d2Ndu2(basisDerivs2.size(),1,1);
  d2Ndu2.fillColumn(1,1,basisDerivs2);

  // Find nodal indices of element containing the given point
  int first, last;
  curv->basis().coefsAffectingParam(x.u,first,last);
  std::vector<int> ip(last-first+1);
  std::iota(ip.begin(),ip.end(),first);

  // Find nodal coordinates of element containing the given point
  Matrix Velm, Xelm(nsd,ip.size());
  RealArray::const_iterator cit = curv->coefs_begin();
  for (int inod = first; inod <= last; inod++)
  {
    int iip = inod*curv->dimension();
    for (size_t i = 0; i < nsd; i++)
      Xelm(i+1,inod-first+1) = *(cit+(iip+i));
  }

  // Evaluate the field hessian at the given point
  utl::gather(ip,nf,values,Velm);
  if (!utl::Jacobian(Jac,dNdX,Xelm,dNdu))
    return false;
  if (!utl::Hessian(Hess,d2NdX2,Jac,Xelm,d2Ndu2,dNdX))
    return false;
  return H.multiply(Velm,d2NdX2); // H = Xelm * d2NdX2
}
