// $Id$
//==============================================================================
//!
//! \file TensorFunction.C
//!
//! \date May 10 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Spatial tensor-valued functions.
//!
//==============================================================================

#include "TensorFunction.h"


utl::matrix3d<Real> TensorFunc::gradient(const Vec3& X) const
{
  const size_t nsd = sqrt(ncmp);
  utl::matrix3d<Real> result(nsd, nsd, nsd);
  result.fill(this->evalGradient(X).data());
  return result;
}


Tensor TensorFunc::timeDerivative(const Vec3& X) const
{
  const size_t nsd = sqrt(ncmp);
  Tensor result(nsd);
  result = this->evalTimeDerivative(X);
  return result;
}


size_t STensorFunc::index(size_t nsd, size_t i, size_t j) const
{
  if (i == j)
    return i-1; // diagonal term
  else if (nsd == 2)
    return ncmp-1; // off-diagonal term (2D)

  if (i == j+1 || i+2 == j)
    std::swap(i,j);

  return i+2; // upper triangular term (3D)
}


utl::matrix3d<Real> STensorFunc::gradient(const Vec3& X) const
{
  const size_t nsd = ncmp > 5 ? 3 : (ncmp > 2 ? 2 : 1);
  utl::matrix3d<Real> result(nsd,nsd,nsd);
  const std::vector<Real> temp = this->evalGradient(X);

  for (size_t d = 1; d <= nsd; ++d)
    for (size_t i = 1; i <= nsd; ++i)
      for (size_t j = 1; j <= nsd; ++j)
        result(i,j,d) = temp[index(nsd,i,j) + (d-1)*ncmp];

  return result;
}


SymmTensor STensorFunc::timeDerivative(const Vec3& X) const
{
  const size_t nsd = ncmp > 5 ? 3 : (ncmp > 2 ? 2 : 1);
  SymmTensor result(nsd, ncmp == 4);
  result = this->evalTimeDerivative(X);
  return result;
}
