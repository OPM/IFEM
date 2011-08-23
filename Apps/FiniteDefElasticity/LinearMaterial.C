// $Id$
//==============================================================================
//!
//! \file LinearMaterial.C
//!
//! \date Mar 08 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief General linear elastic material with push-forward transformations.
//!
//==============================================================================

#include "LinearMaterial.h"
#include "Tensor.h"


double LinearMaterial::getMassDensity (const Vec3& X) const
{
  return material->getMassDensity(X);
}


bool LinearMaterial::evaluate (Matrix& C, SymmTensor& sigma, double& U,
			       const Vec3& X, const Tensor& Fbar,
			       const SymmTensor& eps, char iop,
			       const TimeDomain* prm, const Tensor* Fpf) const
{
  // Evaluate the constitutive matrix and the stress tensor at this point
  if (!material->evaluate(C,sigma,U,X,Fbar,eps,iop,prm))
    return false;
  else if (iop == 2)
    return true;

  const Tensor& F = Fpf ? *Fpf : Fbar;
  double J = F.det();
  if (J == 0.0)
  {
    std::cerr <<" *** LinearMaterial::evaluate: "
	      <<" Singular/zero deformation gradient\n"<< Fbar;
    return false;
  }

  if (iop > 0)
  {
    // Push-forward the stress tensor to current configuration
    // sigma = 1/J * F * sigma * F^t
    sigma.transform(F);
    sigma *= 1.0/J;
#ifdef INT_DEBUG
    std::cout <<"LinearMaterial::sigma =\n"<< sigma;
#endif
    if (iop > 1) return true;
  }

  // Push-forward the constitutive matrix to current configuration

  size_t i, j, k, l;
  size_t n = sigma.dim();
  size_t m = C.rows();
  Matrix T(m,m), Ctmp;

  for (i = 1; i <= n; i++)
  {
    for (j = 1; j <= n; j++)
      T(i,j) = F(j,i)*F(j,i);
    for (j = n+1, k = 1; j <= m; j++, k++)
      T(i,j) = 0.5*(F(k,i)*F(k%3+1,i) + F(k%3+1,i)*F(k,i));
  }

  for (i = n+1, k = 1; i <= m; i++, k++)
  {
    for (j = 1; j <= n; j++)
      T(i,j) = F(j,k)*F(j,k%3+1) + F(j,k%3+1)*F(j,k);
    for (j = n+1, l = 1; j <= m; j++, l++)
      T(i,j) = 0.5*(F(l,k)*F(l%3+1,k%3+1) + F(l%3+1,k)*F(l,k%3+1) +
		    F(l,k%3+1)*F(l%3+1,k) + F(l%3+1,k%3+1)*F(l,k));
  }

  // C = 1/J * T^t * C * T
  C.multiply(T,Ctmp.multiply(C,T),true);
  C *= 1.0/J;
#ifdef INT_DEBUG
  std::cout <<"LinearMaterial::C ="<< C;
#endif

  return true;
}
