// $Id$
//==============================================================================
//!
//! \file NeoHookeMaterial.C
//!
//! \date Mar 08 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Neo-Hookean hyperelastic material model.
//!
//==============================================================================

#include "NeoHookeMaterial.h"
#include "Tensor.h"

#ifdef USE_FTNMAT
extern "C" {
  //! \brief Interface to 2D nonlinear material routines (FORTRAN-77 code).
  void cons2d_(const int& ipsw, const int& iter, const int& iwr,
	       const int& lfirst, const int& mTYP, const int& mVER,
	       const int& nSig, const int& nDF, const double& detF,
	       const double* F, const double* pMAT,
	       double& Engy, const double* Sig, double* Cst, int& ierr);
  //! \brief Interface to 3D nonlinear material routines (FORTRAN-77 code).
  void cons3d_(const int& ipsw, const int& iter, const int& iwr,
	       const int& lfirst, const int& mTYP, const int& mVER,
	       const int& nHV, const int& nSig, const int& nTM,
	       const double& detF, const double* F,
	       const double* pMAT, double* HV,
	       double& Engy, const double* Sig, double* Cst, int& ierr);
}
#ifndef INT_DEBUG
#define INT_DEBUG 0
#endif
#endif


NeoHookeMaterial::NeoHookeMaterial () : LinIsotropic(false)
{
  mTYP = 10;
  mVER = 0;
  pmat[0] = Emod;
  pmat[1] = nu;
  pmat[2] = pmat[3] = pmat[4] = 0.0;
  pmat[5] = 0.5 * Emod / (1.0 + nu);
  if (nu == 0.5) // incompressible material
    pmat[4] = 0.0;
  else
  {
    pmat[4] = Emod / (3.0 - 6.0*nu);
    if (mVER == 1)
      pmat[4] -= pmat[5]*2.0/3.0;
  }
  pmat[6] = 1;
}


NeoHookeMaterial::NeoHookeMaterial (double E, double v, double d, int ver)
  : LinIsotropic(E,v,d,false)
{
  mTYP = 10;
  mVER = ver / 10;
  pmat[0] = Emod;
  pmat[1] = nu;
  pmat[2] = pmat[3] = pmat[4] = 0.0;
  if (nu > 0.5)
  {
    mTYP = -10; // Assume Lame' parameters (kappa and mu are specified)
    double kappa = Emod;
    double mu = nu;
    // Calculate Young's modulus and Poisson's ratio
    Emod = 9.0*kappa*mu/(3.0*kappa + mu);
    nu   = (1.5*kappa - mu)/(3.0*kappa + mu);
  }
  else
  {
    pmat[5] = 0.5 * Emod / (1.0 + nu);
    if (nu == 0.5) // incompressible material
      pmat[4] = 0.0;
    else
    {
      pmat[4] = Emod / (3.0 - 6.0*nu);
      if (mVER == 1)
	pmat[4] -= pmat[5]*2.0/3.0;
    }
  }
  pmat[6] = ver % 10;
}


void NeoHookeMaterial::print (std::ostream& os) const
{
  this->LinIsotropic::print(os);
  std::cout <<"NeoHookeMaterial: mVER = "<< mVER << std::endl;
}


bool NeoHookeMaterial::evaluate (Matrix& C, SymmTensor& sigma, double& U,
				 const Vec3&, const Tensor& F,
				 const SymmTensor& eps, char iop,
				 const TimeDomain*, const Tensor* Fpf) const
{
  double J = F.det();
  if (J == 0.0)
  {
    std::cerr <<" *** NeoHookeMaterial::evaluate: "
	      <<" Singular/zero deformation gradient\n"<< F;
    return false;
  }

  size_t ndim = sigma.dim();
  size_t ncmp = ndim*(ndim+1)/2;
  C.resize(ncmp,ncmp);
  int ierr = -99;
#ifdef USE_FTNMAT
  // Invoke the FORTRAN routine for Neo-Hookean hyperelastic material models.
  if (ndim == 2)
    cons2d_(INT_DEBUG,0,6,0,mTYP,mVER,sigma.size(),F.dim(),J,F.ptr(),pmat,
	    U,sigma.ptr(),C.ptr(),ierr);
  else
    cons3d_(INT_DEBUG,0,6,0,mTYP,mVER,0,sigma.size(),0,J,F.ptr(),pmat,&U,
	    U,sigma.ptr(),C.ptr(),ierr);
#else
  std::cerr <<" *** NeoHookeMaterial::evaluate: Not included."<< std::endl;
#endif

  if (iop == 2 || (iop == 3 && U == 0.0))
  {
    // Transform to 2nd Piola-Kirchhoff stresses,
    // via pull-back to reference configuration
    Tensor Fi(Fpf ? *Fpf : F);
    J = Fi.inverse();
    if (iop == 2)
    {
      sigma.transform(Fi); // sigma = F^-1 * sigma * F^-t
      sigma *= J;
      //TODO: Also pull-back the C-matrix (Total Lagrange formulation)
    }
    else
    {
      SymmTensor S(sigma); // Make a copy of sigma since it should be Cauchy stress when iop=3
      U = S.transform(Fi).innerProd(eps)*J; //TODO: Replace this by proper path integral
    }
  }

#if INT_DEBUG > 0
  if (iop > 0)
    std::cout <<"NeoHookeMaterial::sigma =\n"<< sigma;
  std::cout <<"NeoHookeMaterial::C ="<< C;
#endif

  return ierr == 0;
}
