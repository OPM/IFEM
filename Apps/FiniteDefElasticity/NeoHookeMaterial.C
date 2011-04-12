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
	       const int& nDF, const double& detF, const double* F,
	       const double* pMAT,
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
}


NeoHookeMaterial::NeoHookeMaterial (double E, double v, double d, int ver)
  : LinIsotropic(E,v,d,false)
{
  mTYP = 10;
  mVER = ver;
  pmat[0] = Emod;
  pmat[1] = nu;
  if (nu > 0.5)
  {
    mTYP = -10; // Assume Lame' parameters (kappa and mu are specified)
    double kappa = Emod;
    double mu = nu;
    // Calculate Young's modulus and Poisson's ratio
    Emod = 9.0*kappa*mu/(3.0*kappa + mu);
    nu   = (1.5*kappa - mu)/(3.0*kappa + mu);
  }
}


void NeoHookeMaterial::print (std::ostream& os) const
{
  this->LinIsotropic::print(os);
  std::cout <<"NeoHookeMaterial: mVER = "<< mVER << std::endl;
}


bool NeoHookeMaterial::evaluate (Matrix& C, SymmTensor& sigma, double& U,
				 const Vec3&, const Tensor& F,
				 const SymmTensor& eps, char iop,
				 const TimeDomain*) const
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
    cons2d_(INT_DEBUG,0,6,0,mTYP,mVER,F.dim(),J,F.ptr(),pmat,
	    U,sigma.ptr(),C.ptr(),ierr);
  else
    cons3d_(INT_DEBUG,0,6,0,mTYP,mVER,0,6,0,J,F.ptr(),pmat,&U,
	    U,sigma.ptr(),C.ptr(),ierr);
#else
  std::cerr <<" *** NeoHookeMaterial::evaluate: Not included."<< std::endl;
#endif

  if (iop > 1)
  {
    // Transform to 2nd Piola-Kirchhoff stresses,
    // via pull-back to reference configuration
    Tensor Fi(F);
    Fi.inverse();
    sigma.transform(Fi); // sigma = F^-1 * sigma * F^-t
    sigma *= J;

    //TODO: If invoked with iop=2, also pull-back the C-matrix
    //TODO: When mixed formulation, we need the standard F, not Fbar here
  }

  if (iop == 3 && U == 0.0)
    U = sigma.innerProd(eps); //TODO: Replace this by proper path integral

#if INT_DEBUG > 0
  if (iop > 0)
    std::cout <<"NeoHookeMaterial::sigma =\n"<< sigma;
  std::cout <<"NeoHookeMaterial::C ="<< C;
#endif

  return ierr == 0;
}
