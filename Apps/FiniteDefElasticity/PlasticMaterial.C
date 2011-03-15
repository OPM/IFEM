// $Id$
//==============================================================================
//!
//! \file PlasticMaterial.C
//!
//! \date Mar 16 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Plasticity material models.
//!
//==============================================================================

#include "PlasticMaterial.h"
#include "TimeDomain.h"

#ifdef USE_FTNMAT
extern "C" {
  //! \brief Interface to 2D plastic material routines (FORTRAN-77 code).
  void plas2d_(const int& ipsw, const int& iwr, const int& iter,
	       const int& lfirst, const int& ntm, const double* pMAT,
	       const double& detF, const double* Fn1, const double* Fn,
	       const double* be, const double& Epp, const double* Epl,
	       const double* Sig, double* Cst, int& ierr);
  //! \brief Interface to 3D plastic material routines (FORTRAN-77 code).
  void plas3d_(const int& ipsw, const int& iwr, const int& iter,
	       const int& lfirst, const int& ntm, const double* pMAT,
	       const double& detF, const double* Fn1, const double* Fn,
	       const double* be, const double& Epp, const double* Epl,
	       const double* Sig, double* Cst, int& ierr);
}
#ifndef INT_DEBUG
#define INT_DEBUG 0
#endif
#endif


PlasticMaterial::PlasticMaterial (const PlasticPrm* prm, unsigned short int n)
  : pMAT(prm->pMAT), Fp(n)
{
  memset(HVc,0,sizeof(HVc));
  memset(HVp,0,sizeof(HVc));
}


void PlasticPrm::print (std::ostream& os) const
{
  std::cout <<"PlasticMaterial: pMAT =";
  for (size_t i = 0; i < pMAT.size(); i++)
    std::cout <<" "<< pMAT[i];
  if (rho > 0.0)
    std::cout <<"\nPlasticMaterial: rho ="<< rho;
  std::cout << std::endl;
}


bool PlasticMaterial::evaluate (Matrix& C, SymmTensor& sigma, double&,
				const Vec3&, const Tensor& F, const SymmTensor&,
				char, const TimeDomain* prm) const
{
  double J = F.det();
  if (J == 0.0)
  {
    std::cerr <<" *** PlasticMaterial::evaluate: "
	      <<" Singular/zero deformation gradient\n"<< F;
    return false;
  }
  else if (!prm)
  {
    std::cerr <<" *** PlasticMaterial::evaluate: No time domain."<< std::endl;
    return false;
  }

  size_t ndim = F.dim();
  size_t ncmp = ndim*(ndim+1)/2;
  C.resize(ncmp,ncmp);
  int ierr = -99;
#ifdef USE_FTNMAT
  // Invoke the FORTRAN routine for plasticity material models.
  memcpy(const_cast<double*>(HVc),HVp,sizeof(HVc));
  if (ndim == 2)
    plas2d_(INT_DEBUG,6,prm->it,prm->first,4,&pMAT.front(),J,F.ptr(),Fp.ptr(),
	    HVc,HVc[6],HVc+7,sigma.ptr(),C.ptr(),ierr);
  else
    plas3d_(INT_DEBUG,6,prm->it,prm->first,6,&pMAT.front(),J,F.ptr(),Fp.ptr(),
	    HVc,HVc[6],HVc+7,sigma.ptr(),C.ptr(),ierr);
#if INT_DEBUG > 0
  std::cout <<"PlasticMaterial::sigma =\n"<< sigma;
  std::cout <<"PlasticMaterial::C ="<< C;
#endif
#else
  std::cerr <<" *** PlasticMaterial::evaluate: Not included."<< std::endl;
#endif

  return ierr == 0;
}


void PlasticMaterial::updateHistoryVars ()
{
  memcpy(HVp,HVc,sizeof(HVc));
}
