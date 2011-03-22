// $Id$
//==============================================================================
//!
//! \file PlasticMaterial.C
//!
//! \date Mar 16 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Elasto-plastic material models.
//!
//==============================================================================

#include "PlasticMaterial.h"
#include "TimeDomain.h"

#ifdef USE_FTNMAT
extern "C" {
  //! \brief Interface to 2D elasto-plastic material routines (FORTRAN-77 code).
  void plas2d_(const int& ipsw, const int& iwr, const int& iter,
	       const int& lfirst, const int& ntm, const double* pMAT,
	       const double& detF, const double* Fn1, const double* Fn,
	       const double* be, const double& Epp, const double* Epl,
	       const double* Sig, double* Cst, int& ierr);
  //! \brief Interface to 3D elasto-plastic material routines (FORTRAN-77 code).
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


/*!
  See the Fortran routine Material/plas3d.f for interpretation of \a pMAT.
*/

PlasticPrm::PlasticPrm (const RealArray& p) : pMAT(p)
{
  if (pMAT.size() < 11) pMAT.resize(11,0.0);

  double Emod = pMAT[0];
  double nu   = pMAT[1];
  double Bmod, Smod;
  if (nu > 0.5)
  {
    // Calculate E and nu from the Bulk and Shear moduli
    Bmod = Emod;
    Smod = nu;
    Emod = 9.0*Bmod*Smod/(3.0*Bmod + Smod);
    nu   = (1.5*Bmod - Smod)/(3.0*Bmod + Smod);
  }
  else if (nu < 0.5)
  {
    // Calculate the Bulk and Shear moduli from E and nu
    Smod = Emod / (2.0 + 2.0*nu);
    Bmod = Emod / (3.0 - 6.0*nu);
  }
  else
  {
    Smod = Emod / 3.0;
    Bmod = 0.0;
  }

  pMAT[0] = Emod;
  pMAT[1] = nu;
  pMAT.insert(pMAT.begin()+4,Bmod);
  pMAT.insert(pMAT.begin()+5,Smod);
}


void PlasticPrm::print (std::ostream& os) const
{
  std::cout <<"PlasticMaterial: pMAT =";
  for (size_t i = 0; i < pMAT.size(); i++)
    std::cout <<" "<< pMAT[i];
  std::cout << std::endl;
}


PlasticMaterial::PlasticMaterial (const PlasticPrm* prm, unsigned short int n)
  : pMAT(prm->pMAT), Fp(n)
{
  memset(HVc,0,sizeof(HVc));
  memset(HVp,0,sizeof(HVc));
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
