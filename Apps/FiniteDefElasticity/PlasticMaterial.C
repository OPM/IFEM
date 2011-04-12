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
	       const int& lfirst, const double* pMAT,
	       const int& nDF, const double& detF,
	       const double* Fn1, const double* Fn,
	       const double* be, const double& Epp, const double* Epl,
	       const double* Sig, double* Cst, int& ierr);
  //! \brief Interface to 3D elasto-plastic material routines (FORTRAN-77 code).
  void plas3d_(const int& ipsw, const int& iwr, const int& iter,
	       const int& lfirst, const int& nSig, const int& nTM,
	       const double* pMAT,
	       const double& detF, const double* Fn1, const double* Fn,
	       const double* be, const double& Epp, const double* Epl,
	       const double* Sig, double* Cst, int& ierr);
}
#ifndef INT_DEBUG
#define INT_DEBUG 0
#endif
#endif


/*!
  See the Fortran the routine Material/plas3d.f for interpretation of \a pMAT.
*/

PlasticMaterial::PlasticMaterial (const RealArray& p) : pMAT(p), iP1(0), iP2(0)
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
    Bmod = Emod / (3.0 - 6.0*nu);
    Smod = Emod / (2.0 + 2.0*nu);
  }
  else
  {
    Bmod = 0.0;
    Smod = Emod / 3.0;
  }

  pMAT[0] = Emod;
  pMAT[1] = nu;
  pMAT.insert(pMAT.begin()+4,Bmod);
  pMAT.insert(pMAT.begin()+5,Smod);
}


PlasticMaterial::~PlasticMaterial ()
{
  size_t i;
  for (i = 0; i < itgPoints.size(); i++)
    delete itgPoints[i];
  for (i = 0; i < resPoints.size(); i++)
    delete resPoints[i];
}


void PlasticMaterial::print (std::ostream& os) const
{
  std::cout <<"PlasticMaterial: pMAT =";
  for (size_t i = 0; i < pMAT.size(); i++)
    std::cout <<" "<< pMAT[i];
  std::cout << std::endl;
}


void PlasticMaterial::initIntegration (const TimeDomain& prm)
{
#if INT_DEBUG > 0
  std::cout <<"PlasticMaterial::initIntegration: "<< iP1 << std::endl;
#endif

  iP1 = 0;
  iAmIntegrating = true;

  if (prm.it > 0 || prm.first) return;

  int nUpdated = 0;
  for (size_t i = 0; i < itgPoints.size(); i++)
    if (itgPoints[i]->updateHistoryVars())
      nUpdated++;

#if INT_DEBUG > 0
  std::cout <<"PlasticMaterial::initIntegration: History updated "
	    << nUpdated << std::endl;
#endif
}


void PlasticMaterial::initResultPoints ()
{
#if INT_DEBUG > 0
  std::cout <<"PlasticMaterial::initResultPoints: "<< iP1 << std::endl;
#endif

  iP2 = 0;
  iAmIntegrating = false;

  int nUpdated = 0;
  for (size_t i = 0; i < resPoints.size(); i++)
    if (resPoints[i]->updateHistoryVars())
      nUpdated++;

#if INT_DEBUG > 0
  std::cout <<"PlasticMaterial::initResultPoints: History updated "
	    << nUpdated << std::endl;
#endif
}


bool PlasticMaterial::evaluate (Matrix& C, SymmTensor& sigma, double& U,
				const Vec3&, const Tensor& F,
				const SymmTensor& eps, char iop,
				const TimeDomain* prm) const
{
  if (iAmIntegrating)
  {
    if (!prm)
    {
      std::cerr <<" *** PlasticMaterial::evaluate: No time domain."<< std::endl;
      return false;
    }

    while (itgPoints.size() <= iP1)
      itgPoints.push_back(new PlasticPoint(this,F.dim()));

    if (prm->it == 0 && !prm->first)
      itgPoints[iP1]->Fp = F;

    if (!itgPoints[iP1++]->evaluate(C,sigma,F,*prm))
      return false;
  }
  else // Result evaluation
  {
    while (resPoints.size() <= iP2)
      resPoints.push_back(new PlasticPoint(this,F.dim()));

    if (!resPoints[iP2]->evaluate(C,sigma,F,TimeDomain(0,false)))
      return false;

    // Assume only one evaluation per increment; always update Fp
    resPoints[iP2++]->Fp = F;
  }

  if (iop > 1)
  {
    // Transform to 2nd Piola-Kirchhoff stresses,
    // via pull-back to reference configuration
    Tensor Fi(F);
    double J = Fi.inverse();
    sigma.transform(Fi); // sigma = F^-1 * sigma * F^-t
    sigma *= J;

    //TODO: When mixed formulation, we need the standard F here (not Fbar)
    //TODO: If invoked with iop=2, also pull-back the C-matrix
  }

  if (iop == 3)
  {
    if (iAmIntegrating)
      U = itgPoints[iP1-1]->energyIntegral(sigma,eps);
    else
      U = resPoints[iP2-1]->energyIntegral(sigma,eps);
  }

  return true;
}


PlasticMaterial::PlasticPoint::PlasticPoint (const PlasticMaterial* prm,
					     unsigned short int n)
  : pMAT(prm->pMAT), updated(false), Ep(0), Sp(0), Fp(n)
{
  // Initialize the history variables
  memset(HVc,0,sizeof(HVc));
  memset(HVp,0,sizeof(HVp));
  Fp = HVp[0] = HVp[1] = HVp[2] = 1.0;
  Up = 0.0;
}


bool PlasticMaterial::PlasticPoint::updateHistoryVars ()
{
  if (!updated) return false;

  // Update history variables with values of the new converged solution
  memcpy(HVp,HVc,sizeof(HVc));
  updated = false;

  return true;
}


bool PlasticMaterial::PlasticPoint::evaluate (Matrix& C,
					      SymmTensor& sigma,
					      const Tensor& F,
					      const TimeDomain& prm) const
{
  double J = F.det();
  if (J == 0.0)
  {
    std::cerr <<" *** PlasticMaterial::evaluate: "
	      <<" Singular/zero deformation gradient\n"<< F;
    return false;
  }

  // Restore history variables from the previous, converged configuration
  memcpy(const_cast<double*>(HVc),HVp,sizeof(HVc));

  size_t ndim = sigma.dim();
  size_t ncmp = ndim*(ndim+1)/2;
  C.resize(ncmp,ncmp);

  int ierr = -99;
#ifdef USE_FTNMAT
  // Invoke the FORTRAN routine for plasticity material models
#if INT_DEBUG > 0
  std::cout <<"PlasticMaterial::Fp =\n"<< Fp;
  std::cout <<"PlasticMaterial::HV =";
  for (int i = 0; i < 10; i++) std::cout <<" "<< HVc[i];
  std::cout << std::endl;
#endif
  if (ndim == 2)
    plas2d_(INT_DEBUG,6,prm.it,prm.first,&pMAT.front(),F.dim(),J,
	    F.ptr(),Fp.ptr(),HVc,HVc[6],HVc+7,sigma.ptr(),C.ptr(),ierr);
  else
    plas3d_(INT_DEBUG,6,prm.it,prm.first,6,6,&pMAT.front(),J,
	    F.ptr(),Fp.ptr(),HVc,HVc[6],HVc+7,sigma.ptr(),C.ptr(),ierr);
#if INT_DEBUG > 0
  std::cout <<"PlasticMaterial::sigma =\n"<< sigma;
  std::cout <<"PlasticMaterial::C ="<< C;
#endif
#else
  std::cerr <<" *** PlasticMaterial::evaluate: Not included."<< std::endl;
#endif

  const_cast<PlasticPoint*>(this)->updated = true;
  return ierr == 0;
}


double PlasticMaterial::PlasticPoint::energyIntegral (const SymmTensor& S,
						      const SymmTensor& E)
{
  Up += 0.5*((S+Sp)*(E-Ep));
  Sp.copy(S);
  Ep.copy(E);
  return Up;
}
