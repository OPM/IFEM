// $Id$
//==============================================================================
//!
//! \file PlasticityUL.C
//!
//! \date Mar 16 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for finite deformation plasticity problems.
//!
//==============================================================================

#include "PlasticityUL.h"
#include "PlasticMaterial.h"
#include "ElmMats.h"
#include "TimeDomain.h"


PlasticityUL::PlasticityUL (unsigned short int n, char lop)
  : NonlinearElasticityUL(n,lop)
{
  iP = 0;

  // Both the current and previous solutions are needed
  primsol.resize(2);
}


PlasticityUL::~PlasticityUL ()
{
  for (size_t i = 0; i < pBuf.size(); i++)
    delete pBuf[i];
}


void PlasticityUL::print (std::ostream& os) const
{
  material->print(os);
  std::cout <<"PlasticityUL: Updated Lagrangian formulation"<< std::endl;
}


void PlasticityUL::initIntegration (const TimeDomain& prm)
{
  iP = 0;

  if (prm.it == 0 && !prm.first)
    for (size_t i = 0; i < pBuf.size(); i++)
      pBuf[i]->updateHistoryVars();
}


bool PlasticityUL::evalInt (LocalIntegral*& elmInt,
			    const TimeDomain& prm, double detJW,
			    const Vector& N, const Matrix& dNdX,
			    const Vec3& X) const
{
  while (pBuf.size() <= iP)
    pBuf.push_back(new PlasticMaterial(static_cast<PlasticPrm*>(material),nsd));
  PlasticMaterial* ptData = pBuf[iP++];

  if (prm.it == 0)
  {
    // Evaluate the deformation gradient, Fp, at previous configuration
    const_cast<PlasticityUL*>(this)->eV = &myMats->b[2];
    if (!this->formDefGradient(dNdX,ptData->defGrad()))
      return false;
  }

  // Evaluate the deformation gradient, F, and the Green-Lagrange strains, E
  Tensor F(nsd);
  SymmTensor E(nsd);
  const_cast<PlasticityUL*>(this)->eV = &myMats->b[1];
  if (!this->kinematics(dNdX,F,E))
    return false;

  bool lHaveStrains = !E.isZero();
  if (lHaveStrains)
  {
    // Invert the deformation gradient ==> Fi
    Matrix Fi(nsd,nsd);
    Fi.fill(F.ptr());
    double J = Fi.inverse();
    if (J == 0.0) return false;

    // Scale with J=|F| since we are integrating over current configuration
    detJW *= J;

    if (eKm || iS)
    {
      // Push-forward the basis function gradients to current configuration
      dNdx.multiply(dNdX,Fi); // dNdx = dNdX * F^-1
      // Compute the small-deformation strain-displacement matrix B from dNdx
      if (!this->formBmatrix(dNdx)) return false;
#if INT_DEBUG > 0
      std::cout <<"PlasticityUL::dNdx ="<< dNdx;
      std::cout <<"PlasticityUL::B ="<< Bmat;
#endif
    }
  }
  else if (eKm || iS)
    // Initial state, no deformation yet
    if (!this->formBmatrix(dNdX)) return false;

  // Evaluate the constitutive relation
  SymmTensor sigma(nsd);
  if (eKm || eKg || iS)
  {
    double U = 0.0;
    if (!ptData->evaluate(Cmat,sigma,U,X,F,E,true,&prm))
      return false;
  }

  if (eKm)
  {
    // Integrate the material stiffness matrix
    CB.multiply(Cmat,Bmat).multiply(detJW); // CB = C*B*|J|*w
    eKm->multiply(Bmat,CB,true,false,true); // EK += B^T * CB
  }

  if (eKg && lHaveStrains)
    // Integrate the geometric stiffness matrix
    this->formKG(*eKg,dNdx,sigma,detJW);

  if (eM)
    // Integrate the mass matrix
    this->formMassMatrix(*eM,N,X,detJW);

  if (iS && lHaveStrains)
  {
    // Integrate the internal forces
    sigma *= -detJW;
    if (!Bmat.multiply(sigma,*iS,true,true)) // ES -= B^T*sigma
      return false;
  }

  if (eS)
    // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(*eS,N,X,detJW);

  return this->getIntegralResult(elmInt);
}
