// $Id$
//==============================================================================
//!
//! \file LinearElasticity.C
//!
//! \date Nov 12 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for linear elasticity problems.
//!
//==============================================================================

#include "LinearElasticity.h"
#include "FiniteElement.h"
#include "MaterialBase.h"
#include "Tensor.h"
#include "Vec3Oper.h"


LinearElasticity::LinearElasticity (unsigned short int n) : Elasticity(n)
{
  // Only the current solution is needed
  primsol.resize(1);
}


bool LinearElasticity::evalInt (LocalIntegral*& elmInt, const FiniteElement& fe,
				const Vec3& X) const
{
  bool lHaveStrains = false;
  SymmTensor eps(nsd), sigma(nsd);

  if (eKm || eKg || iS)
  {
    // Compute the strain-displacement matrix B from dNdX
    // and evaluate the symmetric strain tensor if displacements are available
    if (!this->kinematics(fe.dNdX,eps,eps)) return false;
    if (!eps.isZero(1.0e-16)) lHaveStrains = true;

    // Evaluate the constitutive matrix and the stress tensor at this point
    double U;
    if (!material->evaluate(Cmat,sigma,U,X,eps,eps))
      return false;
  }

  if (eKm)
  {
    // Integrate the material stiffness matrix
    CB.multiply(Cmat,Bmat).multiply(fe.detJxW); // CB = C*B*|J|*w
    eKm->multiply(Bmat,CB,true,false,true); // EK += B^T * CB
  }

  if (eKg && lHaveStrains)
    // Integrate the geometric stiffness matrix
    this->formKG(*eKg,fe.dNdX,sigma,fe.detJxW);

  if (eM)
    // Integrate the mass matrix
    this->formMassMatrix(*eM,fe.N,X,fe.detJxW);

  if (iS && lHaveStrains)
  {
    // Integrate the internal forces
    sigma *= -fe.detJxW;
    if (!Bmat.multiply(sigma,*iS,true,true)) // ES -= B^T*sigma
      return false;
  }

  if (eS)
    // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(*eS,fe.N,X,fe.detJxW);

  return this->getIntegralResult(elmInt);
}


bool LinearElasticity::evalBou (LocalIntegral*& elmInt, const FiniteElement& fe,
				const Vec3& X, const Vec3& normal) const
{
  if (!tracFld)
  {
    std::cerr <<" *** LinearElasticity::evalBou: No tractions."<< std::endl;
    return false;
  }
  else if (!eS)
  {
    std::cerr <<" *** LinearElasticity::evalBou: No load vector."<< std::endl;
    return false;
  }

  // Evaluate the surface traction
  Vec3 T = (*tracFld)(X,normal);

  // Store the traction value for vizualization
  if (!T.isZero()) tracVal[X] = T;

  // Integrate the force vector
  Vector& ES = *eS;
  for (size_t a = 1; a <= fe.N.size(); a++)
    for (unsigned short int i = 1; i <= nsd; i++)
      ES(nsd*(a-1)+i) += T[i-1]*fe.N(a)*fe.detJxW;

  return this->getIntegralResult(elmInt);
}
