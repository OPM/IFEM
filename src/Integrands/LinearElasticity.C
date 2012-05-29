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
#include "ElmMats.h"
#include "Tensor.h"
#include "Vec3Oper.h"


LinearElasticity::LinearElasticity (unsigned short int n, bool axS)
  : Elasticity(n,axS)
{
}


bool LinearElasticity::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
				const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  bool lHaveStrains = false;
  SymmTensor eps(nsd,axiSymmetry), sigma(nsd,axiSymmetry);

  Matrix Bmat, Cmat;
  if (eKm || eKg || iS)
  {
    // Compute the strain-displacement matrix B from N, dNdX and r = X.x,
    // and evaluate the symmetric strain tensor if displacements are available
    if (!this->kinematics(elMat.vec.front(),fe.N,fe.dNdX,X.x,Bmat,eps,eps))
      return false;
    else if (!eps.isZero(1.0e-16))
      lHaveStrains = true;

    // Evaluate the constitutive matrix and the stress tensor at this point
    double U;
    if (!material->evaluate(Cmat,sigma,U,fe.iGP,X,eps,eps))
      return false;
  }

  // Axi-symmetric integration point volume; 2*pi*r*|J|*w
  const double detJW = axiSymmetry ? 2.0*M_PI*X.x*fe.detJxW : fe.detJxW;

  if (eKm)
  {
    // Integrate the material stiffness matrix
    Matrix CB;
    CB.multiply(Cmat,Bmat).multiply(detJW); // CB = C*B*|J|*w
    elMat.A[eKm-1].multiply(Bmat,CB,true,false,true); // EK += B^T * CB
  }

  if (eKg && lHaveStrains)
  {
    // Integrate the geometric stiffness matrix
    double r = axiSymmetry ? X.x + elMat.vec.front().dot(fe.N,0,nsd) : 0.0;
    this->formKG(elMat.A[eKg-1],fe.N,fe.dNdX,r,sigma,detJW);
  }

  if (eM)
    // Integrate the mass matrix
    this->formMassMatrix(elMat.A[eM-1],fe.N,X,detJW);

  if (iS && lHaveStrains)
  {
    // Integrate the internal forces
    sigma *= -detJW;
    if (!Bmat.multiply(sigma,elMat.b[iS-1],true,true)) // ES -= B^T*sigma
      return false;
  }

  if (eS)
    // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(elMat.b[eS-1],fe.N,X,detJW);

  return true;
}


bool LinearElasticity::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
				const Vec3& X, const Vec3& normal) const
{
  if (!tracFld && !fluxFld)
  {
    std::cerr <<" *** LinearElasticity::evalBou: No tractions."<< std::endl;
    return false;
  }
  else if (!eS)
  {
    std::cerr <<" *** LinearElasticity::evalBou: No load vector."<< std::endl;
    return false;
  }

  // Axi-symmetric integration point volume; 2*pi*r*|J|*w
  const double detJW = axiSymmetry ? 2.0*M_PI*X.x*fe.detJxW : fe.detJxW;

  // Evaluate the surface traction
  Vec3 T = this->getTraction(X,normal);

  // Store traction value for visualization
  if (fe.iGP < tracVal.size() && !T.isZero())
  {
    tracVal[fe.iGP].first = X;
    tracVal[fe.iGP].second += T;
  }

  // Integrate the force vector
  Vector& ES = static_cast<ElmMats&>(elmInt).b[eS-1];
  for (size_t a = 1; a <= fe.N.size(); a++)
    for (unsigned short int i = 1; i <= nsd; i++)
      ES(nsd*(a-1)+i) += T[i-1]*fe.N(a)*detJW;

  return true;
}
