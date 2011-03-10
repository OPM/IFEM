// $Id$
//==============================================================================
//!
//! \file NeoHookeElasticity.C
//!
//! \date Jul 1 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for Neo-Hooke elasticity problems.
//!
//==============================================================================

#include "NeoHookeElasticity.h"


void NeoHookeElasticity::print (std::ostream& os) const
{
  this->NonlinearElasticity::print(os);
  std::cout <<"NeoHookeElasticity: Neo-Hooke material formulation"<< std::endl;
}


void NeoHookeElasticity::setMaterial (double Emod, double Poiss, double Density)
{
  //this->Elasticity::setMaterial(Emod,Poiss,Density);

  // Define the Lame elasticity parameters
  lambda = Emod*Poiss / ((1.0+Poiss)*(1.0-Poiss-Poiss));
  mu     = lambda * (0.5/Poiss-1.0);
  kappa  = lambda + (mu+mu)/3.0;
}


bool NeoHookeElasticity::formStressTensor (const Matrix& dNdX, const Vec3&,
					   SymmTensor& S) const
{
  if (!eV || eV->empty())
  {
    // Initial state (zero stresses)
    S.zero();
    return true;
  }

  // Evaluate the kinematic quantities, F, E, Ci, and J at this point
  Tensor _F(nsd);
  if (!this->kinematics(dNdX,_F,E))
    return false;

  // Set up the 2nd Piola-Kirchhoff stress tensor, S
  this->formStressTensor(S);

  return true;
}


bool NeoHookeElasticity::formTangent (Matrix& Ctan, SymmTensor& S,
				      const Vec3& X) const
{
  if (!this->formCmatrix(Ctan,X))
    return false;
  else if (eV && !eV->empty())
    this->formStressTensor(S);

  return true;
}


bool NeoHookeElasticity::kinematics (const Matrix& dNdX,
				     Tensor& _F, SymmTensor& _E) const
{
  // Form the deformation gradient, F
  if (!this->NonlinearElasticity::kinematics(dNdX,_F,_E))
    return false;
  else
    J = _F.det();

  // Form the right Cauchy-Green tensor, C = F^T*F
  C.rightCauchyGreen(_F);

  // Invert the right Cauchy-Green tensor, Ci = C^-1
  Ci = C;
  return Ci.inverse() > 0.0;
}


void NeoHookeElasticity::formStressTensor (SymmTensor& S) const
{
  // Set up the 2nd Piola-Kirchhoff stress tensor.
  // See equation (3.119) in P. Wriggers book:
  // S = mu*(1-Ci) + 0.5*lambda*(J^2-1)*Ci

  S  =  0.5*lambda*(J*J-1.0)*Ci;
  S += -mu*(Ci-1.0);
}


bool NeoHookeElasticity::formCmatrix (Matrix& Ctan, const Vec3&, bool) const
{
  // Set up the incremental consitutive tensor
  size_t nst = nsd*(nsd+1)/2;
  Ctan.resize(nst,nst,true);
  SymmTensor4 D(Ctan,nsd); // fourth-order material tensor

  double a = lambda*J*J;
  double b = mu - 0.5*(a-lambda);

  unsigned short int i, j, k, l;
  for (i = 1; i <= nsd; i++)
    for (j = i; j <= nsd; j++)
      for (k = 1; k <= nsd; k++)
        for (l = k; l <= nsd; l++)
          // See equation (3.268) in P. Wriggers book
          D(i,j,k,l) = a*Ci(i,j)*Ci(k,l) + b*(Ci(i,k)*Ci(j,l) + Ci(i,l)*C(j,k));

  return true;
}


void NeoHookeElasticityIV::print (std::ostream& os) const
{
  this->NonlinearElasticity::print(os);
  std::cout <<"NeoHookeElasticity: Neo-Hooke material formulation"
	    <<" with isochoric/volumetric decomposition"<< std::endl;
}


bool NeoHookeElasticityIV::formTangent (Matrix& Ctan, SymmTensor& S,
					const Vec3& X) const
{
  // Set up the stress tensor first because Ctan also depends on Siso
  if (eV && !eV->empty())
    this->formStressTensor(S);
  else
  {
    // Initial state (zero stresses)
    S.zero();
    Siso.zero();
  }

  return this->formCmatrix(Ctan,X);
}


void NeoHookeElasticityIV::formStressTensor (SymmTensor& S) const
{
  // Set up the 2nd Piola-Kirchhoff stress tensor.
  // See equation (3.127) or (3.274) in P. Wriggers book:
  // S = mu*J^(-2/3)*(1-tr(C)/3)*Ci

  S  = -(C.trace()/3.0)*Ci;
  S +=  1.0;
  S *=  mu*pow(J,-2.0/3.0);
  Siso = S; // Isochoric part
  S += 0.5*kappa*(J*J-1.0)*Ci; // Volumetric part
}


bool NeoHookeElasticityIV::formCmatrix (Matrix& Ctan, const Vec3&, bool) const
{
  // Set up the incremental consitutive tensor
  size_t nst = nsd*(nsd+1)/2;
  Ctan.resize(nst,nst,true);
  SymmTensor4 D(Ctan,nsd); // fourth-order material tensor

  double t = mu*pow(J,-2.0/3.0)*C.trace();
  double a = kappa*J*J + (t+t)/9.0;
  double b = 0.5*kappa*(1.0-J*J) + t/3.0;
  double c = (t+t)/3.0; // Yuri
  //double c = 2.0/3.0; // Wriggers

  unsigned short int i, j, k, l;
  for (i = 1; i <= nsd; i++)
    for (j = i; j <= nsd; j++)
      for (k = 1; k <= nsd; k++)
	for (l = k; l <= nsd; l++)
	{
	  // See equation (3.275) in P. Wriggers book
	//D(i,j,k,l) = a*Ci(i,j)*Ci(k,l) + b*(Ci(i,k)*Ci(j,l) + Ci(i,l)*C(j,k))
	//  	     - c*(Ci(k,l)*Siso(i,j) + Ci(i,j)*Siso(k,l));
	  // See equation (81) In Yuri's Comput. Mech. paper
	  D(i,j,k,l) = a*Ci(i,j)*Ci(k,l) + b*(Ci(i,k)*Ci(j,l) + Ci(i,l)*C(j,k));
	  if (i == j) D(i,j,k,l) += c*Ci(k,l);
	  if (k == l) D(i,j,k,l) += c*Ci(i,j);
	}

  return true;
}
