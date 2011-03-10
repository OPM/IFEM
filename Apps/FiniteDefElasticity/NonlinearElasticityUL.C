// $Id$
//==============================================================================
//!
//! \file NonlinearElasticityUL.C
//!
//! \date Sep 21 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for nonlinear elasticity problems.
//!
//==============================================================================

#include "NonlinearElasticityUL.h"
#include "MaterialBase.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Tensor.h"
#include "Vec3Oper.h"


NonlinearElasticityUL::NonlinearElasticityUL (unsigned short int n, char lop)
  : Elasticity(n)
{
  loadOp = lop;

  // Only the current solution is needed
  primsol.resize(1);
}


void NonlinearElasticityUL::print (std::ostream& os) const
{
  material->print(os);
  std::cout <<"NonlinearElasticityUL: Updated Lagrangian formulation"
	    << std::endl;
}


void NonlinearElasticityUL::setMode (SIM::SolutionMode mode)
{
  if (!myMats) return;

  size_t nvec = 1 + primsol.size();
  myMats->rhsOnly = false;
  eM = eKm = eKg = 0;
  iS = eS  = eV  = 0;

  switch (mode)
    {
    case SIM::STATIC:
      myMats->resize(1,nvec);
      eKm = &myMats->A[0];
      eKg = &myMats->A[0];
      break;

    case SIM::DYNAMIC:
      myMats->resize(2,nvec);
      eKm = &myMats->A[0];
      eKg = &myMats->A[0];
      eM  = &myMats->A[1];
      break;

    case SIM::RHS_ONLY:
      if (myMats->A.empty())
	myMats->resize(1,nvec);
      else
        myMats->b.resize(nvec);
      eKm = &myMats->A[0];
      eKg = &myMats->A[0];
      myMats->rhsOnly = true;
      break;

    default:
      this->Elasticity::setMode(mode);
      return;
    }

  // We always need the force+displacement vectors in nonlinear simulations
  iS = &myMats->b[0];
  eS = &myMats->b[0];
  eV = &myMats->b[1];
  tracVal.clear();
}


bool NonlinearElasticityUL::evalInt (LocalIntegral*& elmInt, double detJW,
				     const Vector& N, const Matrix& dNdX,
				     const Vec3& X) const
{
  // Evaluate the deformation gradient, F, and the Green-Lagrange strains, E
  Tensor F(nsd);
  SymmTensor E(nsd);
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
#ifdef INT_DEBUG
      std::cout <<"NonlinearElasticityUL::dNdx ="<< dNdx;
      std::cout <<"NonlinearElasticityUL::B ="<< Bmat;
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
    if (!material->evaluate(Cmat,sigma,U,X,F,E, (eKg || iS) && lHaveStrains))
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


bool NonlinearElasticityUL::evalBou (LocalIntegral*& elmInt, double detJW,
				     const Vector& N, const Matrix& dNdX,
				     const Vec3& X, const Vec3& normal) const
{
  if (!tracFld)
  {
    std::cerr <<" *** NonlinearElasticityUL::evalBou: No tractions."
	      << std::endl;
    return false;
  }
  else if (!eS)
  {
    std::cerr <<" *** NonlinearElasticityUL::evalBou: No load vector."
	      << std::endl;
    return false;
  }

  // Evaluate the surface traction
  Vec3 T = (*tracFld)(X,normal);

  // Store the traction value for vizualization
  if (!T.isZero()) tracVal[X] = T;

  unsigned short int i, j;
  if (loadOp == 1)
  {
    // Compute the deformation gradient, F
    Tensor F(nsd);
    if (!this->formDefGradient(dNdX,F)) return false;

    // Check for with-rotated pressure load
    if (tracFld->isNormalPressure())
    {
      // Compute its inverse and determinant, J
      double J = F.inverse();
      if (J == 0.0) return false;

      // Pull-back the normal traction to the initial configuration.
      // See equation (3.4.5) on page 102 in T. Belytschko et. al (2000):
      //    p*n*dS ==> J*p*n0*F^-1*dS0
      Vec3 t = J*T; T = 0.0;
      for (i = 1; i <= nsd; i++)
	for (j = 1; j <= nsd; j++)
	  T[i-1] += t[j-1]*F(j,i);
    }
    else
      // Scale with J=|F| since we are integrating over current configuration
      detJW *= F.det();
  }

  // Integrate the force vector
  Vector& ES = *eS;
  for (size_t a = 1; a <= N.size(); a++)
    for (i = 1; i <= nsd; i++)
      ES(nsd*(a-1)+i) += T[i-1]*N(a)*detJW;

  return this->getIntegralResult(elmInt);
}


bool NonlinearElasticityUL::formDefGradient (const Matrix& dNdX,
					     Tensor& F) const
{
  SymmTensor dummy(0);
  return this->kinematics(dNdX,F,dummy);
}


bool NonlinearElasticityUL::kinematics (const Matrix& dNdX,
					Tensor& F, SymmTensor& E) const
{
  if (!eV || eV->empty())
  {
    // Initial state, unit deformation gradient and zero strains
    F = 1.0;
    E.zero();
    return true;
  }

  const size_t nenod = dNdX.rows();
  if (eV->size() != nenod*nsd || dNdX.cols() < nsd)
  {
    std::cerr <<" *** NonlinearElasticityUL::kinematics: Invalid dimension,"
	      <<" dNdX("<< nenod <<","<< dNdX.cols() <<")"<< std::endl;
    return false;
  }

  // Compute the deformation gradient, [F] = [I] + [dudX] = [I] + [dNdX]*[u],
  // and the Green-Lagrange strains, E_ij = 0.5(F_ij+F_ji+F_ki*F_kj).

  // Notice that the matrix multiplication method used here treats the
  // element displacement vector, *eV, as a matrix whose number of columns
  // equals the number of rows in the matrix dNdX.
  Matrix dUdX;
  if (dUdX.multiplyMat(*eV,dNdX)) // dUdX = Grad{u} = eV*dNdX
    F = dUdX;
  else
    return false;

  unsigned short int i, j, k;

  // Now form the Green-Lagrange strain tensor.
  // Note that for the shear terms (i/=j) we actually compute 2*E_ij
  // to be consistent with the engineering strain style constitutive matrix.
  for (i = 1; i <= E.dim(); i++)
    for (j = 1; j <= i; j++)
    {
      double Eij = F(i,j) + F(j,i);
      for (k = 1; k <= nsd; k++)
	Eij += F(k,i)*F(k,j);
      E(i,j) = i == j ? 0.5*Eij : Eij;
    }

  // Add the unit tensor to F to form the deformation gradient
  F += 1.0;

#ifdef INT_DEBUG
  std::cout <<"NonlinearElasticityUL::eV ="<< *eV;
  std::cout <<"NonlinearElasticityUL::F =\n"<< F;
#endif

  return true;
}


NormBase* NonlinearElasticityUL::getNormIntegrand (AnaSol*) const
{
  return new ElasticityNormUL(*const_cast<NonlinearElasticityUL*>(this));
}


bool ElasticityNormUL::evalInt (LocalIntegral*& elmInt, double detJW,
				const Vector&, const Matrix& dNdX,
				const Vec3& X) const
{
  ElmNorm& pnorm = ElasticityNorm::getElmNormBuffer(elmInt);

  NonlinearElasticityUL* ulp = static_cast<NonlinearElasticityUL*>(&problem);

  // Evaluate the deformation gradient, F, and the Green-Lagrange strains, E
  Tensor F(dNdX.cols());
  SymmTensor E(dNdX.cols());
  if (!ulp->kinematics(dNdX,F,E))
    return false;

  // Evaluate the 2nd Piola Kirchhoff stresses, S
  Matrix C;
  SymmTensor S(E.dim());
  double U = 0.0;
  if (!ulp->material->evaluate(C,S,U,X,F,E,2))
    return false;

  if (U == 0.0)
  {
    // Integrate the energy norm a(u^h,u^h) = Int_Omega0 (S:E) dV0
    const RealArray& sigma = S;
    const RealArray& epsil = E;
    for (size_t i = 0; i < sigma.size(); i++)
      U += sigma[i]*epsil[i];
  }

  pnorm[0] += U*detJW;
  return true;
}
