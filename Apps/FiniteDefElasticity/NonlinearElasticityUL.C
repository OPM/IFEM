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
#include "FiniteElement.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "TimeDomain.h"
#include "Tensor.h"
#include "Vec3Oper.h"


NonlinearElasticityUL::NonlinearElasticityUL (unsigned short int n, char lop)
  : Elasticity(n), loadOp(lop), plam(-999.9)
{
  // Only the current solution is needed
  primsol.resize(1);
}


void NonlinearElasticityUL::print (std::ostream& os) const
{
  std::cout <<"NonlinearElasticityUL: Updated Lagrangian formulation"
	    << std::endl;

  this->Elasticity::print(os);
}


void NonlinearElasticityUL::setMode (SIM::SolutionMode mode)
{
  if (!myMats) return;

  myMats->rhsOnly = false;
  eM = eKm = eKg = 0;
  iS = eS  = eV  = 0;

  switch (mode)
    {
    case SIM::STATIC:
      myMats->resize(1,1);
      eKm = &myMats->A[0];
      eKg = &myMats->A[0];
      break;

    case SIM::DYNAMIC:
      myMats->resize(2,1);
      eKm = &myMats->A[0];
      eKg = &myMats->A[0];
      eM  = &myMats->A[1];
      break;

    case SIM::RHS_ONLY:
      if (myMats->A.empty())
	myMats->resize(1,1);
      else
        myMats->b.resize(1);
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
  mySols.resize(1);
  eV = &mySols[0];
  tracVal.clear();
}


void NonlinearElasticityUL::initIntegration (const TimeDomain& prm)
{
  if (material)
    material->initIntegration(prm);
}


void NonlinearElasticityUL::initResultPoints (double lambda)
{
  if (material && lambda > plam)
  {
    material->initResultPoints();
    plam = lambda;
  }
}


bool NonlinearElasticityUL::evalInt (LocalIntegral*& elmInt,
				     const FiniteElement& fe,
				     const TimeDomain& prm,
				     const Vec3& X) const
{
  // Evaluate the deformation gradient, F, and the Green-Lagrange strains, E
  Tensor F(nsd);
  SymmTensor E(nsd);
  if (!this->kinematics(fe.dNdX,F,E))
    return false;

  bool lHaveStrains = !E.isZero(1.0e-16);
  if (lHaveStrains)
  {
    // Invert the deformation gradient ==> Fi
    Matrix Fi(nsd,nsd);
    Fi.fill(F.ptr());
    double J = Fi.inverse();
    if (J == 0.0) return false;

    // Scale with J=|F| since we are integrating over current configuration
    const_cast<double&>(fe.detJxW) *= J;

    if (eKm || iS)
    {
      // Push-forward the basis function gradients to current configuration
      dNdx.multiply(fe.dNdX,Fi); // dNdx = dNdX * F^-1
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
    if (!this->formBmatrix(fe.dNdX)) return false;

  // Evaluate the constitutive relation
  SymmTensor sigma(nsd);
  if (eKm || eKg || iS)
  {
    double U = 0.0;
    if (!material->evaluate(Cmat,sigma,U,X,F,E,(eKg || iS),&prm))
      return false;
  }

  if (eKm)
  {
    // Integrate the material stiffness matrix
    CB.multiply(Cmat,Bmat).multiply(fe.detJxW); // CB = C*B*|J|*w
    eKm->multiply(Bmat,CB,true,false,true);     // EK += B^T * CB
  }

  if (eKg && lHaveStrains)
    // Integrate the geometric stiffness matrix
    this->formKG(*eKg,dNdx,sigma,fe.detJxW);

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


bool NonlinearElasticityUL::evalBou (LocalIntegral*& elmInt,
				     const FiniteElement& fe,
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
    if (!this->formDefGradient(fe.dNdX,F)) return false;

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
      const_cast<double&>(fe.detJxW) *= F.det();
  }

  // Integrate the force vector
  Vector& ES = *eS;
  for (size_t a = 1; a <= fe.N.size(); a++)
    for (i = 1; i <= nsd; i++)
      ES(nsd*(a-1)+i) += T[i-1]*fe.N(a)*fe.detJxW;

  return this->getIntegralResult(elmInt);
}


bool NonlinearElasticityUL::formDefGradient (const Matrix& dNdX,
					     Tensor& F) const
{
  static SymmTensor dummy(0);
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
  if (!dUdX.multiplyMat(*eV,dNdX)) // dUdX = Grad{u} = eV*dNdX
    return false;

  unsigned short int i, j, k;
  if (dUdX.rows() < F.dim())
    F.zero();

  // Cannot use operator= here, in case F is of higher dimension than dUdX
  for (i = 1; i <= dUdX.rows(); i++)
    for (j = 1; j <= dUdX.cols(); j++)
      F(i,j) = dUdX(i,j);

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
  return new ElasticityNormUL(const_cast<NonlinearElasticityUL&>(*this));
}


void ElasticityNormUL::initIntegration (const TimeDomain& prm)
{
  this->ElasticityNorm::initIntegration(prm);
  iP = 0;
}


bool ElasticityNormUL::evalInt (LocalIntegral*& elmInt,
				const FiniteElement& fe,
				const TimeDomain& prm,
				const Vec3& X) const
{
  NonlinearElasticityUL& ulp = static_cast<NonlinearElasticityUL&>(myProblem);

  ElmNorm& pnorm = NormBase::getElmNormBuffer(elmInt);

  // Evaluate the deformation gradient, F, and the Green-Lagrange strains, E
  Tensor F(fe.dNdX.cols());
  SymmTensor E(fe.dNdX.cols());
  if (!ulp.kinematics(fe.dNdX,F,E))
    return false;

  // Compute the strain energy density, U(E) = Int_E (Sig:Eps) dEps
  double U = 0.0;
  SymmTensor S(E.dim());
  if (!ulp.material->evaluate(ulp.Cmat,S,U,X,F,E,3,&prm))
    return false;

  // Integrate the energy norm a(u^h,u^h) = Int_Omega0 U(E) dV0
  pnorm[0] += U*fe.detJxW;
  // Integrate the L2-norm ||Sig|| = Int_Omega0 Sig:Sig dV0
  pnorm[2] += S.L2norm(false)*fe.detJxW;
  // Integrate the von Mises stress norm
  pnorm[3] += S.vonMises(false)*fe.detJxW;

  return true;
}


bool ElasticityNormUL::evalBou (LocalIntegral*& elmInt,
				const FiniteElement& fe,
				const Vec3& X, const Vec3& normal) const
{
  NonlinearElasticityUL& ulp = static_cast<NonlinearElasticityUL&>(myProblem);
  if (!ulp.haveLoads()) return true;

  ElmNorm& pnorm = NormBase::getElmNormBuffer(elmInt);

  // Evaluate the current surface traction
  Vec3 t = ulp.getTraction(X,normal);
  // Evaluate the current displacement field
  Vec3 u = ulp.evalSol(fe.N);

  // Integrate the external energy (path integral)
  if (iP < Ux.size())
  {
    // New load step, update integration point values
    Ux[iP] += 0.5*(t+tp[iP])*(u-up[iP]);
    tp[iP] = t;
    up[iP] = u;
  }
  else
  {
    // This is the first load step at this integration point
    Ux.push_back(0.5*t*u);
    tp.push_back(t);
    up.push_back(u);
  }

  pnorm[1] += Ux[iP++]*fe.detJxW;
  return true;
}
