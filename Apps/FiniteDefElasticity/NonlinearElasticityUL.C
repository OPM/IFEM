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
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Tensor.h"
#include "Vec3Oper.h"

#ifdef USE_FTNMAT
extern "C" {
  //! \brief Interface to 2D nonlinear material routines (FORTRAN-77 code).
  void cons2d_(const int& mTYP, const int& mVER, const double& detF,
	       const double* F, const double* pmat, double& Engy,
	       const double* Sig, double* Cst,
	       const int& ipsw, const int& iwr, int& ierr);
  //! \brief Interface to 3D nonlinear material routines (FORTRAN-77 code).
  void cons3d_(const int& mTYP, const int& mVER, const double& detF,
	       const double* F, const double* pmat, double& Engy,
	       const double* Sig, double* Cst,
	       const int& ipsw, const int& iwr, int& ierr);
}
#ifndef INT_DEBUG
#define INT_DEBUG 0
#endif
#endif


NonlinearElasticityUL::NonlinearElasticityUL (unsigned short int n,
					      bool ps, int mver, char lop)
  : Elasticity(n,ps)
{
  mVER = -1;
#ifdef USE_FTNMAT
  if (n == 3 || (n == 2 && !ps)) mVER = mver;
#endif
  loadOp = lop;

  // Only the current solution is needed
  primsol.resize(1);
}


void NonlinearElasticityUL::print (std::ostream& os) const
{
  this->Elasticity::print(os);
  std::cout <<"NonlinearElasticityUL: Updated Lagrangian formulation";
  if (mVER >= 0) std::cout <<", mVER="<< mVER;
  if (mVER >= 1) std::cout <<" (Neo-Hooke)";
  std::cout << std::endl;
}


void NonlinearElasticityUL::setMaterial (double Young, double Poiss,
					 double Density)
{
  mTYP = 10;
  pmat[0] = Young;
  pmat[1] = Poiss;
  if (Poiss > 0.5)
  {
    mTYP *= -1; // Assume Lame' parameters (kappa and mu are specified)
    double kappa = Young;
    double mu = Poiss;
    // Calculate Young's modulus and Poisson's ratio
    Young = 9.0*kappa*mu/(3.0*kappa + mu);
    Poiss = (1.5*kappa - mu)/(3.0*kappa + mu);
  }
  this->Elasticity::setMaterial(Young,Poiss,Density);
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
      myMats->resize(1,2);
      eKm = &myMats->A[0];
      eKg = &myMats->A[0];
      break;

    case SIM::DYNAMIC:
      myMats->resize(2,2);
      eKm = &myMats->A[0];
      eKg = &myMats->A[0];
      eM  = &myMats->A[1];
      break;

    case SIM::RHS_ONLY:
      if (myMats->A.empty())
	myMats->resize(1,2);
      else
        myMats->b.resize(2);
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
  elmInt = myMats;

  // Evaluate the deformation gradient, F, and the Green-Lagrange strains, E
  SymmTensor E(nsd);
  if (!this->kinematics(dNdX,E))
    return false;

  bool lHaveStrains = !E.isZero();
  if (lHaveStrains)
  {
    // Invert the deformation gradient ==> Fi
    Matrix Fi(F);
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
  double U = 0.0;
  if (eKm || eKg || iS)
    if (!this->constitutive(Cmat,sigma,U,E,X, (eKg || iS) && lHaveStrains))
      return false;

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

  return true;
}


bool NonlinearElasticityUL::evalBou (LocalIntegral*& elmInt, double detJW,
				     const Vector& N, const Matrix& dNdX,
				     const Vec3& X, const Vec3& normal) const
{
  elmInt = myMats;
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
    if (!this->formDefGradient(dNdX)) return false;

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

  return true;
}


bool NonlinearElasticityUL::formDefGradient (const Matrix& dNdX) const
{
  static SymmTensor dummy(0);
  return this->kinematics(dNdX,dummy);
}


bool NonlinearElasticityUL::kinematics (const Matrix& dNdX, SymmTensor& E) const
{
  if (!eV || eV->empty())
  {
    // Initial state, unit deformation gradient and zero strains
    F.diag(1.0,nsd);
    E.zero();
    return true;
  }

  const size_t nenod = dNdX.rows();
  const size_t nstrc = nsd*(nsd+1)/2;
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
  if (!F.multiplyMat(*eV,dNdX)) // F = Grad{u} = eV*dNdX
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
  for (i = 1; i <= nsd; i++)
    F(i,i) += 1.0;

#if INT_DEBUG > 0
  std::cout <<"NonlinearElasticityUL::eV ="<< *eV;
  std::cout <<"NonlinearElasticityUL::F ="<< F;
#endif

  return true;
}


bool NonlinearElasticityUL::constitutive (Matrix& C, SymmTensor& sigma,
					  double& U, const SymmTensor& eps,
					  const Vec3& X,
					  char calcStress) const
{
  double J = F.det();
  if (J == 0.0)
  {
    std::cerr <<" *** NonlinearElasticityUL::constitutive: "
	      <<" Singular/zero deformation gradient"<< F;
    return false;
  }

  size_t ndim = F.rows();
#ifdef USE_FTNMAT
  if (mVER >= 0)
  {
    // Invoke the FORTRAN routine for constitutive relations.
    // This also includes the Neo-Hookean hyperelastic material models.
    size_t ncmp = ndim*(ndim+1)/2;
    C.resize(ncmp,ncmp);
    int ipsw = INT_DEBUG > 1 ? 9 : 0;
    int ierr = 0;
    if (ndim == 2)
      cons2d_(mTYP,mVER,J,F.ptr(),pmat,U,sigma.ptr(),C.ptr(),ipsw,6,ierr);
    else
      cons3d_(mTYP,mVER,J,F.ptr(),pmat,U,sigma.ptr(),C.ptr(),ipsw,6,ierr);
    if (calcStress == 2)
    {
      // Transform to 2nd Piola-Kirchhoff stresses,
      // via pull-back to reference configuration
      F.inverse();
      Tensor Fi(ndim);
      sigma.transform(Fi=F); // sigma = F^-1 * sigma * F^-t
      sigma *= J;
    }
#if INT_DEBUG > 0
    if (calcStress)
      std::cout <<"NonlinearElasticityUL::sigma =\n"<< sigma;
    std::cout <<"NonlinearElasticityUL::C ="<< C;
#endif
    return ierr == 0;
  }
#endif

  // We are only doing pure isotrophic linear elastic material in this method.
  // Evaluate the constitutive matrix at this point, but first reset nsd to ndim
  // since it is used to set the dimension on C.
  unsigned short int tmp = nsd;
  const_cast<NonlinearElasticityUL*>(this)->nsd = ndim;
  if (!this->formCmatrix(C,X)) return false;
  const_cast<NonlinearElasticityUL*>(this)->nsd = tmp;

  if (calcStress)
    // Evaluate the symmetric stress tensor
    C.multiply(eps,sigma); // sigma = C*eps

  if (calcStress == 1)
  {
    // Push-forward the stress tensor to current configuration
    Tensor dUdX(ndim);
    sigma.transform(dUdX=F); // sigma = F * sigma * F^t
    sigma *= 1.0/J;
#if INT_DEBUG > 0
    std::cout <<"NonlinearElasticityUL::sigma =\n"<< sigma;
#endif
  }
  else if (calcStress == 2)
    return true;

  // Push-forward the constitutive matrix to current configuration

  size_t i, j, k, l, m = C.rows();
  Matrix T(m,m), Ctmp;

  for (i = 1; i <= ndim; i++)
  {
    for (j = 1; j <= ndim; j++)
      T(i,j) = F(j,i)*F(j,i);
    for (j = ndim+1, k = 1; j <= m; j++, k++)
      T(i,j) = 0.5*(F(k,i)*F(k%3+1,i) + F(k%3+1,i)*F(k,i));
  }

  for (i = ndim+1, k = 1; i <= m; i++, k++)
  {
    for (j = 1; j <= ndim; j++)
      T(i,j) = F(j,k)*F(j,k%3+1) + F(j,k%3+1)*F(j,k);
    for (j = ndim+1, l = 1; j <= m; j++, l++)
      T(i,j) = 0.5*(F(l,k)*F(l%3+1,k%3+1) + F(l%3+1,k)*F(l,k%3+1) +
		    F(l,k%3+1)*F(l%3+1,k) + F(l%3+1,k%3+1)*F(l,k));
  }

  // C = 1/J * T^t * C * T
  Ctmp.multiply(C,T);
  C.multiply(T,Ctmp,true);
  C *= 1.0/J;
#if INT_DEBUG > 0
  std::cout <<"NonlinearElasticityUL::C ="<< C;
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
  SymmTensor E(dNdX.cols());
  if (!ulp->kinematics(dNdX,E))
    return false;

  // Evaluate the 2nd Piola Kirchhoff stresses, S
  Matrix C;
  SymmTensor S(E.dim());
  double U = 0.0;
  if (!ulp->constitutive(C,S,U,E,X,2))
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
