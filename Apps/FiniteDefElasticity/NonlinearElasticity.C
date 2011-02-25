// $Id$
//==============================================================================
//!
//! \file NonlinearElasticity.C
//!
//! \date May 25 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for nonlinear elasticity problems.
//!
//==============================================================================

#include "NonlinearElasticity.h"
#include "Utilities.h"
#include "Profiler.h"

#ifndef NO_FORTRAN
#if defined(_WIN32)
#define stiff_tl2d_       STIFF_TL2D
#define stiff_tl3d_       STIFF_TL3D
#define stiff_tl3d_isoel_ STIFF_TL3D_ISOEL
#define stiff_tl3d_isoel_ STIFF_TL3D_ISOEL
#elif defined(_AIX)
#define stiff_tl2d_       stiff_tl2d
#define stiff_tl3d_       stiff_tl3d
#define stiff_tl2d_isoel_ stiff_tl2d_isoel
#define stiff_tl3d_isoel_ stiff_tl3d_isoel
#endif

extern "C" {
  //! \brief Calculates material stiffness contributions for 2D problems.
  void stiff_tl2d_(const int& nenod, const double& detJW,
		   const double* dNdX, const double* F, const double* Cmat,
		   double* EM);
  //! \brief Calculates material stiffness contributions for 3D problems.
  void stiff_tl3d_(const int& nenod, const double& detJW,
		   const double* dNdX, const double* F, const double* Cmat,
		   double* EM);
  //! \brief Optimized for isotropic linear elastic materials in 2D.
  void stiff_tl2d_isoel_(const int& nenod, const double& detJW,
			 const double* dNdX, const double* F,
			 const double& C1, const double& C2,const double& C3,
			 double* EM);
  //! \brief Optimized for isotropic linear elastic materials in 3D.
  void stiff_tl3d_isoel_(const int& nenod, const double& detJW,
			 const double* dNdX, const double* F,
			 const double& C1, const double& C2,const double& C3,
			 double* EM);
}
#endif


NonlinearElasticity::NonlinearElasticity (unsigned short int n, bool ps)
  : NonlinearElasticityTL(n,ps), E(n)
{
  fullCmat = false;
}


void NonlinearElasticity::print (std::ostream& os) const
{
  this->Elasticity::print(os);
  std::cout <<"NonlinearElasticity: Total Lagrangian formulation "
	    <<" (tensorial form)"<< std::endl;
}


void NonlinearElasticity::setMode (SIM::SolutionMode mode)
{
  this->NonlinearElasticityTL::setMode(mode);
  formB = false; // We don't need the B-matrix in the tensor formulation
}


bool NonlinearElasticity::evalInt (LocalIntegral*& elmInt, double detJW,
				   const Vector& N, const Matrix& dNdX,
				   const Vec3& X) const
{
  PROFILE3("NonlinearEl::evalInt");

  // Evaluate the kinematic quantities, F and E, at this point
  if (!this->kinematics(dNdX,E))
    return false;

  // Evaluate current tangent at this point, that is
  // the incremental constitutive matrix, Cmat, and
  // the 2nd Piola-Kirchhoff stress tensor, S
  static SymmTensor S(nsd);
  if (!this->formTangent(Cmat,S,X))
    return false;

  bool haveStrains = !E.isZero(1.0e-16);

  size_t             a, b;
  unsigned short int i, j, k;

  if (eKm)
  {
    // Integrate the material stiffness matrix
#ifndef NO_FORTRAN
    // Using Fortran routines optimized for symmetric constitutive tensors
    PROFILE4("stiff_TL_");
    if (nsd == 3)
      if (fullCmat)
	stiff_tl3d_(N.size(),detJW,dNdX.ptr(),F.ptr(),
		    Cmat.ptr(),eKm->ptr());
      else
	stiff_tl3d_isoel_(N.size(),detJW,dNdX.ptr(),F.ptr(),
			  Cmat(1,1),Cmat(1,2),Cmat(4,4),eKm->ptr());
    else if (nsd == 2)
      if (fullCmat)
	stiff_tl2d_(N.size(),detJW,dNdX.ptr(),F.ptr(),
		    Cmat.ptr(),eKm->ptr());
      else
	stiff_tl2d_isoel_(N.size(),detJW,dNdX.ptr(),F.ptr(),
			  Cmat(1,1),Cmat(1,2),Cmat(3,3),eKm->ptr());
    else if (nsd == 1)
      for (a = 1; a <= N.size(); a++)
	for (b = 1; b <= N.size(); b++)
	  (*eKm)(a,b) += dNdX(a,1)*F(1,1)*Cmat(1,1)*F(1,1)*dNdX(b,1);
#else
    // This is too costly, but is basically what is done in the fortran routines
    PROFILE4("dNdX^t*F^t*C*F*dNdX");
    unsigned short int l, m, n;
    SymmTensor4 D(Cmat,nsd); // fourth-order material tensor
    Matrix& EM = *eKm;
    for (a = 1; a <= N.size(); a++)
      for (b = 1; b <= N.size(); b++)
	for (m = 1; m <= nsd; m++)
	  for (n = 1; n <= nsd; n++)
	  {
	    double km = 0.0;
	    for (i = 1; i <= nsd; i++)
	      for (j = 1; j <= nsd; j++)
		for (k = 1; k <= nsd; k++)
		  for (l = 1; l <= nsd; l++)
		    km += dNdX(a,i)*F(m,j)*D(i,j,k,l)*F(n,k)*dNdX(b,l);

	    EM(nsd*(a-1)+m,nsd*(b-1)+n) += km*detJW;
	  }
#endif
  }

  if (eKg && haveStrains)
    // Integrate the geometric stiffness matrix
    this->formKG(*eKg,dNdX,S,detJW);

  if (eM)
    // Integrate the mass matrix
    this->formMassMatrix(*eM,N,X,detJW);

  if (iS && haveStrains)
  {
    // Integrate the internal forces
    Vector& ES = *iS;
    for (a = 1; a <= N.size(); a++)
      for (k = 1; k <= nsd; k++)
      {
	double f = 0.0;
	for (i = 1; i <= nsd; i++)
	  for (j = 1; j <= nsd; j++)
	    f -= dNdX(a,i)*F(k,j)*S(i,j);
	ES(nsd*(a-1)+k) += f*detJW;
      }
  }

  if (eS)
    // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(*eS,N,X,detJW);

  return this->getIntegralResult(elmInt);
}


bool NonlinearElasticity::evalSol (Vector& s, const Vector&,
				   const Matrix& dNdX, const Vec3& X,
				   const std::vector<int>& MNPC) const
{
  PROFILE3("NonlinearEl::evalSol");

  int ierr = 0;
  if (eV && !primsol.front().empty())
    if ((ierr = utl::gather(MNPC,nsd,primsol.front(),*eV)))
    {
      std::cerr <<" *** NonlinearElasticity::evalSol: Detected "
		<< ierr <<" node numbers out of range."<< std::endl;
      return false;
    }

  // Evaluate the stress state at this point
  static SymmTensor Sigma(nsd);
  if (!this->formStressTensor(dNdX,X,Sigma))
    return false;

  // Congruence transformation to local coordinate system at current point
  if (locSys) Sigma.transform(locSys->getTmat(X));

  s = Sigma;
  s.push_back(Sigma.vonMises());
  return true;
}


bool NonlinearElasticity::evalSol (Vector& s,
				   const Matrix& dNdX, const Vec3& X) const
{
  PROFILE3("NonlinearEl::evalSol");

  static SymmTensor Sigma(nsd);
  if (!this->formStressTensor(dNdX,X,Sigma))
    return false;

  s = Sigma;
  return true;
}


bool NonlinearElasticity::formStressTensor (const Matrix& dNdX, const Vec3& X,
					    SymmTensor& S) const
{
  if (!eV || eV->empty())
  {
    // Initial state (zero stresses)
    S.zero();
    return true;
  }

  // Evaluate the kinematic quantities, F and E, at this point
  if (!this->kinematics(dNdX,E))
    return false;

  // Evaluate the constitutive matrix, C, at this point
  if (!this->formCmatrix(Cmat,X))
    return false;

  return Cmat.multiply(E,S); // S = C*E
}


bool NonlinearElasticity::formTangent (Matrix& Ctan, SymmTensor& S,
				       const Vec3& X) const
{
  if (!this->formCmatrix(Ctan,X))
    return false;
  else if (eV && !eV->empty())
    return Ctan.multiply(E,S); // S = C*E
  else
    return true; // Initial state (no stresses)
}
