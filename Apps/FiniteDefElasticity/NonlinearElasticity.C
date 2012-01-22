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
#include "MaterialBase.h"
#include "FiniteElement.h"
#include "ElmMats.h"
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


NonlinearElasticity::NonlinearElasticity (unsigned short int n)
  : NonlinearElasticityTL(n)
{
  fullCmat = false;
}


void NonlinearElasticity::print (std::ostream& os) const
{
  std::cout <<"NonlinearElasticity: Total Lagrangian formulation "
	    <<" (tensorial form)"<< std::endl;

  this->Elasticity::print(os);
}


void NonlinearElasticity::setMode (SIM::SolutionMode mode)
{
  this->NonlinearElasticityTL::setMode(mode);
  formB = false; // We don't need the B-matrix in the tensor formulation
}


bool NonlinearElasticity::evalInt (LocalIntegral& elmInt,
				   const FiniteElement& fe,
				   const Vec3& X) const
{
  PROFILE3("NonlinearEl::evalInt");

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  // Evaluate the kinematic quantities, F and E, at this point
  Matrix Bmat;
  Tensor F(nsd);
  SymmTensor E(nsd);
  if (!this->kinematics(elMat.vec.front(),fe.N,fe.dNdX,X.x,Bmat,F,E))
    return false;

  // Evaluate current tangent at this point, that is
  // the incremental constitutive matrix, Cmat, and
  // the 2nd Piola-Kirchhoff stress tensor, S
  Matrix Cmat;
  SymmTensor S(nsd);
  if (!this->formTangent(Cmat,S,fe.iGP,X,F,E))
    return false;

  bool haveStrains = !E.isZero(1.0e-16);

  size_t             a, b;
  unsigned short int i, j, k;

  if (eKm)
  {
    // Integrate the material stiffness matrix
    Matrix& EM = elMat.A[eKm-1];
#ifndef NO_FORTRAN
    // Using Fortran routines optimized for symmetric constitutive tensors
    PROFILE4("stiff_TL_");
    if (nsd == 3)
      if (fullCmat)
	stiff_tl3d_(fe.N.size(),fe.detJxW,fe.dNdX.ptr(),F.ptr(),
		    Cmat.ptr(),EM.ptr());
      else
	stiff_tl3d_isoel_(fe.N.size(),fe.detJxW,fe.dNdX.ptr(),F.ptr(),
			  Cmat(1,1),Cmat(1,2),Cmat(4,4),EM.ptr());
    else if (nsd == 2)
      if (fullCmat)
	stiff_tl2d_(fe.N.size(),fe.detJxW,fe.dNdX.ptr(),F.ptr(),
		    Cmat.ptr(),EM.ptr());
      else
	stiff_tl2d_isoel_(fe.N.size(),fe.detJxW,fe.dNdX.ptr(),F.ptr(),
			  Cmat(1,1),Cmat(1,2),Cmat(3,3),EM.ptr());
    else if (nsd == 1)
      for (a = 1; a <= fe.N.size(); a++)
	for (b = 1; b <= fe.N.size(); b++)
	  EM(a,b) += fe.dNdX(a,1)*F(1,1)*Cmat(1,1)*F(1,1)*fe.dNdX(b,1);
#else
    // This is too costly, but is basically what is done in the fortran routines
    PROFILE4("dNdX^t*F^t*C*F*dNdX");
    unsigned short int l, m, n;
    SymmTensor4 D(Cmat,nsd); // fourth-order material tensor
    for (a = 1; a <= fe.N.size(); a++)
      for (b = 1; b <= fe.N.size(); b++)
	for (m = 1; m <= nsd; m++)
	  for (n = 1; n <= nsd; n++)
	  {
	    double km = 0.0;
	    for (i = 1; i <= nsd; i++)
	      for (j = 1; j <= nsd; j++)
		for (k = 1; k <= nsd; k++)
		  for (l = 1; l <= nsd; l++)
		    km += fe.dNdX(a,i)*F(m,j)*D(i,j,k,l)*F(n,k)*fe.dNdX(b,l);

	    EM(nsd*(a-1)+m,nsd*(b-1)+n) += km*fe.detJxW;
	  }
#endif
  }

  if (eKg && haveStrains)
    // Integrate the geometric stiffness matrix
    this->formKG(elMat.A[eKg-1],fe.N,fe.dNdX,X.x,S,fe.detJxW);

  if (eM)
    // Integrate the mass matrix
    this->formMassMatrix(elMat.A[eM-1],fe.N,X,fe.detJxW);

  if (iS && haveStrains)
  {
    // Integrate the internal forces
    Vector& ES = elMat.b[iS-1];
    for (a = 1; a <= fe.N.size(); a++)
      for (k = 1; k <= nsd; k++)
      {
	double f = 0.0;
	for (i = 1; i <= nsd; i++)
	  for (j = 1; j <= nsd; j++)
	    f -= fe.dNdX(a,i)*F(k,j)*S(i,j);
	ES(nsd*(a-1)+k) += f*fe.detJxW;
      }
  }

  if (eS)
    // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(elMat.b[eS-1],fe.N,X,fe.detJxW);

  return true;
}


bool NonlinearElasticity::evalSol (Vector& s, const Vector&,
				   const Matrix& dNdX, const Vec3& X,
				   const std::vector<int>& MNPC) const
{
  PROFILE3("NonlinearEl::evalSol");

  // Extract element displacements
  Vector eV;
  int ierr = 0;
  if (!primsol.front().empty())
    if ((ierr = utl::gather(MNPC,nsd,primsol.front(),eV)))
    {
      std::cerr <<" *** NonlinearElasticity::evalSol: Detected "<< ierr
		<<" node numbers out of range."<< std::endl;
      return false;
    }

  // Evaluate the stress state at this point
  SymmTensor Sigma(nsd);
  if (!this->formStressTensor(eV,dNdX,X,Sigma))
    return false;

  // Congruence transformation to local coordinate system at current point
  if (locSys) Sigma.transform(locSys->getTmat(X));

  s = Sigma;
  s.push_back(Sigma.vonMises());
  return true;
}


bool NonlinearElasticity::evalSol (Vector& s, const Vector& eV,
				   const Matrix& dNdX, const Vec3& X) const
{
  PROFILE3("NonlinearEl::evalSol");

  SymmTensor Sigma(nsd);
  if (!this->formStressTensor(eV,dNdX,X,Sigma))
    return false;

  s = Sigma;
  return true;
}


bool NonlinearElasticity::formStressTensor (const Vector& eV,
					    const Matrix& dNdX, const Vec3& X,
					    SymmTensor& S) const
{
  if (eV.empty())
  {
    // Initial state (zero stresses)
    S.zero();
    return true;
  }

  // Evaluate the kinematic quantities, F and E, at this point
  Matrix B;
  Tensor F(nsd);
  SymmTensor E(nsd);
  if (!this->kinematics(eV,Vector(),dNdX,X.x,B,F,E))
    return false;

  // Evaluate the 2nd Piola-Kirchhoff stress tensor, S, at this point
  Matrix Cmat;
  double U;
  return material->evaluate(Cmat,S,U,0,X,F,E,2);
}


bool NonlinearElasticity::formTangent (Matrix& Ctan, SymmTensor& S, size_t iGP,
				       const Vec3& X, const Tensor& F,
				       const SymmTensor& E) const
{
  double U;
  return material->evaluate(Ctan,S,U,iGP,X,F,E, E.isZero(1.0e-16) ? 0 : 2);
}
