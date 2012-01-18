// $Id$
//==============================================================================
//!
//! \file NonlinearElasticityTL.C
//!
//! \date May 25 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for nonlinear elasticity problems.
//!
//==============================================================================

#include "NonlinearElasticityTL.h"
#include "MaterialBase.h"
#include "FiniteElement.h"
#include "ElmMats.h"
#include "Tensor.h"
#include "Vec3Oper.h"

#ifndef epsR
//! \brief Zero tolerance for the radial coordinate.
#define epsR 1.0e-16
#endif


NonlinearElasticityTL::NonlinearElasticityTL (unsigned short int n, bool axS)
  : Elasticity(n,axS)
{
  // Only the current solution is needed
  primsol.resize(1);
}


void NonlinearElasticityTL::print (std::ostream& os) const
{
  std::cout <<"NonlinearElasticityTL: Total Lagrangian formulation"<< std::endl;

  this->Elasticity::print(os);
}


void NonlinearElasticityTL::setMode (SIM::SolutionMode mode)
{
  formB = true;
  m_mode = mode;
  eM = eKm = eKg = 0;
  eS = iS  = 0;

  switch (mode)
    {
    case SIM::STATIC:
    case SIM::RHS_ONLY:
      eKm = eKg = 1;
      break;

    case SIM::DYNAMIC:
      eKm = eKg = 1;
      eM  = 2;
      break;

    default:
      formB = false;
      this->Elasticity::setMode(mode);
      return;
    }

  // We always need the force vectors in nonlinear simulations
  iS = eS = 1;
  tracVal.clear();
}


LocalIntegral* NonlinearElasticityTL::getLocalIntegral (size_t nen, size_t,
							bool neumann) const
{
  ElmMats* result = new ElmMats;
  switch (m_mode)
  {
    case SIM::STATIC:
      result->rhsOnly = neumann;
      result->withLHS = !neumann;
      result->resize(neumann?0:1,1);
      break;

    case SIM::DYNAMIC:
      result->rhsOnly = neumann;
      result->withLHS = !neumann;
      result->resize(neumann?0:2,1);
      break;

    case SIM::VIBRATION:
    case SIM::BUCKLING:
      result->withLHS = true;
      result->resize(2,0);
      break;

    case SIM::STIFF_ONLY:
    case SIM::MASS_ONLY:
      result->withLHS = true;
      result->resize(1,0);
      break;

    case SIM::RHS_ONLY:
      result->rhsOnly = true;
      result->resize(neumann?0:1,1);
      break;

    case SIM::RECOVERY:
      result->rhsOnly = true;
      break;

    default:
      ;
  }

  for (size_t i = 0; i < result->A.size(); i++)
    result->A[i].resize(nsd*nen,nsd*nen);

  if (result->b.size())
    result->b.front().resize(nsd*nen);

  return result;
}


bool NonlinearElasticityTL::evalInt (LocalIntegral& elmInt,
				     const FiniteElement& fe,
				     const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  // Evaluate the deformation gradient, F, and the Green-Lagrange strains, E,
  // and compute the nonlinear strain-displacement matrix, B, from dNdX and F
  Matrix Bmat;
  Tensor F(nDF);
  SymmTensor E(nsd,axiSymmetry), S(nsd,axiSymmetry);
  if (!this->kinematics(elMat.vec,fe.N,fe.dNdX,X.x,Bmat,F,E))
    return false;

  // Evaluate the constitutive relation
  Matrix Cmat;
  bool lHaveStrains = !E.isZero(1.0e-16);
  if (eKm || eKg || iS)
  {
    double U;
    if (!material->evaluate(Cmat,S,U,X,F,E, (eKg || iS) && lHaveStrains ? 2:0))
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
    // Integrate the geometric stiffness matrix
    this->formKG(elMat.A[eKg-1],fe.N,fe.dNdX,X.x,S,detJW);

  if (eM)
    // Integrate the mass matrix
    this->formMassMatrix(elMat.A[eM-1],fe.N,X,detJW);

  if (iS && lHaveStrains)
  {
    // Integrate the internal forces
    S *= -detJW;
    if (!Bmat.multiply(S,elMat.b[iS-1],true,true)) // ES -= B^T*S
      return false;
  }

  if (eS)
    // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(elMat.b[eS-1],fe.N,X,detJW);

  return true;
}


bool NonlinearElasticityTL::evalBou (LocalIntegral& elmInt,
				     const FiniteElement& fe,
				     const Vec3& X, const Vec3& normal) const
{
  if (!tracFld)
  {
    std::cerr <<" *** NonlinearElasticityTL::evalBou: No tractions."
	      << std::endl;
    return false;
  }
  else if (!eS)
  {
    std::cerr <<" *** NonlinearElasticityTL::evalBou: No load vector."
	      << std::endl;
    return false;
  }

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  // Evaluate the surface traction
  Vec3 T = (*tracFld)(X,normal);

  // Store the traction value for vizualization
  if (!T.isZero()) tracVal[X] = T;

  // Check for with-rotated pressure load
  unsigned short int i, j;
  if (tracFld->isNormalPressure())
  {
    // Compute the deformation gradient, F
    Matrix B;
    Tensor F(nDF);
    SymmTensor dummy(0);
    if (!this->kinematics(elMat.vec,fe.N,fe.dNdX,X.x,B,F,dummy)) return false;

    // Compute its inverse and determinant, J
    double J = F.inverse();
    if (J <= 0.0) return false;

    // Pull-back the normal traction to the initial configuration.
    // See equation (3.4.5) on page 102 in T. Belytschko et. al (2000):
    //    p*n*dS ==> J*p*n0*F^-1*dS0
    Vec3 t = J*T; T = 0.0;
    for (i = 1; i <= nsd; i++)
      for (j = 1; j <= nsd; j++)
	T[i-1] += t[j-1]*F(j,i);
  }

  // Axi-symmetric integration point volume; 2*pi*r*|J|*w
  const double detJW = axiSymmetry ? 2.0*M_PI*X.x*fe.detJxW : fe.detJxW;

  // Integrate the force vector
  Vector& ES = elMat.b[eS-1];
  for (size_t a = 1; a <= fe.N.size(); a++)
    for (i = 1; i <= nsd; i++)
      ES(nsd*(a-1)+i) += T[i-1]*fe.N(a)*detJW;

  return true;
}


bool NonlinearElasticityTL::kinematics (const Vectors& eV,
					const Vector& N, const Matrix& dNdX,
					double r, Matrix& Bmat, Tensor& F,
					SymmTensor& E) const
{
  if (eV.empty() || eV.front().empty())
  {
    // Initial state, unit deformation gradient and linear B-matrix
    F = 1.0;
    E.zero();
    if (axiSymmetry)
      return this->formBmatrix(Bmat,N,dNdX,r);
    else
      return this->formBmatrix(Bmat,dNdX);
  }

  const size_t nenod = dNdX.rows();
  const size_t nstrc = axiSymmetry ? 4 : nsd*(nsd+1)/2;
  if (eV.front().size() != nenod*nsd || dNdX.cols() < nsd)
  {
    std::cerr <<" *** NonlinearElasticityTL::kinematics: Invalid dimension,"
	      <<" dNdX("<< nenod <<","<< dNdX.cols() <<")"<< std::endl;
    return false;
  }

  // Compute the deformation gradient, [F] = [I] + [dudX] = [I] + [dNdX]*[u],
  // and the Green-Lagrange strains, E_ij = 0.5(F_ij+F_ji+F_ki*F_kj).

  // Notice that the matrix multiplication method used here treats the
  // element displacement vector, *eV, as a matrix whose number of columns
  // equals the number of rows in the matrix dNdX.
  Matrix dUdX;
  if (!dUdX.multiplyMat(eV.front(),dNdX)) // dUdX = Grad{u} = eV*dNd
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
  // TODO: How is this for axisymmetric problems?
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
  // Add the dU/r term to the F(3,3)-term for axisymmetric problems
  if (axiSymmetry && r > epsR) F(3,3) += eV.front().dot(N,0,nsd)/r;

#ifdef INT_DEBUG
  std::cout <<"NonlinearElasticityTL::F =\n"<< F;
#endif

  if (!formB || E.dim() < nsd) return true;

  // Form the nonlinear B-matrix

  Bmat.resize(nstrc*nsd,nenod,true);

#define INDEX(i,j) i+nstrc*(j-1)

  // Normal strain part
  size_t a;
  for (a = 1; a <= nenod; a++)
    for (i = 1; i <= nsd; i++)
      for (j = 1; j <= nsd; j++)
	Bmat(INDEX(j,i),a) = F(i,j)*dNdX(a,j);

  // Shear strain part
  if (nsd == 3)
    for (a = 1; a <= nenod; a++)
      for (i = 1; i <= nsd; i++)
      {
	Bmat(INDEX(4,i),a) = F(i,1)*dNdX(a,2) + F(i,2)*dNdX(a,1);
	Bmat(INDEX(5,i),a) = F(i,2)*dNdX(a,3) + F(i,3)*dNdX(a,2);
	Bmat(INDEX(6,i),a) = F(i,3)*dNdX(a,1) + F(i,1)*dNdX(a,3);
      }

  else if (nsd == 2)
    for (a = 1; a <= nenod; a++)
      if (axiSymmetry)
      {
	for (i = 1; i <= nsd; i++)
	  Bmat(INDEX(4,i),a) = F(i,1)*dNdX(a,2) + F(i,2)*dNdX(a,1);
	// Hoop strain part for axisymmetry (TODO: check this)
	Bmat(INDEX(3,1),i) = F(3,3) * (r <= epsR ? dNdX(a,1) : N(a)/r);
      }
      else
	for (i = 1; i <= nsd; i++)
	  Bmat(INDEX(3,i),a) = F(i,1)*dNdX(a,2) + F(i,2)*dNdX(a,1);

  Bmat.resize(nstrc,nsd*nenod);
#ifdef INT_DEBUG
  std::cout <<"NonlinearElasticityTL::B ="<< Bmat;
#endif
  return true;
}
