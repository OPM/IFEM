// $Id: NonlinearElasticityTL.C,v 1.4 2011-02-08 09:06:02 kmo Exp $
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

#include "NonlinearElasticity.h"
#include "ElmMats.h"
#include "Tensor.h"
#include "Vec3Oper.h"


NonlinearElasticityTL::NonlinearElasticityTL (unsigned short int n, bool ps)
  : Elasticity(n,ps)
{
  // Only the current solution is needed
  primsol.resize(1);
}


void NonlinearElasticityTL::print (std::ostream& os) const
{
  this->Elasticity::print(os);
  std::cout <<"NonlinearElasticityTL: Total Lagrangian formulation"<< std::endl;
}


void NonlinearElasticityTL::setMode (SIM::SolutionMode mode)
{
  if (!myMats) return;

  formB = true;
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
      formB = false;
      this->Elasticity::setMode(mode);
      return;
    }

  // We always need the force+displacement vectors in nonlinear simulations
  iS = &myMats->b[0];
  eS = &myMats->b[0];
  eV = &myMats->b[1];
  tracVal.clear();
}


bool NonlinearElasticityTL::evalBou (LocalIntegral*& elmInt, double detJW,
				     const Vector& N, const Matrix& dNdX,
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

  // Evaluate the surface traction
  Vec3 T = (*tracFld)(X,normal);

  // Store the traction value for vizualization
  if (!T.isZero()) tracVal[X] = T;

  // Check for with-rotated pressure load
  unsigned short int i, j;
  if (tracFld->isNormalPressure())
  {
    // Compute the deformation gradient, F
    SymmTensor dummy(0);
    if (!this->kinematics(dNdX,dummy)) return false;

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

  // Integrate the force vector
  Vector& ES = *eS;
  for (size_t a = 1; a <= N.size(); a++)
    for (i = 1; i <= nsd; i++)
      ES(nsd*(a-1)+i) += T[i-1]*N(a)*detJW;

  return this->getIntegralResult(elmInt);
}


bool NonlinearElasticityTL::kinematics (const Matrix& dNdX, SymmTensor& E) const
{
  if (!eV || eV->empty())
  {
    // Initial state, unit deformation gradient and linear B-matrix
    F.diag(1.0,nsd);
    return this->Elasticity::kinematics(dNdX,E);
  }

  const size_t nenod = dNdX.rows();
  const size_t nstrc = nsd*(nsd+1)/2;
  if (eV->size() != nenod*nsd || dNdX.cols() < nsd)
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

#ifdef INT_DEBUG
  std::cout <<"NonlinearElasticityTL::F ="<< F;
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
      for (i = 1; i <= nsd; i++)
	Bmat(INDEX(3,i),a) = F(i,1)*dNdX(a,2) + F(i,2)*dNdX(a,1);

  Bmat.resize(nstrc,nsd*nenod);
#ifdef INT_DEBUG
  std::cout <<"NonlinearElasticityTL::B ="<< Bmat;
#endif
  return true;
}
