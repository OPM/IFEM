// $Id$
//==============================================================================
//!
//! \file AxSymmElasticity.C
//!
//! \date Apr 28 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for axial-symmetric elasticity problems.
//!
//==============================================================================

#include "AxSymmElasticity.h"
#include "FiniteElement.h"
#include "MaterialBase.h"
#include "Utilities.h"
#include "Tensor.h"
#include "Vec3Oper.h"

#ifndef epsR
//! \brief Zero tolerance for the radial coordinate.
#define epsR 1.0e-16
#endif


AxSymmElasticity::AxSymmElasticity () : Elasticity(2)
{
  // Only the current solution is needed
  primsol.resize(1);
}


void AxSymmElasticity::print (std::ostream& os) const
{
  std::cout <<"Axial-symmetric Elasticity problem\n";
  this->Elasticity::print(os);
}


/*!
  The strain-displacement matrix for an axiallly symmetric 3D continuum element
  is formally defined as:
  \f[
  [B] = \left[\begin{array}{cc}
  \frac{\partial}{\partial r} &                0            \\
                 0            & \frac{\partial}{\partial z} \\
         \frac{1}{r}          &                0            \\
  \frac{\partial}{\partial z} & \frac{\partial}{\partial r}
  \end{array}\right] [N]
  \f]
  where
  [\a N ] is the element basis functions arranged in a [2][2*NENOD] matrix.
*/

bool AxSymmElasticity::kinematics (const Vector& N, const Matrix& dNdX,
				   double r, SymmTensor& eps) const
{
  const size_t nenod = N.size();
  Bmat.resize(8,nenod,true);
  if (dNdX.cols() < 2)
  {
    std::cerr <<" *** AxSymmElasticity::kinematics: Invalid dimension on dNdX, "
	      << dNdX.rows() <<"x"<< dNdX.cols() <<"."<< std::endl;
    return false;
  }
  else if (r < -epsR)
  {
    std::cerr <<" *** AxSymmElasticity::kinematics: Invalid point r < 0, "
	      << r << std::endl;
    return false;
  }

#define INDEX(i,j) i+4*(j-1)

  // Strain-displacement matrix for 3D axisymmetric elements:
  //
  //         | d/dr   0   |
  //   [B] = |  0    d/dz | * [N]
  //         | 1/r    0   |
  //         | d/dz  d/dr |

  for (size_t i = 1; i <= nenod; i++)
  {
    // Normal strain part
    Bmat(INDEX(1,1),i) = dNdX(i,1);
    Bmat(INDEX(2,2),i) = dNdX(i,2);
    // Hoop strain part
    Bmat(INDEX(3,1),i) = r <= epsR ? dNdX(i,1) : N(i)/r;
    // Shear strain part
    Bmat(INDEX(4,1),i) = dNdX(i,2);
    Bmat(INDEX(4,2),i) = dNdX(i,1);
  }

  Bmat.resize(4,2*nenod);

  // Evaluate the strains
  if (eV && !eV->empty() && eps.dim() > 0)
    return Bmat.multiply(*eV,eps); // eps = B*eV

  return true;
}


bool AxSymmElasticity::evalInt (LocalIntegral*& elmInt, const FiniteElement& fe,
				const Vec3& X) const
{
  bool lHaveStrains = false;
  SymmTensor eps(2,true), sigma(2,true);

  if (eKm || eKg || iS)
  {
    // Compute the strain-displacement matrix B from N, dNdX and r = X.x,
    // and evaluate the symmetric strain tensor if displacements are available
    if (!this->kinematics(fe.N,fe.dNdX,X.x,eps)) return false;
    if (!eps.isZero(1.0e-16)) lHaveStrains = true;

    // Evaluate the constitutive matrix and the stress tensor at this point
    double U;
    if (!material->evaluate(Cmat,sigma,U,X,eps,eps))
      return false;
  }

  // Axi-symmetric integration point volume; 2*pi*r*|J|*w
  const double detJW = 2.0*M_PI*X.x*fe.detJxW;

  if (eKm)
  {
    // Integrate the material stiffness matrix
    CB.multiply(Cmat,Bmat).multiply(detJW); // CB = C*B*|J|*w
    eKm->multiply(Bmat,CB,true,false,true); // EK += B^T * CB
  }

  if (eKg && lHaveStrains)
    // Integrate the geometric stiffness matrix
    this->formKG(*eKg,fe.dNdX,sigma,detJW);

  if (eM)
    // Integrate the mass matrix
    this->formMassMatrix(*eM,fe.N,X,detJW);

  if (iS && lHaveStrains)
  {
    // Integrate the internal forces
    sigma *= -detJW;
    if (!Bmat.multiply(sigma,*iS,true,true)) // ES -= B^T*sigma
      return false;
  }

  if (eS)
    // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(*eS,fe.N,X,detJW);

  return this->getIntegralResult(elmInt);
}


bool AxSymmElasticity::evalBou (LocalIntegral*& elmInt, const FiniteElement& fe,
				const Vec3& X, const Vec3& normal) const
{
  if (!tracFld)
  {
    std::cerr <<" *** AxSymmElasticity::evalBou: No tractions."<< std::endl;
    return false;
  }
  else if (!eS)
  {
    std::cerr <<" *** AxSymmElasticity::evalBou: No load vector."<< std::endl;
    return false;
  }

  // Axi-symmetric integration point volume; 2*pi*r*|J|*w
  const double detJW = 2.0*M_PI*X.x*fe.detJxW;

  // Evaluate the surface traction
  Vec3 T = (*tracFld)(X,normal);

  // Store the traction value for vizualization
  if (!T.isZero()) tracVal[X] = T;

  // Integrate the force vector
  Vector& ES = *eS;
  for (size_t a = 1; a <= fe.N.size(); a++)
  {
    ES(2*a-1) += T.x*fe.N(a)*detJW;
    ES(2*a  ) += T.y*fe.N(a)*detJW;
  }

  return this->getIntegralResult(elmInt);
}


bool AxSymmElasticity::evalSol (Vector& s, const Vector& N,
				const Matrix& dNdX, const Vec3& X,
				const std::vector<int>& MNPC) const
{
  int ierr = 0;
  if (eV && !primsol.front().empty())
    if ((ierr = utl::gather(MNPC,2,primsol.front(),*eV)))
    {
      std::cerr <<" *** AxSymmElasticity::evalSol: Detected "
		<< ierr <<" node numbers out of range."<< std::endl;
      return false;
    }

  // Evaluate the strain tensor, eps
  SymmTensor eps(2,true);
  if (!this->kinematics(N,dNdX,X.x,eps))
    return false;

  // Calculate the stress tensor through the constitutive relation
  SymmTensor sigma(2,true); double U;
  if (!material->evaluate(Cmat,sigma,U,X,eps,eps))
    return false;

  // Congruence transformation to local coordinate system at current point
  if (locSys) sigma.transform(locSys->getTmat(X));

  s = sigma;
  s.push_back(sigma.vonMises());
  return true;
}


bool AxSymmElasticity::evalSol (Vector& s, const Vector& N,
				const Matrix& dNdX, const Vec3& X) const
{
  if (!eV || eV->empty())
  {
    std::cerr <<" *** AxSymmElasticity::evalSol: No displacement vector."
	      << std::endl;
    return false;
  }
  else if (eV->size() != dNdX.rows()*2)
  {
    std::cerr <<" *** AxSymmElasticity::evalSol: Invalid displacement vector."
	      <<"\n     size(eV) = "<< eV->size() <<"   size(dNdX) = "
	      << dNdX.rows() <<","<< dNdX.cols() << std::endl;
    return false;
  }

  // Evaluate the strain tensor, eps
  SymmTensor eps(2,true);
  if (!this->kinematics(N,dNdX,X.x,eps))
    return false;

  // Calculate the stress tensor through the constitutive relation
  SymmTensor sigma(2,true); double U;
  if (!material->evaluate(Cmat,sigma,U,X,eps,eps))
    return false;

  s = sigma;
  return true;
}


const char* AxSymmElasticity::getField1Name (size_t i, const char* prefix) const
{
  if (i > 2)
    return "displacement";
  else if (!prefix)
    return i == 0 ? "u_r" : "u_z";

  static std::string name;
  name = prefix + std::string(i == 0 ? " u_r" : " u_z");

  return name.c_str();
}


const char* AxSymmElasticity::getField2Name (size_t i, const char* prefix) const
{
  if (i > 4) return 0;

  static const char* s[5] = { "s_rr","s_zz","s_tt","s_zr", "von Mises stress" };
  if (!prefix) return s[i];

  static std::string name;
  name = prefix + std::string(" ") + s[i];

  return name.c_str();
}
