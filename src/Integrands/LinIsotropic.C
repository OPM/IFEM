// $Id$
//==============================================================================
//!
//! \file LinIsotropic.C
//!
//! \date Mar 01 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Isotropic linear elastic material model.
//!
//==============================================================================

#include "LinIsotropic.h"
#include "Tensor.h"
#include "Vec3.h"


LinIsotropic::LinIsotropic (bool ps, bool ax) : planeStress(ps), axiSymmetry(ax)
{
  // Default material properties - typical values for steel (SI units)
  Efunc = 0;
  Emod = 2.05e11;
  nu = 0.29;
  rho = 7.85e3;
}


LinIsotropic::LinIsotropic (RealFunc* E, double v, double den, bool ps, bool ax)
  : Efunc(E), nu(v), rho(den), planeStress(ps), axiSymmetry(ax)
{
  Emod = -1.0; // Should not be referenced
}


void LinIsotropic::print (std::ostream& os) const
{
  std::cout <<"LinIsotropic: ";
  if (axiSymmetry)
    std::cout <<"axial-symmetric, ";
  else if (planeStress)
    std::cout <<"plane stress, ";
  std::cout <<"E = "<< Emod <<", nu = "<< nu <<", rho = "<< rho << std::endl;
}


/*!
  The consitutive matrix for Isotropic linear elastic problems
  is defined as follows:

  For 2D plain stress: \f[
  [C] = \frac{E}{(1-\nu^2)} \left[\begin{array}{ccc}
  1 & \ \ \nu & 0 \\
  \nu & \ \ 1 & 0 \\
  0 & \ \ 0 & \frac{1}{2}(1-\nu)
  \end{array}\right] \f]

  For 2D plain strain: \f[
  [C] = \frac{E}{(1+\nu)(1-2\nu)} \left[\begin{array}{ccc}
  1-\nu & \nu & 0 \\
  \nu & 1-\nu & 0 \\
  0 & 0 & \frac{1}{2}-\nu
  \end{array}\right] \f]

  For 3D axisymmetric solids: \f[
  [C] = \frac{E}{(1+\nu)(1-2\nu)} \left[\begin{array}{cccc}
  1-\nu & \nu & \nu & 0 \\
  \nu & 1-\nu & \nu & 0 \\
  \nu & \nu & 1-\nu & 0 \\
  0 & 0 & 0 & \frac{1}{2}-\nu
  \end{array}\right] \f]

  For 3D: \f[
  [C] = \frac{E}{(1+\nu)(1-2\nu)} \left[\begin{array}{cccccc}
  1-\nu & \nu & \nu & 0 & 0 & 0 \\
  \nu & 1-\nu & \nu & 0 & 0 & 0 \\
  \nu & \nu & 1-\nu & 0 & 0 & 0 \\
  0 & 0 & 0 & \frac{1}{2}-\nu & 0 & 0 \\
  0 & 0 & 0 & 0 & \frac{1}{2}-\nu & 0 \\
  0 & 0 & 0 & 0 & 0 & \frac{1}{2}-\nu
  \end{array}\right] \f]
*/

bool LinIsotropic::evaluate (Matrix& C, SymmTensor& sigma, double& U, size_t,
			     const Vec3& X, const Tensor&,
			     const SymmTensor& eps, char iop,
			     const TimeDomain*, const Tensor*) const
{
  const size_t nsd = sigma.dim();
  const size_t nst = nsd == 2 && axiSymmetry ? 4 : nsd*(nsd+1)/2;
  C.resize(nst,nst,true);

  // Evaluate the scalar stiffness function, if defined
  const double E = Efunc ? (*Efunc)(X) : Emod;

  if (nsd == 1)
  {
    // Special for 1D problems
    C(1,1) = iop < 0 ? 1.0/E : E;
    if (iop > 0)
    {
      sigma = eps; sigma *= E;
      if (iop == 3)
	U = 0.5*sigma(1,1)*eps(1,1);
    }
    return true;
  }
  else if (nu < 0.0 || nu >= 0.5)
  {
    std::cerr <<" *** LinIsotropic::evaluate: Poisson's ratio "<< nu
              <<" out of range [0,0.5>."<< std::endl;
    return false;
  }

  if (iop < 0) // The inverse C-matrix is wanted
    if (nsd == 3 || (nsd == 2 && (planeStress || axiSymmetry)))
    {
      C(1,1) = 1.0 / E;
      C(2,1) = -nu / E;
    }
    else // 2D plain strain
    {
      C(1,1) = (1.0 - nu*nu) / E;
      C(2,1) = (-nu - nu*nu) / E;
    }

  else
    if (nsd == 2 && planeStress && !axiSymmetry)
    {
      C(1,1) = E / (1.0 - nu*nu);
      C(2,1) = C(1,1) * nu;
    }
    else // 2D plain strain, axisymmetric or 3D
    {
      double fact = E / ((1.0 + nu) * (1.0 - nu - nu));
      C(1,1) = fact * (1.0 - nu);
      C(2,1) = fact * nu;
    }

  C(1,2) = C(2,1);
  C(2,2) = C(1,1);

  const double G = E / (2.0 + nu + nu);
  C(nsd+1,nsd+1) = iop < 0 ? 1.0 / G : G;

  if (nsd == 2 && axiSymmetry)
  {
    C(4,4) = C(3,3);
    C(3,1) = C(2,1);
    C(3,2) = C(2,1);
    C(1,3) = C(2,1);
    C(2,3) = C(2,1);
    C(3,3) = C(1,1);
  }
  else if (nsd > 2)
  {
    C(3,1) = C(2,1);
    C(3,2) = C(2,1);
    C(1,3) = C(2,1);
    C(2,3) = C(2,1);
    C(3,3) = C(1,1);
    C(5,5) = C(4,4);
    C(6,6) = C(4,4);
  }

  if (iop > 0)
  {
    // Calculate the stress tensor, sigma = C*eps
    Vector sig; // Use a local variable to avoid redimensioning of sigma
    if (eps.dim() != sigma.dim())
    {
      // Account for non-matching tensor dimensions
      SymmTensor epsil(sigma.dim(), nsd == 2 && axiSymmetry);
      if (!C.multiply(epsil=eps,sig))
	return false;
    }
    else
      if (!C.multiply(eps,sig))
	return false;

    sigma = sig; // Add sigma_zz in case of plane strain
    if (!planeStress && ! axiSymmetry && nsd == 2 && sigma.size() == 4)
      sigma(3,3) = nu * (sigma(1,1)+sigma(2,2));
  }

  if (iop == 3) // Calculate strain energy density, // U = 0.5*sigma:eps
    U = 0.5*sigma.innerProd(eps);

  return true;
}
