// $Id: AnalyticSolutions.C,v 1.11 2011-02-08 12:19:52 rho Exp $
//==============================================================================
//!
//! \file AnalyticSolutions.C
//!
//! \date Jul 1 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Analytic solutions for some linear elasticity and Poisson problems.
//!
//==============================================================================

#include "AnalyticSolutions.h"
#include "Vec3.h"
#include "Vec3Oper.h"
#include <math.h>


/*!
  \class Hole

  Smooth plane strain problem.

  Reference: Zienkiewicz & Zhu, "The superconvergent patch recovery...",
  IJNME, 33, 1331-1364, 1992, page 1355 (eq. 36).
*/

Tensor Hole::evaluate (const Vec3& X) const
{
  double R  = hypot(X.x,X.y);
  double th = atan2(X.y,X.x);
  double C2 = cos(2.0*th);
  double C4 = cos(4.0*th);
  double S2 = sin(2.0*th);
  double S4 = sin(4.0*th);
  double R2 = R <= a ? 1.0 : a*a/(R*R);
  double R4 = R <= a ? 1.0 : R2*R2;

  SymmTensor sigma(is3D ? 3 : 2);
  sigma(1,1) = F0 * (1.0 - R2*(1.5*C2 + C4) + 1.5*R4*C4);
  sigma(2,2) = F0 * (    - R2*(0.5*C2 - C4) - 1.5*R4*C4);
  sigma(1,2) = F0 * (    - R2*(0.5*S2 + S4) + 1.5*R4*S4);
  if (is3D) sigma(3,3) = F0 * nu*(1.0 - 2.0*R2*C2);

  return sigma;
}


/*!
  \class Lshape

  Symmetric (mode 1) loading. Plane strain problem with a singular point.

  Reference: Szabo & Babuska, "Finite Element Analysis", 1991, page 188-192.

  \code
                   y ^
                     |

     (-a,a) +--------+(0,a)
            |        |
            |        |
            |        |
            |        +--------+ (a,0) - - -> x
            |   (0,0)         |
            |                 |
            |                 |
    (-a,-a) +-----------------+ (a,-a)

   \endcode
*/

Lshape::Lshape (double r, double f, double P, bool use3D)
  : a(r), F0(f), nu(P), is3D(use3D), T(use3D ? 3 : 2)
{
  // Set up the local-to-global transformation tensor
  T(1,1) = T(2,2) = T(2,1) = -sqrt(0.5);
  T(1,2) = sqrt(0.5);

  if (is3D) T(3,3) = 1.0;
}


Tensor Lshape::evaluate (const Vec3& X) const
{
  // Some constants (see Szabo & Babuska for an elaboration)
  const double lambda = 0.544483737;
  const double q      = 0.543075579;
  const double lp1    = lambda + 1.0;
  const double lm1    = lambda - 1.0;
  const double lm3    = lambda - 3.0;
  const double tol    = a*1.0e-32;

  // Find local (polar) coordinates
  double x = X.x*T(1,1) + X.y*T(2,1);
  double y = X.x*T(1,2) + X.y*T(2,2);
  double r = hypot(x,y);
  if (r < tol) r = tol; // truncate the singularity to avoid NaN values
  double theta = atan2(y,x);

  // Set up the stress tensor in local system
  SymmTensor sigma(is3D ? 3 : 2);
  double c0  = F0*lambda*pow(r,lm1);
  sigma(1,1) = c0*((2.0-q*lp1)*cos(lm1*theta) - lm1*cos(lm3*theta));
  sigma(2,2) = c0*((2.0+q*lp1)*cos(lm1*theta) + lm1*cos(lm3*theta));
  sigma(1,2) = c0*(     q*lp1 *sin(lm1*theta) + lm1*sin(lm3*theta));
  if (is3D) sigma(3,3) = nu * (sigma(1,1)+sigma(2,2));

  // Transform to global coordinates
  return sigma.transform(T);
}


Tensor CanTS::evaluate (const Vec3& X) const
{
  double x = X.x/L;
  double y = (is3D ? X.z : X.y)/H - 0.5;
  double I = H*H*H / 12.0;
  size_t n = is3D ? 3 : 2;

  SymmTensor sigma(n);
  sigma(1,1) = F0*L*H/I * (x-1.0)*y;
  sigma(1,n) = F0*H*H/I * 0.5*(0.25-y*y);

  return sigma;
}


Tensor CanTM::evaluate (const Vec3& X) const
{
  double y = (is3D ? X.z : X.y)/H - 0.5;
  double I = H*H*H / 12.0;

  SymmTensor sigma(is3D ? 3 : 2);
  sigma(1,1) = M0*H/I * y;

  return sigma;
}


/*!
  \class CurvedBeam

  Smooth plane stress problem with a prescribed edge displacement.

  Reference: Zienkiewicz & Taylor, "The Finite Element Method",
  pages 42-45.
*/

CurvedBeam::CurvedBeam (double u0, double Ri, double Ro, double E, bool use3D)
  : a(Ri), b(Ro), is3D(use3D), T(use3D ? 3 : 2)
{
  if (is3D) T(3,3) = 1.0;

  PN = -u0*E/(M_PI*(a*a+b*b));
}


Tensor CurvedBeam::evaluate (const Vec3& X) const
{
  // Find polar coordinates
  double r     = hypot(X.x,X.y);
  double theta = atan2(X.y,X.x);

  // Set up the stress tensor in polar coordinates
  SymmTensor sigma(is3D ? 3 : 2);
  double c1 = a*a*b*b/(r*r*r);
  double c2 = (a*a + b*b)/r;
  double ct = cos(theta);
  double st = sin(theta);
  sigma(1,1) = PN*(    r + c1 - c2)*st;
  sigma(2,2) = PN*(3.0*r - c1 - c2)*st;
  sigma(1,2) = PN*(   -r - c1 + c2)*ct;

  // Transform to Cartesian coordinates
  T(1,1) =  ct;
  T(2,1) =  st;
  T(1,2) = -st;
  T(2,2) =  ct;
  return sigma.transform(T);
}


Vec3 Square2D::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;
  double pi = M_PI;

  Vec3 temp;
  temp.x = -pi*sin(pi*x)*(2.0-y);
  temp.y = -cos(pi*x);

  return temp;
}


double Square2DHeat::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;
  double pi = M_PI;

  return -pi*pi*cos(pi*x)*(2.0-y);
}


Vec3 LshapePoisson::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;
  double R = hypot(x,y);
  double pi = M_PI;
  double Rp = pow(R,-1.0/3.0);
  double th = x > 0.0 ? asin(y/R) : pi-asin(y/R);
  double frac = 2.0/3.0;

  Vec3 temp;
  temp.x =  frac*Rp*sin(th/3.0);
  temp.y = -frac*Rp*cos(th/3.0);

  return temp;
}


Vec3 PoissonCube::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;
  double z = X.z;
  double pi = M_PI;

  Vec3 temp;
  temp.x = -pi*cos(x*pi)*sin(y*pi)*sin(z*pi);
  temp.y = -pi*sin(x*pi)*cos(y*pi)*sin(z*pi);
  temp.z = -pi*sin(x*pi)*sin(y*pi)*cos(z*pi);

  return temp;
}


double PoissonCubeSource::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;
  double z = X.z;
  double pi = M_PI;

  return 3.0*pi*pi*sin(x*pi)*sin(y*pi)*sin(z*pi);
}


Vec3 PoissonLine::evaluate (const Vec3& X) const
{
  double x = X.x;
  double pi = M_PI;

  Vec3 temp;
  temp.x = -(pi/L)*cos(x*pi/L);

  return temp;
}


double PoissonLineSource::evaluate (const Vec3& X) const
{
  double x = X.x;
  double pi = M_PI;

  return (pi*pi)/(L*L)*sin(pi*x/L);
}
