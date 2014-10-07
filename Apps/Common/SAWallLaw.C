//==============================================================================
//!
//! \file SAWallLaw.C
//!
//! \date Jun 19 2013
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Implementation of Spalding parametrization of a turbulent
//!        boundary layer. Computes the mean velocity uplus parallel
//!        to a solid wall given the distance y and the tangential
//!        velocity component ut.
//!
//==============================================================================

#include "SAWallLaw.h"

SAWallLaw::f_SA_wall_law::f_SA_wall_law() : B(5.03339087905055799), a1(8.148221580024245),
					    a2(-6.9287093849022945), b1(7.4600876082527945), b2(7.468145790401841),
					    c1(2.5496773539754747), c2(1.3301651588535228), c3(3.599459109332379),
					    c4(3.6397531868684494)
{}


bool SAWallLaw::computeYplus(double y, double nu, double utan, double& yplus) const
{
  // Initial guess
  yplus = y;

  // yplus is zero for zero tangential velocity
  if (utan < 1.0e-12)
    return true;
  
  // Residual
  double uplus = y*utan/(nu*yplus); 
  double r = uplus - SA.up(yplus);
  
  int it = 0;
  double drdyplus;
  while ((fabs(r) > rtol) && (it < maxit)) {
    drdyplus = -y*utan/(nu*yplus*yplus) - SA.dupdyp(yplus);
    yplus -= r/drdyplus;
    yplus = fabs(yplus);
    uplus = y*utan/(nu*yplus); 
    r = uplus - SA.up(yplus);
    it++;
  }

  return (it < maxit);
}


bool SAWallLaw::computeUstar(double y, double nu, double utan, double& ustar) const
{
  double yplus;
  if (!this->computeYplus(y,nu,utan,yplus))
    return false;

  double uplus = SA.up(yplus);
  ustar = utan/uplus;
  return true;
}


bool SAWallLaw::computeTauB(double y, double nu, double utan, double& tauB) const
{
  // Constant
  const double CbI = 4.0;

  // Initial guess 
  tauB = CbI*nu/y;

  if (utan < 1.0e-12)
    return true;

  double yplus;
  if (!this->computeYplus(y,nu,utan,yplus))
    return false;

  double uplus = SA.up(yplus);
  double ustar = utan/uplus;
  tauB = ustar*ustar/utan;

  return true;
}


bool SAWallLaw::computeNu(double y, double nu, double utan, double& nuwall) const
{  
  double yplus;
  if (!this->computeYplus(y,nu,utan,yplus))
    return false;

  //! Wall BC for Spalart-Allmaras
  double uplus = SA.up(yplus);
  double ustar = utan/uplus;;
  nuwall = K*ustar*y;
  
  return true;
}
