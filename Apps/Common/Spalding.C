//==============================================================================
//!
//! \file Spalding.C
//!
//! \date Jan 25 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Implementation of Spalding parametrization of a boundary layer.
//!
//==============================================================================

#include "Spalding.h"
#include <iostream>

bool Spalding::computeTauB(double hb, double CbI, double nu, double ut, double& tauB) const
{
  // Initial guess 
  tauB = CbI*nu/hb;

  double yplus = sqrt(ut/tauB);
  double uplus = yplus;
  double r     = yplus - f(uplus);

  double dr;
  int it = 0;
  while ((fabs(r) > rtol) && (it < maxit)) {
    dr = drdtauB(uplus,ut,nu,hb,CbI,tauB);
    tauB -= r/dr;
    yplus = hb/(nu*CbI)*sqrt(tauB*ut);
    uplus = sqrt(ut/tauB);
    r = yplus - f(uplus);
    it ++;
  }

  if ((it == maxit) && (fabs(r) > rtol)) {
    std::cout << "Spalding::computeTauB: Newton iteration did not converge: |r| = " 
	      << fabs(r) << std::endl;

    tauB = CbI*nu/hb;
    return false;
  }

  return true;
}
