// $Id$
//==============================================================================
//!
//! \file BDF.C
//!
//! \date Oct 29 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Helper functions for BDF based time stepping.
//!
//==============================================================================

#include "BDF.h"


void TimeIntegration::BDF::setOrder (int order)
{
  if (order >= 1)
    coefs1.resize(2,-1.0);

  if (order <= 1)
    coefs = coefs1;
  else {
    coefs.resize(3);
    coefs[0] =  1.5;
    coefs[1] = -2.0;
    coefs[2] =  0.5;
  }
}


int TimeIntegration::BDF::getOrder() const
{
  int degree = this->getDegree();
  return step < degree+1 ? coefs1.size()-1 : coefs.size()-degree;
}


void TimeIntegration::BDF::advanceStep (double dt, double dtn)
{
  if (++step > 2 && coefs.size() == 3 && dt > 0.0 && dtn > 0.0)
  {
    double tau = dt/dtn;
    double taup1 = tau + 1.0;

    coefs[0] = (1.0+2.0*tau)/taup1;
    coefs[1] = -taup1;
    coefs[2] = tau*tau/taup1;
  }
}


const std::vector<double>& TimeIntegration::BDF::getCoefs () const
{
  return step < 2 ? coefs1 : coefs;
}


void TimeIntegration::BDFD2::setOrder (int order)
{
  if (order >= 1) {
    // Assume zero time derivative at t = 0
    coefs1.resize(2);
    coefs1[0] =  2.0;
    coefs1[1] = -2.0;
  }

  if (order >= 2) {
    // Use second order for second step
    coefs2.resize(3);
    coefs2[0] =  2.5;
    coefs2[1] = -8.0;
    coefs2[2] =  5.5;
  }

  coefs.resize(order <= 2 ? order+2 : 4, 1.0);
  if (order == 1) {
    coefs[0] =  1.0;
    coefs[1] = -2.0;
    coefs[2] =  1.0;
  }
  else if (order >= 2) {
    coefs[0] =  2.0;
    coefs[1] = -5.0;
    coefs[2] =  4.0;
    coefs[3] = -1.0;
  }
}


const std::vector<double>& TimeIntegration::BDFD2::getCoefs () const
{
  if (step < 2)
    return coefs1;
  else if (step == 2 && this->getActualOrder() >= 2)
    return coefs2;
  else
    return coefs;
}
