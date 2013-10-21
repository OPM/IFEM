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


TimeIntegration::BDF::BDF (int order) : step(0)
{
  coefs1.resize(order > 0 ? 2 : 1, 1.0);
  if (order > 0)
    coefs1[1] = -1.0;

  this->setOrder(order);
}


void TimeIntegration::BDF::setOrder (int order)
{
  if (order <= 1)
    coefs = coefs1;
  else if (order == 2) {
    coefs.resize(3);
    coefs[0] =  1.5;
    coefs[1] = -2.0;
    coefs[2] =  0.5;
  }
}


void TimeIntegration::BDF::advanceStep (double dt, double dtn)
{
  if (++step > 2 && coefs.size() == 3) {
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


double TimeIntegration::BDF::extrapolate (const double* values) const
{
  if (step > 1 && coefs.size() == 3) // second order
    return 2.0*values[0] - values[1];
  else // first order
    return values[0];
}
