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


TimeIntegration::BDF::BDF (int order) : degree(1), step(0)
{
  if (order >= 0)
    this->setOrder(order);
}


void TimeIntegration::BDF::setOrder (int order)
{
  if (degree > 1)
    return; // Not for 2nd order problems

  if (order < 1)
    coefs1 = { 1.0 };
  else
    coefs1 = { 1.0, -1.0 };

  if (order < 2)
    coefs = coefs1;
  else
    coefs = { 1.5, -2.0, 0.5 };
}


int TimeIntegration::BDF::getOrder() const
{
  return step < degree+1 ? coefs1.size()-1 : coefs.size()-degree;
}


bool TimeIntegration::BDF::advanceStep (double dt, double dtn)
{
  if (coefs1.size() < 2)
    return false; // stationary problem

  if (++step > 2 && degree == 1 && coefs.size() == 3 && dt > 0.0 && dtn > 0.0)
  {
    double tau = dt/dtn;
    double taup1 = tau + 1.0;

    coefs[0] = (1.0+2.0*tau)/taup1;
    coefs[1] = -taup1;
    coefs[2] = tau*tau/taup1;
  }

  return true;
}


const std::vector<double>& TimeIntegration::BDF::getCoefs () const
{
  return step < 2 ? coefs1 : coefs;
}


TimeIntegration::BDFD2::BDFD2 (int order, int step_) : BDF(-1)
{
  degree = 2;
  step = step_;

  if (order < 1)
    coefs1 = { 1.0 };
  else
    coefs1 = { 2.0, -2.0 }; // Assume zero time derivative at t = 0

  if (order > 1)
    coefs2 = { 2.5, -8.0, 5.5 }; // Use second order for second step

  coefs.resize(order <= 2 ? order+2 : 4, 1.0);
  if (order == 1)
    coefs[1] = -2.0;
  else if (order >= 2)
    coefs = { 2.0, -5.0,  4.0, -1.0 };
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
