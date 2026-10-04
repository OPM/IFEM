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

  ++step;
  tau = dt > 0.0 && dtn > 0.0 ? dt/dtn : 1.0;
  const double dtnp = dtnn;
  dtnn = dtn;
  if (dt <= 0.0 || dtn <= 0.0)
    return true; // constant step size

  // The variable step coefficients replace those of a constant step for the
  // steps using all the points of the scheme, that is, from the second step
  // with a first derivative and from the third with a second one, the second
  // step of the second order scheme for second derivatives being special
  if (degree == 1 && coefs.size() == 3 && step >= 2)
    this->variableCoefs(coefs,{dt,dtn});
  else if (degree == 2 && coefs.size() == 3 && step >= 2)
    this->variableCoefs(coefs,{dt,dtn});
  else if (degree == 2 && coefs.size() == 4 && step >= 3 && dtnp > 0.0)
    this->variableCoefs(coefs,{dt,dtn,dtnp});

  return true;
}


/*!
  The points are at the times \f$t_0 = 0\f$ of the new step and
  \f$t_k = t_{k-1} - \Delta t_{k-1}\f$ before it, and the coefficient of
  point \a k is the derivative of its Lagrange polynomial at \f$t_0\f$.
*/

void TimeIntegration::BDF::variableCoefs (std::vector<double>& c,
                                          const std::vector<double>& dts) const
{
  const size_t n = dts.size() + 1;
  std::vector<double> t(n,0.0);
  for (size_t k = 1; k < n; k++)
    t[k] = t[k-1] - dts[k-1];

  // The product of (0 - t_m)/(t_k - t_m) over the points m other than k, j
  // and l
  auto product = [&t,n](size_t k, size_t j, size_t l)
  {
    double p = 1.0;
    for (size_t m = 0; m < n; m++)
      if (m != k && m != j && m != l)
        p *= -t[m]/(t[k]-t[m]);
    return p;
  };

  c.assign(n,0.0);
  for (size_t k = 0; k < n; k++)
    for (size_t j = 0; j < n; j++)
      if (j != k) {
        if (degree == 1)
          c[k] += product(k,j,j)/(t[k]-t[j]);
        else
          for (size_t l = 0; l < n; l++)
            if (l != k && l != j)
              c[k] += product(k,j,l)/((t[k]-t[j])*(t[k]-t[l]));
      }

  const double scale = degree == 1 ? dts.front() : dts.front()*dts.front();
  for (double& ck : c)
    ck *= scale;
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
