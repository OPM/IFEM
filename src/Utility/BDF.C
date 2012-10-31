#include "BDF.h"

BDF::BDF(int order)
{
  coefs1.resize(2);
  coefs1[0] = 1.0; coefs1[1] = -1.0;
  if (order == 1)
    coefs = coefs1;
  else if (order == 2) {
    coefs.resize(3);
    coefs[0] = 3.0/2.0, coefs[1] = -2.0; coefs[2] = 1.0/2.0;
  }

  step = 0;
}

const std::vector<double>& BDF::getCoefs() const
{
  return step<2?coefs1:coefs;
}

double BDF::extrapolate(const double* values) const
{
  double result = values[0];
  if (coefs.size() == 3 && step > 1) { // second order
    result *= 2.0;
    result -= values[1];
  }

  return result;
}
