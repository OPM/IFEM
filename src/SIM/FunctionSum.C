// $Id$
//==============================================================================
//!
//! \file FunctionSum.C
//!
//! \date Apr 16 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Unary spatial function as a sum of other spatial functions.
//!
//==============================================================================

#include "FunctionSum.h"
#include "matrix.h"


FunctionSum::FunctionSum (FunctionBase* f, double w)
{
  comps.push_back(std::make_pair(f,w));
  ncmp = f->dim();
}


bool FunctionSum::add (FunctionBase* f, double w)
{
  if (f->dim() == comps.front().first->dim())
  {
    comps.push_back(std::make_pair(f,w));
    return true;
  }

  std::cerr <<" *** FunctionSum::add: Inconsistent dimensions "
            << f->dim() <<" != "<< comps.front().first->dim() << std::endl;
  return false;
}


bool FunctionSum::inDomain (const Vec3& X) const
{
  for (const WeightedFunc& cmp : comps)
    if (cmp.first->inDomain(X))
      return true;

  return false;
}


bool FunctionSum::initPatch (size_t idx)
{
  bool affected = false;
  for (WeightedFunc& cmp : comps)
    affected |= cmp.first->initPatch(idx);

  return affected;
}


std::vector<double> FunctionSum::getValue (const Vec3& X) const
{
  utl::vector<double> sum(ncmp);
  for (size_t i = 0; i < comps.size(); i++)
    if (comps[i].second > 0.0)
      sum.add(comps[i].first->getValue(X),comps[i].second);
    else if (i == 0)
      sum = comps[i].first->getValue(X);
    else
    {
      // Find the max value
      std::vector<double> val = comps[i].first->getValue(X);
      for (size_t j = 0; j < val.size(); j++)
        if (val[j] > sum[j]) sum[j] = val[j];
    }

  return sum;
}


double FunctionSum::getScalarValue (const Vec3& X) const
{
  double sum = 0.0;
  for (size_t i = 0; i < comps.size(); i++)
    if (comps[i].second > 0.0)
      sum += comps[i].first->getScalarValue(X)*comps[i].second;
    else
    {
      // Find the max value
      double val = comps[i].first->getScalarValue(X);
      if (i == 0 || val > sum) sum = val;
    }

  return sum;
}
