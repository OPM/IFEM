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
  for (const WeightedFunc& cmp : comps)
    sum.add(cmp.first->getValue(X),cmp.second);

  return sum;
}
