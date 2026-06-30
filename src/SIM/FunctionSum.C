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
#include "Functions.h"
#include "IFEM.h"

#include <sstream>
#include <cstring>


FunctionSum::~FunctionSum ()
{
  if (ownFunc)
    for (WeightedFunc& fn : comps)
      delete fn.first;
}


bool FunctionSum::add (FunctionBase* f, double w)
{
  if (comps.empty())
  {
    comps.push_back(std::make_pair(f,w));
    ncmp = f->dim();
    return true;
  }
  else if (f->dim() == comps.front().first->dim())
  {
    comps.push_back(std::make_pair(f,w));
    return true;
  }

  std::cerr <<" *** FunctionSum::add: Inconsistent dimensions "
            << f->dim() <<" != "<< comps.front().first->dim() << std::endl;
  return false;
}


unsigned char FunctionSum::getType () const
{
  if (comps.empty()) return 0;

  unsigned char myType = comps.front().first->getType();
  for (const WeightedFunc& cmp : comps)
    if (cmp.first->getType() != myType)
      return 0;

  return myType;
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


DiracSum::DiracSum (const char* input, double tol, int nsd)
  : FunctionSum(true)
{
  if (!input || input[0] == 0)
    return; // avoid segfault on empty string

  size_t nc = strlen(input);
  char* cpy = strdup(input);
  // Replace all '\' and '|' characters in the string by newline '\n'
  for (size_t i = 0; i < nc; i++)
    if (cpy[i] == '\\' || cpy[i] == '|')
      cpy[i] = '\n';

  IFEM::cout <<" DiracSum";
  std::stringstream str(cpy);
  char temp[512];
  while (str.getline(temp,512))
    if (temp[0] != '#' && temp[0] != 0)
    {
      std::stringstream sline(temp);
      Vec3 X;
      double value = 0.0;
      for (int i = 0; i < nsd; i++)
        sline >> X[i];
      sline >> value;

      IFEM::cout <<"\n\t\tDirac("<< X.x;
      for (int i = 1; i < nsd; i++)
        IFEM::cout <<", "<< X[i];
      IFEM::cout <<") = "<< value;

      this->add(new DiracSpaceFunc(value,X,tol,nsd));
    }
  IFEM::cout << std::endl;

  free(cpy);
}
