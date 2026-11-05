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
#include "matrix.h"

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
    comps.emplace_back(f,w);
    ncmp = f->dim();
    return true;
  }
  else if (f->dim() == comps.front().first->dim())
  {
    comps.emplace_back(f,w);
    return true;
  }

  std::cerr <<" *** FunctionSum::add: Inconsistent dimensions "
            << f->dim() <<" != "<< comps.front().first->dim() << std::endl;
  if (ownFunc) delete f;
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


void FunctionSum::setParam (const std::string& name, double value)
{
  for (WeightedFunc& func : comps)
    func.first->setParam(name,value);
}


void RealFuncSum::addFuncComp (const char* ampl, RealFunc* f)
{
  if (ampl)
  {
    if (ScalarFunc* tf = utl::parseTimeFunc(ampl); tf)
      f = new SpaceTimeFunc(f,tf);
    else
      IFEM::cout <<"  ** Malformed time expression \""<< ampl <<"\","
                 <<" the function component will remain constant."<< std::endl;
  }

  this->add(f);
}


DiracSum::DiracSum (const char* input, double tol, int nsd)
{
  if (!input || input[0] == 0)
    return; // avoid segfault on empty string

  size_t nc = strlen(input);
  char* cpy = strdup(input);
  // Replace all '\' and '|' characters in the string by newline '\n'
  for (size_t i = 0; i < nc; i++)
    if (cpy[i] == '\\' || cpy[i] == '|')
      cpy[i] = '\n';

  IFEM::cout <<": DiracSum\n";
  std::stringstream str(cpy);
  char temp[512];
  while (str.getline(temp,512))
    if (temp[0] != '#' && temp[0] != 0)
    {
      std::stringstream sline(temp);
      Vec3 X;
      std::string value;
      for (int i = 0; i < nsd; i++)
        sline >> X[i];
      sline >> value;

      IFEM::cout <<"\t\tDirac("<< X.x;
      for (int i = 1; i < nsd; i++)
        IFEM::cout <<", "<< X[i];
      IFEM::cout <<") = ";

      char* endPtr = nullptr;
      const char* ampExpr = value.c_str();
      double ampConst = strtod(ampExpr,&endPtr);
      if (strlen(endPtr) == 0)
      {
        ampExpr = nullptr;
        IFEM::cout << ampConst << std::endl;
      }
      else
        ampConst = 1.0; // An expression was specified
      this->addFuncComp(ampExpr, new DiracSpaceFunc(ampConst,X,tol,nsd));
    }

  free(cpy);
}
