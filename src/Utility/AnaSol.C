// $Id$
//==============================================================================
//!
//! \file AnaSol.C
//!
//! \date Feb 20 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Analytical solution fields (primary and secondary).
//!
//==============================================================================

#include "AnaSol.h"
#include "Functions.h"
#include "Utilities.h"
#include "tinyxml.h"


AnaSol::AnaSol (std::istream& is, const int nlines, bool scalarSol)
  : scalSol(0), scalSecSol(0), vecSol(0), vecSecSol(0), stressSol(0)
{
  size_t pos = 0;
  std::string variables;
  for (int i = 0; i < nlines; i++)
  {
    std::string function = utl::readLine(is);

    if ((pos = function.find("Variables=")) != std::string::npos)
    {
      variables += function.substr(pos+10);
      if (variables[variables.size()-1] != ';') variables += ";";
      std::cout <<"\tVariables="<< variables << std::endl;
    }

    if ((pos = function.find("Primary=")) != std::string::npos)
    {
      std::string primary = function.substr(pos+8);
      std::cout <<"\tPrimary="<< primary << std::endl;
      if (scalarSol)
        scalSol = new EvalFunction((variables+primary).c_str());
      else
        vecSol = new VecFuncExpr(primary,variables);
    }

    if ((pos = function.find("Secondary=")) != std::string::npos)
    {
      std::string secondary = function.substr(pos+10);
      std::cout <<"\tSecondary="<< secondary << std::endl;
      if (scalarSol)
        scalSecSol = new VecFuncExpr(secondary,variables);
      else
        vecSecSol = new TensorFuncExpr(secondary,variables);
    }

    if ((pos = function.find("Stress=")) != std::string::npos)
    {
      std::string stress = function.substr(pos+7);
      std::cout <<"\tStress="<< stress << std::endl;
      stressSol = new STensorFuncExpr(stress,variables);
    }
  }
}


AnaSol::AnaSol (const TiXmlElement* elem, bool scalarSol)
  : scalSol(0), scalSecSol(0), vecSol(0), vecSecSol(0), stressSol(0)
{
  std::string variables;

  const TiXmlElement* var = elem->FirstChildElement("variables");
  if (var && var->FirstChild())
  {
    variables = var->FirstChild()->Value();
    if (variables[variables.size()-1] != ';') variables += ";";
    std::cout <<"\tVariables="<< variables << std::endl;
  }

  const TiXmlElement* prim = elem->FirstChildElement("primary");
  if (prim && prim->FirstChild())
  {
    std::string primary = prim->FirstChild()->Value();
    std::cout <<"\tPrimary="<< primary << std::endl;
    if (scalarSol)
      scalSol = new EvalFunction((variables+primary).c_str());
    else
      vecSol = new VecFuncExpr(primary,variables);
  }

  prim = elem->FirstChildElement("scalarprimary");
  if (prim && prim->FirstChild())
  {
    std::string primary = prim->FirstChild()->Value();
    std::cout <<"\tScalar Primary="<< primary << std::endl;
    scalSol = new EvalFunction((variables+primary).c_str());
  }

  const TiXmlElement* sec = elem->FirstChildElement("secondary");
  if (sec && sec->FirstChild())
  {
    std::string secondary = sec->FirstChild()->Value();
    std::cout <<"\tSecondary="<< secondary << std::endl;
    if (scalarSol)
      scalSecSol = new VecFuncExpr(secondary,variables);
    else
      vecSecSol = new TensorFuncExpr(secondary,variables);
  }

  sec = elem->FirstChildElement("scalarsecondary");
  if (sec && sec->FirstChild())
  {
    std::string secondary = sec->FirstChild()->Value();
    std::cout <<"\tScalar Secondary="<< secondary << std::endl;
    scalSecSol = new VecFuncExpr(secondary,variables);
  }

  const TiXmlElement* stress = elem->FirstChildElement("stress");
  if (stress && stress->FirstChild())
  {
    std::string sigma = stress->FirstChild()->Value();
    std::cout <<"\tStress="<< sigma << std::endl;
    stressSol = new STensorFuncExpr(sigma,variables);
  }
}
