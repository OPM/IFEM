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
#include "IFEM.h"
#include "Utilities.h"
#include "tinyxml.h"


AnaSol::AnaSol(RealFunc* s1, VecFunc* s2, VecFunc* v1,
               TensorFunc* v2, STensorFunc* v3) :
  vecSol(v1), vecSecSol(v2), stressSol(v3)
{
  if (s1) scalSol.push_back(s1);
  if (s2) scalSecSol.push_back(s2);
}


AnaSol::AnaSol (std::istream& is, const int nlines, bool scalarSol)
  : vecSol(0), vecSecSol(0), stressSol(0)
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
      IFEM::cout <<"\tVariables="<< variables << std::endl;
    }

    if ((pos = function.find("Primary=")) != std::string::npos)
    {
      std::string primary = function.substr(pos+8);
      IFEM::cout <<"\tPrimary="<< primary << std::endl;
      if (scalarSol)
        scalSol.push_back(new EvalFunction((variables+primary).c_str()));
      else
        vecSol = new VecFuncExpr(primary,variables);
    }

    if ((pos = function.find("Secondary=")) != std::string::npos)
    {
      std::string secondary = function.substr(pos+10);
      IFEM::cout <<"\tSecondary="<< secondary << std::endl;
      if (scalarSol)
        scalSecSol.push_back(new VecFuncExpr(secondary,variables));
      else
        vecSecSol = new TensorFuncExpr(secondary,variables);
    }

    if ((pos = function.find("Stress=")) != std::string::npos)
    {
      std::string stress = function.substr(pos+7);
      IFEM::cout <<"\tStress="<< stress << std::endl;
      stressSol = new STensorFuncExpr(stress,variables);
    }
  }
}


AnaSol::AnaSol (const TiXmlElement* elem, bool scalarSol)
  : vecSol(0), vecSecSol(0), stressSol(0)
{
  std::string variables;

  const TiXmlElement* var = elem->FirstChildElement("variables");
  if (var && var->FirstChild())
  {
    variables = var->FirstChild()->Value();
    if (variables[variables.size()-1] != ';') variables += ";";
    IFEM::cout <<"\tVariables="<< variables << std::endl;
  }

  const TiXmlElement* prim = elem->FirstChildElement("primary");
  if (prim && prim->FirstChild())
  {
    std::string primary = prim->FirstChild()->Value();
    IFEM::cout <<"\tPrimary="<< primary << std::endl;
    if (scalarSol)
      scalSol.push_back(new EvalFunction((variables+primary).c_str()));
    else
      vecSol = new VecFuncExpr(primary,variables);
  }

  prim = elem->FirstChildElement("scalarprimary");
  while (prim && prim->FirstChild())
  {
    std::string primary = prim->FirstChild()->Value();
    IFEM::cout <<"\tScalar Primary="<< primary << std::endl;
    scalSol.push_back(new EvalFunction((variables+primary).c_str()));
    prim = prim->NextSiblingElement("scalarprimary");
  }

  const TiXmlElement* sec = elem->FirstChildElement("secondary");
  if (sec && sec->FirstChild())
  {
    std::string secondary = sec->FirstChild()->Value();
    IFEM::cout <<"\tSecondary="<< secondary << std::endl;
    if (scalarSol)
      scalSecSol.push_back(new VecFuncExpr(secondary,variables));
    else
      vecSecSol = new TensorFuncExpr(secondary,variables);
  }

  sec = elem->FirstChildElement("scalarsecondary");
  while (sec && sec->FirstChild())
  {
    std::string secondary = sec->FirstChild()->Value();
    IFEM::cout <<"\tScalar Secondary="<< secondary << std::endl;
    scalSecSol.push_back(new VecFuncExpr(secondary,variables));
    sec = sec->NextSiblingElement("scalarsecondary");
  }

  const TiXmlElement* stress = elem->FirstChildElement("stress");
  if (stress && stress->FirstChild())
  {
    std::string sigma = stress->FirstChild()->Value();
    IFEM::cout <<"\tStress="<< sigma << std::endl;
    stressSol = new STensorFuncExpr(sigma,variables);
  }
}
