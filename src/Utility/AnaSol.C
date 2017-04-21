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
#include "ExprFunctions.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml.h"
#ifdef HAS_HDF5
#include "ProcessAdm.h"
#include "HDF5Writer.h"
#include "FieldFunctions.h"
#endif


AnaSol::AnaSol (RealFunc* s1, VecFunc* s2,
                VecFunc* v1, TensorFunc* v2, STensorFunc* v3)
  : vecSol(v1), vecSecSol(v2), stressSol(v3)
{
  if (s1) scalSol.push_back(s1);
  if (s2) scalSecSol.push_back(s2);
}


AnaSol::AnaSol (RealFunc* s, STensorFunc* sigma)
  : vecSol(nullptr), vecSecSol(nullptr), stressSol(sigma)
{
  if (s) scalSol.push_back(s);
}


AnaSol::AnaSol (std::istream& is, const int nlines, bool scalarSol)
  : vecSol(nullptr), vecSecSol(nullptr), stressSol(nullptr)
{
  std::string variables;
  for (int i = 0; i < nlines; i++)
  {
    std::string function = utl::readLine(is);

    size_t pos;
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


/*!
  \brief Static helper to parse the derivative of an expression function.
*/

template <class T> void parseDerivatives (T* func, const std::string& variables,
                                          const TiXmlElement* elem)
{
  const TiXmlElement* child = elem->FirstChildElement("derivative");
  for (; child; child = child->NextSiblingElement("derivative"))
  {
    int dir1 = 0, dir2 = 0;
    if (!utl::getAttribute(child,"dir",dir1))
      if (utl::getAttribute(child,"d1",dir1))
        utl::getAttribute(child,"d2",dir2);
    IFEM::cout <<"\tDerivative_"<< dir1;
    if (dir2 > 0) IFEM::cout << dir2;
    std::string derivative = child->FirstChild()->Value();
    IFEM::cout <<"="<< derivative << std::endl;
    func->addDerivative(derivative,variables,dir1,dir2);
  }
}


AnaSol::AnaSol (const TiXmlElement* elem, bool scalarSol)
  : vecSol(nullptr), vecSecSol(nullptr), stressSol(nullptr)
{
  const char* type = elem->Attribute("type");
  if (type && !strcasecmp(type,"fields"))
    parseFieldFunctions(elem, scalarSol);
  else
    parseExpressionFunctions(elem, scalarSol);
}


AnaSol::~AnaSol ()
{
  for (RealFunc* rf : scalSol)
    delete rf;
  for (VecFunc* vf : scalSecSol)
    delete vf;
  delete vecSol;
  delete vecSecSol;
  delete stressSol;
}


void AnaSol::initPatch(size_t pIdx)
{
  for (RealFunc* rf : scalSol)
    rf->initPatch(pIdx);

  for (VecFunc* rf : scalSecSol)
    rf->initPatch(pIdx);

  if (vecSol)
    vecSol->initPatch(pIdx);

  if (vecSecSol)
    vecSecSol->initPatch(pIdx);

  if (stressSol)
    stressSol->initPatch(pIdx);
}


void AnaSol::parseExpressionFunctions(const TiXmlElement* elem, bool scalarSol)
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
    {
      scalSol.push_back(new EvalFunction((variables+primary).c_str()));
      parseDerivatives(static_cast<EvalFunction*>(scalSol.back()),
                       variables,prim);
    }
    else
    {
      vecSol = new VecFuncExpr(primary,variables);
      parseDerivatives(static_cast<VecFuncExpr*>(vecSol),
                       variables,prim);
    }
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
    {
      scalSecSol.push_back(new VecFuncExpr(secondary,variables));
      parseDerivatives(static_cast<VecFuncExpr*>(scalSecSol.back()),
                       variables,sec);
    }
    else
    {
      vecSecSol = new TensorFuncExpr(secondary,variables);
      parseDerivatives(static_cast<TensorFuncExpr*>(vecSecSol),
                       variables,sec);
    }
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
    parseDerivatives(static_cast<STensorFuncExpr*>(stressSol),
                     variables,stress);
  }
}


void AnaSol::parseFieldFunctions(const TiXmlElement* elem, bool scalarSol)
{
#ifdef HAS_HDF5
  std::string file = elem->Attribute("file");
  int level = 0;
  utl::getAttribute(elem, "level", level);
  const TiXmlElement* prim = elem->FirstChildElement("primary");
  const char* basis = elem->Attribute("file_basis");
  if (!basis) {
    std::cerr << "No basis specified on file" << std::endl;
    return;
  }
  if (prim && prim->FirstChild())
  {
    std::string primary = prim->FirstChild()->Value();
    IFEM::cout <<"\tPrimary="<< primary << std::endl;
    if (scalarSol)
      scalSol.push_back(new FieldFunction(file, basis, primary, level));
    else
      vecSol = new VecFieldFunction(file, basis, primary, level);
  }
  prim = elem->FirstChildElement("scalarprimary");
  if (prim && prim->FirstChild())
  {
    std::string primary = prim->FirstChild()->Value();
    IFEM::cout <<"\tScalar Primary="<< primary << std::endl;
    scalSol.push_back(new FieldFunction(file, basis, primary, level));
  }
  const TiXmlElement* sec = elem->FirstChildElement("secondary");
  if (sec && sec->FirstChild())
  {
    std::string secondary = sec->FirstChild()->Value();
    IFEM::cout <<"\tSecondary="<< secondary << std::endl;
    // first component runs the show
    size_t pos = secondary.find_first_of('|');
    std::string sec;
    if (pos == std::string::npos)
      sec = secondary;
    else
      sec = secondary.substr(0,pos);
    if (scalarSol)
      scalSecSol.push_back(new VecFieldFunction(file, basis, secondary, level));
    else
      vecSecSol = new TensorFieldFunction(file, basis, secondary, level);
  }
  sec = elem->FirstChildElement("scalarsecondary");
  if (sec && sec->FirstChild())
  {
    std::string secondary = sec->FirstChild()->Value();
    IFEM::cout <<"\tScalar Secondary="<< secondary << std::endl;
    // first component runs the show
    size_t pos = secondary.find_first_of('|');
    std::string sec;
    if (pos == std::string::npos)
      sec = secondary;
    else
      sec = secondary.substr(0,pos);
    scalSecSol.push_back(new VecFieldFunction(file, basis, secondary, level));
  }
  sec = elem->FirstChildElement("stress");
  if (sec && sec->FirstChild())
  {
    std::string secondary = sec->FirstChild()->Value();
    IFEM::cout <<"\tStress="<< secondary << std::endl;
    // first component runs the show
    size_t pos = secondary.find_first_of('|');
    std::string sec;
    if (pos == std::string::npos)
      sec = secondary;
    else
      sec = secondary.substr(0,pos);
    stressSol = new STensorFieldFunction(file, basis, secondary, level);
  }
#else
  std::cerr << "AnaSol::parseFieldFunctions(..): Compiled without HDF5 support"
            <<", no fields read" << std::endl;
#endif
}
