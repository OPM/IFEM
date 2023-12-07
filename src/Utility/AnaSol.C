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
#include "tinyxml2.h"

#include <autodiff/reverse/var.hpp>

#ifdef HAS_HDF5
#include "FieldFunctions.h"
#endif


namespace {

//! \brief Function returning the gradient of a scalar or vector function.
template<class FromFunc, class ToFunc, class Ret>
class GradientFunc : public ToFunc
{
public:
  //! \brief The constructor initializes the function to use.
  GradientFunc(const FromFunc& f) : func(f) {}

protected:
  //! \brief Evaluates the gradient in a point.
  Ret evaluate(const Vec3& X) const { return func.gradient(X); }

private:
  const FromFunc& func; //!< Reference to scalar/vector function
};


//! \brief Specialization for symmetric tensors.
//! \details We need to construct a symmtensor from the returned tensor.
template<>
SymmTensor GradientFunc<VecFunc,STensorFunc,SymmTensor>::evaluate(const Vec3& X) const
{
  Tensor tmp = func.gradient(X);
  const size_t nsd = tmp.dim();
  SymmTensor ret(nsd);
  for (size_t i = 1; i <= nsd; ++i)
    for (size_t j = i; j <= nsd; ++j)
      ret(i,j) = tmp(i,j);

  return ret;
}

}


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
      IFEM::cout <<"\tVariables: "<< variables << std::endl;
    }

    if ((pos = function.find("Primary=")) != std::string::npos)
    {
      std::string primary = function.substr(pos+8);
      IFEM::cout <<"\tPrimary: "<< primary << std::endl;
      if (scalarSol)
        scalSol.push_back(new EvalFunction((variables+primary).c_str()));
      else
        vecSol = new VecFuncExpr(primary,variables);
    }

    if ((pos = function.find("Secondary=")) != std::string::npos)
    {
      std::string secondary = function.substr(pos+10);
      IFEM::cout <<"\tSecondary: "<< secondary << std::endl;
      if (scalarSol)
        scalSecSol.push_back(new VecFuncExpr(secondary,variables));
      else
        vecSecSol = new TensorFuncExpr(secondary,variables);
    }

    if ((pos = function.find("Stress=")) != std::string::npos)
    {
      std::string stress = function.substr(pos+7);
      IFEM::cout <<"\tStress: "<< stress << std::endl;
      stressSol = new STensorFuncExpr(stress,variables);
    }
  }
}


AnaSol::AnaSol (const tinyxml2::XMLElement* elem, bool scalarSol)
  : vecSol(nullptr), vecSecSol(nullptr), stressSol(nullptr)
{
  const char* type = elem->Attribute("type");
  if (type && !strcasecmp(type,"fields"))
    this->parseFieldFunctions(elem,scalarSol);
  else {
    bool useAD = false;
    utl::getAttribute(elem, "autodiff", useAD);
    utl::getAttribute(elem, "symmetric", symmetric);
    if (useAD)
      this->parseExpressionFunctions<autodiff::var>(elem,scalarSol);
    else
      this->parseExpressionFunctions<Real>(elem,scalarSol);
  }
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


void AnaSol::initPatch (size_t pIdx)
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


template<class Scalar>
void AnaSol::parseExpressionFunctions (const tinyxml2::XMLElement* elem, bool scalarSol)
{
  using EvalF = EvalFuncSpatial<Scalar>;
  using VecF = EvalMultiFunction<VecFunc,Vec3,Scalar>;
  using TensorF = EvalMultiFunction<TensorFunc,Tensor,Scalar>;
  using STensorF = EvalMultiFunction<STensorFunc,SymmTensor,Scalar>;
  std::string variables;
  const tinyxml2::XMLElement* var = elem->FirstChildElement("variables");
  if (var && var->FirstChild())
  {
    variables = var->FirstChild()->Value();
    if (variables[variables.size()-1] != ';') variables += ";";
    IFEM::cout <<"\tVariables: "<< variables << std::endl;
  }

  // Lambda function for parsing the derivative of an expression function
  auto&& parseDerivatives = [variables](auto* func, const tinyxml2::XMLElement* elem)
  {
    const tinyxml2::XMLElement* child = elem->FirstChildElement("derivative");
    for (; child; child = child->NextSiblingElement("derivative"))
    {
      int dir1 = 0, dir2 = 0;
      if (!utl::getAttribute(child,"dir",dir1))
        if (utl::getAttribute(child,"d1",dir1))
          utl::getAttribute(child,"d2",dir2);
      IFEM::cout <<"\tDerivative_"<< dir1;
      if (dir2 > 0) IFEM::cout << dir2;
      std::string derivative = child->FirstChild()->Value();
      IFEM::cout <<": "<< derivative << std::endl;
      func->addDerivative(derivative,variables,dir1,dir2);
    }
  };

  const tinyxml2::XMLElement* prim = elem->FirstChildElement("primary");
  if (prim && prim->FirstChild())
  {
    std::string type = "expression";
    utl::getAttribute(prim, "type", type);
    std::string prType = (type == "expression" ? "" : " ("+type+")");
    std::string primary = prim->FirstChild()->Value();
    IFEM::cout <<"\tPrimary" << prType << ": " << primary << std::endl;
    if (scalarSol)
    {
      if (type == "expression") {
        scalSol.push_back(new EvalF((variables+primary).c_str()));
        parseDerivatives(static_cast<EvalF*>(scalSol.back()),prim);
      } else
        scalSol.push_back(utl::parseRealFunc(primary,type,false));
    }
    else
    {
      if (type == "expression") {
        vecSol = new VecF(primary,variables);
        parseDerivatives(static_cast<VecF*>(vecSol),prim);
      } else
        vecSol = utl::parseVecFunc(primary, type);
    }
  }

  prim = elem->FirstChildElement("scalarprimary");
  while (prim && prim->FirstChild())
  {
    std::string type = "expression";
    utl::getAttribute(prim, "type", type);
    std::string prType = (type == "expression" ? "" : " ("+type+")");
    std::string primary = prim->FirstChild()->Value();
    IFEM::cout <<"\tScalar Primary" << prType << ": " << primary << std::endl;
    if (type == "expression")
      scalSol.push_back(new EvalF((variables+primary).c_str()));
    else
      scalSol.push_back(utl::parseRealFunc(primary, type, false));
    prim = prim->NextSiblingElement("scalarprimary");
  }

  const tinyxml2::XMLElement* sec = elem->FirstChildElement("secondary");
  if (sec && sec->FirstChild())
  {
    std::string type = "expression";
    utl::getAttribute(sec, "type", type);
    std::string prType = (type == "expression" ? "" : " ("+type+") ");
    std::string secondary = sec->FirstChild()->Value();
    IFEM::cout <<"\tSecondary" << prType << ": " << secondary << std::endl;
    if (scalarSol)
    {
      if (type == "expression") {
        scalSecSol.push_back(new VecF(secondary,variables));
        parseDerivatives(static_cast<VecF*>(scalSecSol.back()),sec);
      } else
        scalSecSol.push_back(utl::parseVecFunc(secondary, type));
    }
    else
    {
      if (type == "expression") {
        vecSecSol = new TensorF(secondary,variables);
        parseDerivatives(static_cast<TensorF*>(vecSecSol),sec);
      } else
        vecSecSol = utl::parseTensorFunc(secondary, type);
    }
  }

  sec = elem->FirstChildElement("scalarsecondary");
  while (sec && sec->FirstChild())
  {
    std::string type = "expression";
    utl::getAttribute(sec, "type", type);
    std::string secondary = sec->FirstChild()->Value();
    std::string prType = (type == "expression" ? "" : " ("+type+") ");
    IFEM::cout <<"\tScalar Secondary" << prType << ": " << secondary << std::endl;
    if (type == "expression")
      scalSecSol.push_back(new VecF(secondary,variables));
    else
      scalSecSol.push_back(utl::parseVecFunc(secondary, type));
    sec = sec->NextSiblingElement("scalarsecondary");
  }

  const tinyxml2::XMLElement* stress = elem->FirstChildElement("stress");
  if (stress && stress->FirstChild())
  {
    symmetric = true;
    std::string type = "expression";
    std::string sigma = stress->FirstChild()->Value();
    std::string prType = (type == "expression" ? "" : " ("+type+") ");
    IFEM::cout <<"\tStress" << prType << ": " << sigma << std::endl;
    stressSol = new STensorF(sigma,variables);
    parseDerivatives(static_cast<STensorF*>(stressSol),stress);
  }
}


void AnaSol::parseFieldFunctions (const tinyxml2::XMLElement* elem, bool scalarSol)
{
#ifdef HAS_HDF5
  const char* file = elem->Attribute("file");
  const char* basis = elem->Attribute("file_basis");
  if (!file || !basis)
  {
    std::cerr <<" *** No file or basis specified."<< std::endl;
    return;
  }

  int level = 0;
  utl::getAttribute(elem, "level", level);

  const tinyxml2::XMLElement* prim = elem->FirstChildElement("primary");
  if (prim && prim->FirstChild())
  {
    std::string primary = prim->FirstChild()->Value();
    IFEM::cout <<"\tPrimary: "<< primary << std::endl;
    if (scalarSol)
      scalSol.push_back(new FieldFunction(file, basis, primary, level));
    else
      vecSol = new VecFieldFunction(file, basis, primary, level);
  }
  prim = elem->FirstChildElement("scalarprimary");
  if (prim && prim->FirstChild())
  {
    std::string primary = prim->FirstChild()->Value();
    IFEM::cout <<"\tScalar Primary: "<< primary << std::endl;
    scalSol.push_back(new FieldFunction(file, basis, primary, level));
  }

  const tinyxml2::XMLElement* sec = elem->FirstChildElement("secondary");
  if (sec && sec->FirstChild())
  {
    std::string secondary = sec->FirstChild()->Value();
    IFEM::cout <<"\tSecondary: "<< secondary << std::endl;
    if (scalarSol)
      scalSecSol.push_back(new VecFieldFunction(file, basis, secondary, level));
    else
      vecSecSol = new TensorFieldFunction(file, basis, secondary, level);
  }
  sec = elem->FirstChildElement("scalarsecondary");
  if (sec && sec->FirstChild())
  {
    std::string secondary = sec->FirstChild()->Value();
    IFEM::cout <<"\tScalar Secondary: "<< secondary << std::endl;
    scalSecSol.push_back(new VecFieldFunction(file, basis, secondary, level));
  }
  sec = elem->FirstChildElement("stress");
  if (sec && sec->FirstChild())
  {
    std::string secondary = sec->FirstChild()->Value();
    IFEM::cout <<"\tStress: "<< secondary << std::endl;
    stressSol = new STensorFieldFunction(file, basis, secondary, level);
  }
#else
  std::cerr <<" *** AnaSol::parseFieldFunctions: Compiled without HDF5 support"
            <<", no fields read."<< std::endl;
#endif
}


void AnaSol::setupSecondarySolutions ()
{
  // if we are given a stress sol, no scalar secondaries should be registered.
  // this is tailored to assumptions in elasticity applications.
  if (!stressSol && !scalSol.empty()) {
    scalSecSol.resize(scalSol.size());
    for (size_t i = 0; i < scalSol.size(); ++i)
      if (!scalSecSol[i])
        scalSecSol[i] = new GradientFunc<RealFunc,VecFunc,Vec3>(*scalSol[i]);
  }

  if (vecSol) {
    if (symmetric && !stressSol)
      stressSol = new GradientFunc<VecFunc,STensorFunc,SymmTensor>(*vecSol);
    else if (!symmetric && !vecSecSol)
      vecSecSol = new GradientFunc<VecFunc,TensorFunc,Tensor>(*vecSol);
  }
}
