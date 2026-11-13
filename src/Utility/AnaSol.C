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
#include <cstring>
#include "FieldFunctions.h"
#endif


namespace
{
  //! \brief Function returning the gradient of a scalar or vector function.
  template<class FromFunc, class ToFunc, class Ret>
  class GradientFunc : public ToFunc
  {
  public:
    //! \brief The constructor initializes the reference to the target function.
    explicit GradientFunc(const FromFunc& f) : func(f) {}

  protected:
    //! \brief Evaluates the gradient in a point.
    Ret evaluate(const Vec3& X) const override { return func.gradient(X); }

  private:
    const FromFunc& func; //!< Reference to a scalar or vector function
  };


  //! \brief Specialization for symmetric tensors.
  template<> SymmTensor
  GradientFunc<VecFunc,STensorFunc,SymmTensor>::evaluate(const Vec3& X) const
  {
    // Construct a SymmTensor from the returned Tensor assuming it is symmetric
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

    if (size_t pos = function.find("Variables="); pos != std::string::npos)
    {
      variables += function.substr(pos+10);
      if (variables[variables.size()-1] != ';') variables += ";";
      IFEM::cout <<"\tVariables: "<< variables << std::endl;
    }

    if (size_t pos = function.find("Primary="); pos != std::string::npos)
    {
      std::string primary = function.substr(pos+8);
      IFEM::cout <<"\tPrimary: "<< primary << std::endl;
      if (scalarSol)
        scalSol.push_back(new EvalFunction((variables+primary).c_str()));
      else
        vecSol = new VecFuncExpr(primary,variables);
    }

    if (size_t pos = function.find("Secondary="); pos != std::string::npos)
    {
      std::string secondary = function.substr(pos+10);
      IFEM::cout <<"\tSecondary: "<< secondary << std::endl;
      if (scalarSol)
        scalSecSol.push_back(new VecFuncExpr(secondary,variables));
      else
        vecSecSol = new TensorFuncExpr(secondary,variables);
    }

    if (size_t pos = function.find("Stress="); pos != std::string::npos)
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
  if (std::string tp; utl::getAttribute(elem,"type",tp,true) && tp == "fields")
    this->parseFieldFunctions(elem,scalarSol);
  else if (bool ad = false; utl::getAttribute(elem,"autodiff",ad) && ad)
    this->parseExpressionFunctions<autodiff::var>(elem,scalarSol);
  else
    this->parseExpressionFunctions<Real>(elem,scalarSol);
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
void AnaSol::parseExpressionFunctions (const tinyxml2::XMLElement* elem,
                                       bool scalarSol)
{
  using SclF     = EvalFuncSpatial<Scalar>;
  using VecF     = EvalMultiFunction<VecFunc,Vec3,Scalar>;
  using TensorF  = EvalMultiFunction<TensorFunc,Tensor,Scalar>;
  using STensorF = EvalMultiFunction<STensorFunc,SymmTensor,Scalar>;

  std::string variables;
  const tinyxml2::XMLElement* var = elem->FirstChildElement("variables");
  if (var && var->FirstChild())
  {
    variables = var->FirstChild()->Value();
    if (variables[variables.size()-1] != ';') variables += ";";
    IFEM::cout <<"\tVariables: "<< variables << std::endl;
  }

  std::string atype("expression");
  utl::getAttribute(elem,"type",atype);

  // Lambda function for parsing the derivative of an expression function
  auto&& parseDerivs = [&variables](auto* f, const tinyxml2::XMLElement* xml)
  {
    for (const tinyxml2::XMLElement* der = xml->FirstChildElement("derivative");
         der && der->FirstChild(); der = der->NextSiblingElement("derivative"))
    {
      const char* derivative = der->FirstChild()->Value();

      int dir1 = 0, dir2 = 0;
      if (!utl::getAttribute(der,"dir",dir1))
        if (utl::getAttribute(der,"d1",dir1))
          utl::getAttribute(der,"d2",dir2);
      IFEM::cout <<"\tDerivative_"<< dir1;
      if (dir2 > 0) IFEM::cout << dir2;
      IFEM::cout <<": "<< derivative << std::endl;
      f->addDerivative(derivative,variables,dir1,dir2);
    }
  };

  // Lambda function for parsing the primary solution.
  auto&& parsePrimary = [&variables,&atype,parseDerivs,this]
    (const tinyxml2::XMLElement* elem, bool scalar,
     const char* keyw = "primary")
  {
    for (const tinyxml2::XMLElement* prim = elem->FirstChildElement(keyw);
         prim && prim->FirstChild(); prim = prim->NextSiblingElement(keyw))
    {
      const char* primary = prim->FirstChild()->Value();

      std::string type(atype);
      utl::getAttribute(prim,"type",type);
      IFEM::cout <<"\tPrimary";
      if (type != "expression")
        IFEM::cout <<" ("<< type <<")";
      if (scalar)
      {
        IFEM::cout <<": "<< primary;
        if (type == "expression")
        {
          scalSol.push_back(new SclF((variables+primary).c_str()));
          parseDerivs(static_cast<SclF*>(scalSol.back()),prim);
        }
        else
          scalSol.push_back(utl::parseRealFunc(primary,type,false));
      }
      else if (!vecSol)
      {
        if (type == "expression")
        {
          IFEM::cout <<": "<< primary;
          vecSol = new VecF(primary,variables);
          parseDerivs(static_cast<VecF*>(vecSol),prim);
        }
        else
          vecSol = utl::parseVecFunc(primary,type);
      }
      IFEM::cout << std::endl;
    }
  };

  // Lambda function for parsing the secondary solution.
  auto&& parseSecondary = [&variables,&atype,parseDerivs,this]
    (const tinyxml2::XMLElement* elem, bool scalar,
     const char* keyw = "secondary")
  {
    for (const tinyxml2::XMLElement* secs = elem->FirstChildElement(keyw);
         secs && secs->FirstChild(); secs = secs->NextSiblingElement(keyw))
    {
      const char* secondary = secs->FirstChild()->Value();

      std::string type(atype);
      utl::getAttribute(secs,"type",type);
      IFEM::cout <<"\tSecondary";
      if (type != "expression")
        IFEM::cout <<" ("<< type <<")";
      if (scalar)
      {
        if (type == "expression")
        {
          IFEM::cout <<": "<< secondary;
          scalSecSol.push_back(new VecF(secondary,variables));
          parseDerivs(static_cast<VecF*>(scalSecSol.back()),secs);
        }
        else
          scalSecSol.push_back(utl::parseVecFunc(secondary,type));
      }
      else if (!vecSecSol)
      {
        IFEM::cout <<": "<< secondary;
        if (type == "expression")
        {
          vecSecSol = new TensorF(secondary,variables);
          parseDerivs(static_cast<TensorF*>(vecSecSol),secs);
        }
        else
          vecSecSol = utl::parseTensorFunc(secondary,type);
      }
      IFEM::cout << std::endl;
    }
  };

  parsePrimary(elem,scalarSol,"primary");
  parsePrimary(elem,true,"scalarprimary");

  parseSecondary(elem,scalarSol,"secondary");
  parseSecondary(elem,true,"scalarsecondary");

  if (const tinyxml2::XMLElement* stress = elem->FirstChildElement("stress");
      stress && stress->FirstChild())
  {
    const char* sigma = stress->FirstChild()->Value();

    symmetric = true;
    IFEM::cout <<"\tStress: "<< sigma << std::endl;
    if (!stressSol)
    {
      stressSol = new STensorF(sigma,variables);
      parseDerivs(static_cast<STensorF*>(stressSol),stress);
    }
  }
  else if (vecSol)
    utl::getAttribute(elem,"symmetric",symmetric);
}


void AnaSol::parseFieldFunctions (const tinyxml2::XMLElement* elem,
                                  bool scalarSol)
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
  utl::getAttribute(elem,"level",level);

  // Lambda function for parsing the primary solution.
  auto&& parsePrimary = [file,basis,level,this]
    (const tinyxml2::XMLElement* elem, bool scalar,
     const char* keyw = "primary")
  {
    for (const tinyxml2::XMLElement* prim = elem->FirstChildElement(keyw);
         prim && prim->FirstChild(); prim = prim->NextSiblingElement(keyw))
    {
      const char* primary = prim->FirstChild()->Value();
      IFEM::cout <<"\tPrimary: "<< primary << std::endl;
      if (scalar)
        scalSol.push_back(new FieldFunction(file,basis,primary,level));
      else if (!vecSol)
        vecSol = new VecFieldFunction(file,basis,primary,level);
    }
  };

  // Lambda function for parsing the secondary solution.
  auto&& parseSecondary = [file,basis,level,this]
    (const tinyxml2::XMLElement* elem, bool scalar,
     const char* keyw = "secondary")
  {
    for (const tinyxml2::XMLElement* secs = elem->FirstChildElement(keyw);
         secs && secs->FirstChild(); secs = secs->NextSiblingElement(keyw))
    {
      const char* secondary = secs->FirstChild()->Value();
      IFEM::cout <<"\tSecondary: "<< secondary << std::endl;
      if (scalar)
        scalSecSol.push_back(new VecFieldFunction(file,basis,secondary,level));
      else if (!vecSecSol)
        vecSecSol = new TensorFieldFunction(file,basis,secondary,level);
    }
  };

  parsePrimary(elem,scalarSol,"primary");
  parsePrimary(elem,true,"scalarprimary");

  parseSecondary(elem,scalarSol,"secondary");
  parseSecondary(elem,true,"scalarsecondary");

  if (const tinyxml2::XMLElement* stress = elem->FirstChildElement("stress");
      stress && stress->FirstChild())
  {
    const char* sigma = stress->FirstChild()->Value();
    IFEM::cout <<"\tStress: "<< sigma << std::endl;
    if (!stressSol)
      stressSol = new STensorFieldFunction(file,basis,sigma,level);
  }

#else
  std::cerr <<" *** AnaSol::parseFieldFunctions: Compiled without HDF5 support"
            <<", no fields read."<< std::endl;
#endif
}


/*!
  If no secondary solution has been explicitly defined,
  this method uses derivation (automatic or finite difference)
  to obtain the gradient of the primary solution, which then
  is assigned as the analytic secondary solution field.
*/

void AnaSol::setupSecondarySolutions ()
{
  if (!stressSol) {
    if (scalSecSol.size() < scalSol.size())
      scalSecSol.resize(scalSol.size(),nullptr);
    for (size_t i = 0; i < scalSol.size(); ++i)
      if (!scalSecSol[i] && scalSol[i])
        scalSecSol[i] = new GradientFunc<RealFunc,VecFunc,Vec3>(*scalSol[i]);
  }
  if (vecSol) {
    if (symmetric && !stressSol)
      stressSol = new GradientFunc<VecFunc,STensorFunc,SymmTensor>(*vecSol);
    else if (!symmetric && !vecSecSol)
      vecSecSol = new GradientFunc<VecFunc,TensorFunc,Tensor>(*vecSol);
  }
}
