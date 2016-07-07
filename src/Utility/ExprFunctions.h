// $Id$
//==============================================================================
//!
//! \file ExprFunctions.h
//!
//! \date Dec 1 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Expression function implementations.
//!
//==============================================================================

#ifndef _EXPR_FUNCTIONS_H
#define _EXPR_FUNCTIONS_H

#include "Function.h"
#include <string>
#include <vector>
#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace ExprEval {
  class Expression;
  class FunctionList;
  class ValueList;
}


/*!
  \brief A scalar-valued function, general expression.
*/

class EvalFunc : public ScalarFunc
{
  std::vector<ExprEval::Expression*> expr; //!< Roots of the expression tree
  std::vector<ExprEval::FunctionList*>  f; //!< Lists of functions
  std::vector<ExprEval::ValueList*>     v; //!< Lists of variables and constants

  std::vector<Real*> arg; //!< Function argument values

public:
  //! \brief The constructor parses the expression string.
  EvalFunc(const char* function, const char* x = "x" );
  //! \brief The destructor frees the dynamically allocated objects.
  virtual ~EvalFunc();

  static int numError; //!< Error counter - set by the exception handler

protected:
  //! \brief Non-implemented copy constructor to disallow copying.
  EvalFunc(const EvalFunc&);
  //! \brief Non-implemented assigment operator to disallow copying.
  EvalFunc& operator=(const EvalFunc&);
  //! \brief Evaluates the function expression.
  virtual Real evaluate(const Real& x) const;
};


/*!
  \brief A scalar-valued spatial function, general function expression.
*/

class EvalFunction : public RealFunc
{
  std::vector<ExprEval::Expression*> expr; //!< Roots of the expression tree
  std::vector<ExprEval::FunctionList*>  f; //!< Lists of functions
  std::vector<ExprEval::ValueList*>     v; //!< Lists of variables and constants

  //! \brief A struct representing a spatial function argument.
  struct Arg
  {
    Real* x; //!< X-coordinate
    Real* y; //!< Y-coordinate
    Real* z; //!< Z-coordinate
    Real* t; //!< Time
  };

  std::vector<Arg> arg; //!< Function argument values

  bool IAmConstant; //!< Indicates whether the time coordinate is given or not

public:
  //! \brief The constructor parses the expression string.
  EvalFunction(const char* function);
  //! \brief The destructor frees the dynamically allocated objects.
  virtual ~EvalFunction();

  //! \brief Returns whether the function is time-independent or not.
  virtual bool isConstant() const { return IAmConstant; }

protected:
  //! \brief Evaluates the function expression.
  virtual Real evaluate(const Vec3& X) const;
  //! \brief Non-implemented copy constructor to disallow copying.
  EvalFunction(const EvalFunction&);
  //! \brief Non-implemented assignment operator to disallow copying.
  EvalFunction& operator=(const EvalFunction&);
};


/*!
  \brief A general spatial expression function of any return type.
  \details The function is implemented as an array of EvalFunction objects.
*/

template <class ParentFunc, class Ret>
class EvalMultiFunction : public ParentFunc
{
  std::vector<EvalFunction*> p; //!< Array of component expressions

public:
  //! \brief The constructor parses the expression string for each component.
  EvalMultiFunction<ParentFunc,Ret>(const std::string& functions,
                                    const std::string& variables = "")
  {
    size_t pos = functions.find("|"), pos2 = 0;
    for (int i = 0; pos2 < functions.size(); i++)
    {
      std::string func(variables);
      if (!func.empty() && func[func.size()-1] != ';')
        func += ';';
      if (pos == std::string::npos)
        func += functions.substr(pos2);
      else
        func += functions.substr(pos2,pos-pos2);
      p.push_back(new EvalFunction(func.c_str()));
      pos2 = pos > 0 && pos < std::string::npos ? pos+1 : pos;
      pos = functions.find("|",pos+1);
    }
  }

  //! \brief The destructor frees the dynamically allocated function components.
  virtual ~EvalMultiFunction<ParentFunc,Ret>()
  {
    for (size_t i = 0; i < p.size(); i++)
      delete p[i];
  }

  //! \brief Returns whether the function is time-independent or not.
  virtual bool isConstant() const
  {
    for (size_t i = 0; i < p.size(); i++)
      if (!p[i]->isConstant()) return false;
    return true;
  }

protected:
  //! \brief Evaluates the function expressions.
  virtual Ret evaluate(const Vec3& X) const;
};

//! Vector-valued function expression
typedef EvalMultiFunction<VecFunc,Vec3>           VecFuncExpr;
//! Tensor-valued function expression
typedef EvalMultiFunction<TensorFunc,Tensor>      TensorFuncExpr;
//! Symmetric tensor-valued function expression
typedef EvalMultiFunction<STensorFunc,SymmTensor> STensorFuncExpr;

//! \brief Specialization for vector functions.
template<> Vec3 VecFuncExpr::evaluate(const Vec3& X) const;

//! \brief Specialization for tensor functions.
template<> Tensor TensorFuncExpr::evaluate(const Vec3& X) const;

//! \brief Specialization for symmetric tensor functions.
template<> SymmTensor STensorFuncExpr::evaluate(const Vec3& X) const;

#endif
