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
  ExprEval::Expression* expr; //!< Pointer to the root of the expression tree
  ExprEval::FunctionList*  f; //!< Pointer to list of function in the expression
  ExprEval::ValueList*     v; //!< Pointer to list of variables and constants

  real* arg; //!< Pointer to the function argument

public:
  //! \brief The constructor parses the expression string.
  EvalFunc(const char* function, const char* x = "x" );
  //! \brief The destructor frees the dynamically allocated objects.
  virtual ~EvalFunc();

protected:
  //! \brief Evaluates the function expression.
  virtual real evaluate(const real& x) const;

private:
#ifdef USE_OPENMP
  omp_lock_t lock;
#endif
};


/*!
  \brief A scalar-valued spatial function, general function expression.
*/

class EvalFunction : public RealFunc
{
  ExprEval::Expression* expr; //!< Pointer to the root of the expression tree
  ExprEval::FunctionList*  f; //!< Pointer to list of function in the expression
  ExprEval::ValueList*     v; //!< Pointer to list of variables and constants

  real* x; //!< Pointer to the X-coordinate of the function argument
  real* y; //!< Pointer to the Y-coordinate of the function argument
  real* z; //!< Pointer to the Z-coordinate of the function argument
  real* t; //!< Pointer to the time coordinate of the function argument

public:
  //! \brief The constructor parses the expression string.
  EvalFunction(const char* function);
  //! \brief The destructor frees the dynamically allocated objects.
  virtual ~EvalFunction();

protected:
  //! \brief Evaluates the function expression.
  virtual real evaluate(const Vec3& X) const;

private:
#ifdef USE_OPENMP
  omp_lock_t lock;
#endif
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
				    const std::string& variables="")
  {
    size_t pos = functions.find("|"), pos2 = 0;
    for (int i = 0; pos2 < functions.size(); i++)
    {
      std::string func(variables);
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
