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
#include "TensorFunction.h"
#include <string>
#include <vector>
#include <array>

namespace ExprEval {
  template<class ArgType> class Expression;
  template<class ArgType> class FunctionList;
  template<class ArgType> class ValueList;
}


/*!
  \brief A scalar-valued function, general expression.
*/

class EvalFunc : public ScalarFunc
{
  using Expression = ExprEval::Expression<Real>;     //!< Type alias for expression tree
  using FunctionList = ExprEval::FunctionList<Real>; //!< Type alias for function list
  using ValueList = ExprEval::ValueList<Real>;       //!< Type alias for value list
  std::vector<Expression*> expr; //!< Roots of the expression tree
  std::vector<FunctionList*>  f; //!< Lists of functions
  std::vector<ValueList*>     v; //!< Lists of variables and constants

  std::vector<Real*> arg; //!< Function argument values

  EvalFunc* gradient; //!< First derivative expression

  Real dx; //!< Domain increment for calculation of numerical derivative

public:
  static int numError; //!< Error counter - set by the exception handler

  //! \brief The constructor parses the expression string.
  explicit EvalFunc(const char* function, const char* x = "x",
                    Real eps = Real(1.0e-8));
  //! \brief The destructor frees the dynamically allocated objects.
  virtual ~EvalFunc();

  //! \brief Adds an expression function for a first derivative.
  void derivative(const std::string& function, const char* x = "x");

  //! \brief Returns whether the function is time-independent or not.
  virtual bool isConstant() const { return false; }

  //! \brief Returns the first-derivative of the function.
  virtual Real deriv(Real x) const;

protected:
  //! \brief Non-implemented copy constructor to disallow copying.
  EvalFunc(const EvalFunc&) = delete;
  //! \brief Non-implemented assignment operator to disallow copying.
  EvalFunc& operator=(const EvalFunc&) = delete;
  //! \brief Evaluates the function expression.
  virtual Real evaluate(const Real& x) const;

  //! \brief Cleans up the allocated data.
  void cleanup();
};


/*!
  \brief A scalar-valued spatial function, general function expression.
*/

class EvalFunction : public RealFunc
{
  using Expression = ExprEval::Expression<Real>;     //!< Type alias for expression tree
  using FunctionList = ExprEval::FunctionList<Real>; //!< Type alias for function list
  using ValueList = ExprEval::ValueList<Real>;       //!< Type alias for value list
  std::vector<Expression*> expr; //!< Roots of the expression tree
  std::vector<FunctionList*>  f; //!< Lists of functions
  std::vector<ValueList*>     v; //!< Lists of variables and constants

  //! \brief A struct representing a spatial function argument.
  struct Arg
  {
    Real* x; //!< X-coordinate
    Real* y; //!< Y-coordinate
    Real* z; //!< Z-coordinate
    Real* t; //!< Time
  };

  std::vector<Arg> arg; //!< Function argument values

  std::array<EvalFunction*,3> gradient;  //!< First derivative expressions
  std::array<EvalFunction*,6> dgradient; //!< Second derivative expressions

  bool IAmConstant; //!< Indicates whether the time coordinate is given or not

  Real dx; //!< Domain increment for calculation of numerical derivative
  Real dt; //!< Domain increment for calculation of numerical time-derivative

public:
  //! \brief The constructor parses the expression string.
  explicit EvalFunction(const char* function,
                        Real epsX = Real(1.0e-8), Real epsT = Real(1.0e-12));
  //! \brief The destructor frees the dynamically allocated objects.
  virtual ~EvalFunction();

  //! \brief Adds an expression function for a first or second derivative.
  void addDerivative(const std::string& function, const std::string& variables,
                     int d1, int d2 = 0);

  //! \brief Returns whether the function is time-independent or not.
  virtual bool isConstant() const { return IAmConstant; }

  //! \brief Returns first-derivative of the function.
  virtual Real deriv(const Vec3& X, int dir) const;
  //! \brief Returns second-derivative of the function.
  virtual Real dderiv(const Vec3& X, int dir1, int dir2) const;

  //! \brief Set an additional parameter in the function.
  void setParam(const std::string& name, double value);

protected:
  //! \brief Non-implemented copy constructor to disallow copying.
  EvalFunction(const EvalFunction&) = delete;
  //! \brief Non-implemented assignment operator to disallow copying.
  EvalFunction& operator=(const EvalFunction&) = delete;
  //! \brief Evaluates the function expression.
  virtual Real evaluate(const Vec3& X) const;

  //! \brief Cleans up the allocated data.
  void cleanup();
};


/*!
  \brief A base class for multi-component expression functions.
*/

class EvalFunctions
{
protected:
  //! \brief The constructor parses the expression string for each component.
  EvalFunctions(const std::string& functions, const std::string& variables);
  //! \brief The destructor frees the dynamically allocated function components.
  virtual ~EvalFunctions();

public:
  //! \brief Adds an expression function for a first or second derivative.
  void addDerivative(const std::string& functions,
                     const std::string& variables, int d1, int d2 = 0);

protected:
  std::vector<EvalFunction*> p; //!< Array of component expressions
};


/*!
  \brief A general spatial expression function of any return type.
  \details The function is implemented as an array of EvalFunction objects.
*/

template <class ParentFunc, class Ret>
class EvalMultiFunction : public ParentFunc, public EvalFunctions
{
  size_t nsd; //!< Number of spatial dimensions

public:
  //! \brief The constructor parses the expression string for each component.
  EvalMultiFunction(const std::string& functions,
                    const std::string& variables = "")
    : EvalFunctions(functions,variables), nsd(0) { this->setNoDims(); }

  //! \brief Empty destructor.
  virtual ~EvalMultiFunction() {}

  //! \brief Returns whether the function is time-independent or not.
  virtual bool isConstant() const
  {
    for (EvalFunction* func : p)
      if (!func->isConstant())
        return false;
    return true;
  }

  //! \brief Returns the function type flag.
  virtual unsigned char getType() const { return 2; }

  //! \brief Returns first-derivative of the function.
  virtual Ret deriv(const Vec3& X, int dir) const;
  //! \brief Returns second-derivative of the function.
  virtual Ret dderiv(const Vec3& X, int dir1, int dir2) const;

  //! \brief Set an additional parameter in the function.
  void setParam(const std::string& name, double value)
  {
    for (EvalFunction* func : p)
      func->setParam(name, value);
  }

protected:
  //! \brief Sets the number of spatial dimensions (default implementation).
  void setNoDims() { ParentFunc::ncmp = nsd = p.size(); }

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
template<> void TensorFuncExpr::setNoDims();
//! \brief Specialization for tensor functions.
template<> Tensor TensorFuncExpr::evaluate(const Vec3& X) const;
//! \brief Specialization for tensor functions.
template<> Tensor TensorFuncExpr::deriv(const Vec3& X, int dir) const;
//! \brief Specialization for tensor functions.
template<> Tensor TensorFuncExpr::dderiv(const Vec3& X, int d1, int d2) const;

//! \brief Specialization for symmetric tensor functions.
template<> void STensorFuncExpr::setNoDims();
//! \brief Specialization for symmetric tensor functions.
template<> SymmTensor STensorFuncExpr::evaluate(const Vec3& X) const;
//! \brief Specialization for symmetric tensor functions.
template<> SymmTensor STensorFuncExpr::deriv(const Vec3& X, int dir) const;
//! \brief Specialization for symmetric tensor functions.
template<> SymmTensor STensorFuncExpr::dderiv(const Vec3& X,
                                              int d1, int d2) const;

#endif
