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

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace ExprEval {
  template<class ArgType> class Expression;
  template<class ArgType> class FunctionList;
  template<class ArgType> class ValueList;
}


/*!
  \brief A scalar-valued function, general expression.
*/

template<class Scalar>
class EvalFuncScalar : public ScalarFunc
{
  using Expression = ExprEval::Expression<Scalar>;     //!< Type alias for expression tree
  using FunctionList = ExprEval::FunctionList<Scalar>; //!< Type alias for function list
  using ValueList = ExprEval::ValueList<Scalar>;       //!< Type alias for value list
  using FuncType = EvalFuncScalar<Scalar>;             //!< Type alias for function
  std::vector<std::unique_ptr<Expression>> expr; //!< Roots of the expression tree
  std::vector<std::unique_ptr<FunctionList>>  f; //!< Lists of functions
  std::vector<std::unique_ptr<ValueList>>     v; //!< Lists of variables and constants

  std::vector<Scalar*> arg; //!< Function argument values

  std::unique_ptr<FuncType> gradient; //!< First derivative expression

  Real dx; //!< Domain increment for calculation of numerical derivative

public:
  static int numError; //!< Error counter - set by the exception handler

  //! \brief The constructor parses the expression string.
  explicit EvalFuncScalar(const char* function, const char* x = "x",
                          Real eps = Real(1.0e-8));
  //! \brief Defaulted destructor.
  //! \details The implementation needs to be in compile unit so we have the
  //!          definition for the types of the unique_ptr's.
  virtual ~EvalFuncScalar();

  //! \brief Adds an expression function for a first derivative.
  void addDerivative(const std::string& function, const char* x = "x");

  //! \brief Returns whether the function is time-independent or not.
  bool isConstant() const override { return false; }

  //! \brief Returns the first-derivative of the function.
  Real deriv(Real x) const override;

protected:
  //! \brief Evaluates the function expression.
  Real evaluate(const Real& x) const override;
};


/*!
  \brief A scalar-valued spatial function, general function expression.
*/

template<class Scalar>
class EvalFuncSpatial : public RealFunc
{
  using Expression = ExprEval::Expression<Scalar>;     //!< Type alias for expression tree
  using FunctionList = ExprEval::FunctionList<Scalar>; //!< Type alias for function list
  using ValueList = ExprEval::ValueList<Scalar>;       //!< Type alias for value list
  using FuncType = EvalFuncSpatial<Scalar>;        //!< Type alias for function
  std::vector<std::unique_ptr<Expression>> expr; //!< Roots of the expression tree
  std::vector<std::unique_ptr<FunctionList>>  f; //!< Lists of functions
  std::vector<std::unique_ptr<ValueList>>     v; //!< Lists of variables and constants

  //! \brief A struct representing a spatial function argument.
  struct Arg
  {
    Scalar* x; //!< X-coordinate
    Scalar* y; //!< Y-coordinate
    Scalar* z; //!< Z-coordinate
    Scalar* t; //!< Time

    //! \brief Returns a const ref to a member.
    //! \param dir One-based index to member
    const Scalar& get(int dir) const
    {
      switch (dir) {
        case  1: return *x;
        case  2: return *y;
        case  3: return *z;
        case  4: return *t;
        default: return *x;
      }
    }
  };

  std::vector<Arg> arg; //!< Function argument values

  std::array<std::unique_ptr<FuncType>,4> derivative1; //!< First order derivative expressions
  std::array<std::unique_ptr<FuncType>,6> derivative2; //!< Second order derivative expressions

  bool IAmConstant; //!< Indicates whether the time coordinate is given or not

  Real dx; //!< Domain increment for calculation of numerical derivative
  Real dt; //!< Domain increment for calculation of numerical time-derivative

public:
  //! \brief The constructor parses the expression string.
  explicit EvalFuncSpatial(const char* function,
                               Real epsX = Real(1.0e-8), Real epsT = Real(1.0e-12));
  //! \brief Defaulted destructor.
  //! \details The implementation needs to be in compile unit so we have the
  //!          definition for the types of the unique_ptr's.
  virtual ~EvalFuncSpatial();

  //! \brief Adds an expression function for a first or second derivative.
  void addDerivative(const std::string& function, const std::string& variables,
                     int d1, int d2 = 0);

  //! \brief Returns whether the function is time-independent or not.
  bool isConstant() const override { return IAmConstant; }

  //! \brief Returns first-derivative of the function.
  Real deriv(const Vec3& X, int dir) const override;
  //! \brief Returns second-derivative of the function.
  Real dderiv(const Vec3& X, int dir1, int dir2) const override;

  //! \brief Set an additional parameter in the function.
  void setParam(const std::string& name, double value);

  //! \brief Evaluates first derivatives of the function.
  Vec3 gradient(const Vec3& X) const override
  {
    return this->RealFunc::gradient(X);
  }

  //! \brief Evaluates first derivatives of the function.
  SymmTensor hessian(const Vec3& X) const override
  {
    return this->RealFunc::hessian(X);
  }

protected:
  //! \brief Evaluates the function expression.
  Real evaluate(const Vec3& X) const override;
};


/*!
  \brief A base class for multi-component expression functions.
*/

template<class Scalar>
class EvalFunctions
{
protected:
  using FuncType = EvalFuncSpatial<Scalar>; //!< Type alias for function

  //! \brief The constructor parses the expression string for each component.
  EvalFunctions(const std::string& functions, const std::string& variables,
                const Real epsX, const Real epsT);
  //! \brief Defaulted destructor.
  //! \details The implementation needs to be in compile unit so we have the
  //!          definition for the types of the unique_ptr's.
  virtual ~EvalFunctions();

public:
  //! \brief Adds an expression function for a first or second derivative.
  void addDerivative(const std::string& functions,
                     const std::string& variables, int d1, int d2 = 0);

protected:
  std::vector<std::unique_ptr<FuncType>> p; //!< Array of component expressions
};


/*!
  \brief A general spatial expression function of any return type.
  \details The function is implemented as an array of EvalFunction objects.
*/

template <class ParentFunc, class Ret, class Scalar>
class EvalMultiFunction : public ParentFunc, public EvalFunctions<Scalar>
{
  size_t nsd; //!< Number of spatial dimensions
  using FuncType = typename EvalFunctions<Scalar>::FuncType; //!< Type alias for function

public:
  //! \brief The constructor parses the expression string for each component.
  explicit EvalMultiFunction(const std::string& functions,
                             const std::string& variables = "",
                             const Real epsX = 1e-8,
                             const Real epsT = 1e-12)
    : EvalFunctions<Scalar>(functions,variables,epsX,epsT), nsd(0) { this->setNoDims(); }

  //! \brief Empty destructor.
  virtual ~EvalMultiFunction() {}

  //! \brief Returns whether the function is time-independent or not.
  bool isConstant() const override
  {
    for (const std::unique_ptr<FuncType>& func : this->p)
      if (!func->isConstant())
        return false;
    return true;
  }

  //! \brief Returns the function type flag.
  unsigned char getType() const override { return 2; }

  //! \brief Returns first-derivative of the function.
  Ret deriv(const Vec3& X, int dir) const override;
  //! \brief Returns second-derivative of the function.
  Ret dderiv(const Vec3& X, int dir1, int dir2) const override;

  //! \brief Set an additional parameter in the function.
  void setParam(const std::string& name, double value)
  {
    for (std::unique_ptr<FuncType>& func : this->p)
      func->setParam(name, value);
  }

  //! \brief Returns number of spatial dimension.
  size_t getNoSpaceDim() const { return nsd; }

protected:
  //! \brief Sets the number of spatial dimensions (default implementation).
  void setNoDims();

  //! \brief Evaluates the function expressions.
  Ret evaluate(const Vec3& X) const override;

  //! \brief Returns the gradient of the function as a 1D array.
  std::vector<Real> evalGradient(const Vec3& X) const override;

  //! \brief Returns the second derivatives of the function as a 1D array.
  std::vector<Real> evalHessian(const Vec3& X) const override;

  //! \brief Returns the time derivatives of the function as a 1D array.
  std::vector<Real> evalTimeDerivative(const Vec3& X) const override;
};

//! Scalar-valued function expression
using EvalFunc = EvalFuncScalar<Real>;
//! Scalar-valued spatial function expression
using EvalFunction = EvalFuncSpatial<Real>;
//! Vector-valued function expression
using VecFuncExpr = EvalMultiFunction<VecFunc,Vec3,Real>;
//! Tensor-valued function expression
using TensorFuncExpr = EvalMultiFunction<TensorFunc,Tensor,Real>;
//! Symmetric tensor-valued function expression
using STensorFuncExpr = EvalMultiFunction<STensorFunc,SymmTensor,Real>;

//! \brief Explicit instantiation of error flag.
template<> int EvalFunc::numError;

#endif
