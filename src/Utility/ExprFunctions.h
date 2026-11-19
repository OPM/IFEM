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
#include <algorithm>
#include <memory>
#include <string>
#include <vector>

namespace ExprEval {
  template<class ArgType> class Expression;
  template<class ArgType> class FunctionList;
  template<class ArgType> class ValueList;
  extern int numError;
}


/*!
  \brief Common holder class for expression functions.
*/

template<class Scalar> class ExpressionHolder
{
protected:
  //! Type alias for expression tree
  using Expression = ExprEval::Expression<Scalar>;
  //! Type alias for function list
  using FunctionList = ExprEval::FunctionList<Scalar>;
  //! Type alias for value list
  using ValueList = ExprEval::ValueList<Scalar>;

  //! Roots of the expression tree
  std::vector< std::unique_ptr<Expression> > expr;
  //! Lists of functions
  std::vector< std::unique_ptr<FunctionList> >  f;
  //! Lists of variables and constants
  std::vector< std::unique_ptr<ValueList> >     v;

  //! \brief The constructor parses the expression string.
  explicit ExpressionHolder(const char* function);
  //! \brief Default destructor.
  virtual ~ExpressionHolder();

  //! \brief Sets an additional parameter in the variables and/or constants.
  void setParameter(const std::string& name, Real value);

  //! \brief Evaluates the function expression.
  Real evaluateExpression(size_t i) const;
};


/*!
  \brief A scalar-valued function, general expression.
*/

template<class Scalar>
class EvalFuncScalar : public ScalarFunc, private ExpressionHolder<Scalar>
{
  using FuncType = EvalFuncScalar<Scalar>; //!< Type alias for the function

  std::vector<Scalar*> arg; //!< Function argument values

  std::unique_ptr<FuncType> gradient; //!< First derivative expression

  Real dx; //!< Domain increment for calculation of numerical derivative

public:
  //! \brief The constructor parses the expression string.
  explicit EvalFuncScalar(const char* function, const char* x = "x",
                          Real eps = Real(1.0e-8));
  //! \brief Defaulted destructor.
  //! \details The implementation needs to be in compile unit so we have the
  //!          definition for the types of the unique_ptr's.
  virtual ~EvalFuncScalar();

  //! \brief Adds an expression function for a first derivative.
  void addDerivative(const std::string& function, const char* x = "x");

  //! \brief Sets an additional parameter in the function.
  void setParam(const std::string& name, Real value) override
  {
    this->setParameter(name,value);
  }

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
class EvalFuncSpatial : public RealFunc, private ExpressionHolder<Scalar>
{
  using FuncType = EvalFuncSpatial<Scalar>; //!< Type alias for the function

  //! \brief A struct representing a spatial function argument.
  struct Arg
  {
    Scalar* x; //!< X-coordinate
    Scalar* y; //!< Y-coordinate
    Scalar* z; //!< Z-coordinate
    Scalar* t; //!< Time

    //! \brief Assignment operator;
    const Arg& operator=(const Vec3& X) const
    {
      const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
      if (x) *x = X.x;
      if (y) *y = X.y;
      if (z) *z = X.z;
      if (t) *t = Xt ? Xt->t : Real(0);
      return *this;
    }

    //! \brief Indexing operator (one-based).
    const Scalar& operator()(int dir) const
    {
      static const Scalar dummy{};
      switch (dir) {
      case 1: return x ? *x : dummy;
      case 2: return y ? *y : dummy;
      case 3: return z ? *z : dummy;
      case 4: return t ? *t : dummy;
      }
      return dummy; // Index out of range, error?
    }

    //! \brief Checks whether a component is defined or not.
    bool validComp(int dir) const
    {
      switch (dir) {
      case 1: return x != nullptr;
      case 2: return y != nullptr;
      case 3: return z != nullptr;
      case 4: return t != nullptr;
      }
      return false;
    }
  };

  std::vector<Arg> arg; //!< Function argument values

  //! First and second order derivative expressions
  std::array<std::unique_ptr<FuncType>,10> derivative;

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

  //! \brief Sets an additional parameter in the function.
  void setParam(const std::string& name, Real value) override
  {
    this->setParameter(name,value);
  }

  //! \brief Returns whether the function is time-independent or not.
  bool isConstant() const override { return arg.empty() || !arg.front().t; }

  //! \brief Returns first-derivative of the function.
  Real deriv(const Vec3& X, int dir) const override;
  //! \brief Returns second-derivative of the function.
  Real dderiv(const Vec3& X, int i, int j) const override;

  //! \brief Evaluates first derivatives of the function.
  Vec3 gradient(const Vec3& X) const override;
  //! \brief Evaluates second derivatives of the function.
  SymmTensor hessian(const Vec3& X) const override;

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
  void addDerivative(const std::string& functions, const std::string& variables,
                     int d1, int d2 = 0);

  //! \brief Returns number of spatial dimension.
  size_t getNoSpaceDim() const { return nsd; }

protected:
  std::vector<std::unique_ptr<FuncType>> p; //!< Array of component expressions
  size_t nsd = 0; //!< Number of spatial dimensions
};


/*!
  \brief A general spatial expression function of any return type.
  \details The function is implemented as an array of EvalFunction objects.
*/

template<class ParentFunc, class Ret, class Scalar>
class EvalMultiFunction : public ParentFunc, public EvalFunctions<Scalar>
{
  //! Type alias for the function
  using FuncType = typename EvalFunctions<Scalar>::FuncType;

public:
  //! \brief The constructor parses the expression string for each component.
  explicit EvalMultiFunction(const std::string& functions,
                             const std::string& variables = "",
                             const Real epsX = 1e-8,
                             const Real epsT = 1e-12);

  //! \brief Returns whether the function is time-independent or not.
  bool isConstant() const override
  {
    return std::all_of(this->p.begin(), this->p.end(),
                       [](const std::unique_ptr<FuncType>& func)
                       { return func->isConstant(); });
  }

  //! \brief Returns the function type flag.
  unsigned char getType() const override { return 2; }

  //! \brief Returns first-derivative of the function.
  Ret deriv(const Vec3& X, int dir) const override;
  //! \brief Returns second-derivative of the function.
  Ret dderiv(const Vec3& X, int i, int j) const override;

  //! \brief Sets an additional parameter in the function.
  void setParam(const std::string& name, Real value) override
  {
    for (std::unique_ptr<FuncType>& func : this->p)
      func->setParam(name,value);
  }

protected:
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

#endif
