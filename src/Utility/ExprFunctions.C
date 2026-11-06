// $Id$
//==============================================================================
//!
//! \file ExprFunctions.C
//!
//! \date Dec 1 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Expression function implementations.
//!
//==============================================================================

#include "ExprFunctions.h"
#include "Functions.h"
#include "expreval.h"
#include <autodiff/reverse/var.hpp>
#ifdef USE_OPENMP
#include <omp.h>
#endif
#include <type_traits>


namespace ExprEval {
  //! The error counter is global, and is incremented by run-time evaluation errors.
  int numError = 0;
}


namespace
{

/*!
  Prints an error message with the exception occured to std::cerr.
*/

void ExprException (const ExprEval::Exception& exc, const char* task,
                    const char* function = nullptr)
{
  std::cerr <<"\n *** Error "<< task <<" function";
  if (function)
    std::cerr <<" \""<< function <<"\"";
  if (!exc.GetValue().empty())
    std::cerr <<", "<< exc.GetValue();

  switch (exc.GetType()) {
  case ExprEval::Exception::Type_NotFoundException:
    std::cerr <<": Not found";
    break;
  case ExprEval::Exception::Type_AlreadyExistsException:
    std::cerr <<": Already exists";
    break;
  case ExprEval::Exception::Type_NullPointerException:
    std::cerr <<": Null pointer";
    break;
  case ExprEval::Exception::Type_MathException:
    std::cerr <<": Math exception, "<< exc.GetError();
    break;
  case ExprEval::Exception::Type_DivideByZeroException:
    std::cerr <<": Division by zero";
    break;
  case ExprEval::Exception::Type_NoValueListException:
    std::cerr <<": No value list";
    break;
  case ExprEval::Exception::Type_NoFunctionListException:
    std::cerr <<": No function list";
    break;
  case ExprEval::Exception::Type_AbortException:
    std::cerr <<": Abort";
    break;
  case ExprEval::Exception::Type_EmptyExpressionException:
    std::cerr <<": Empty expression";
    break;
  case ExprEval::Exception::Type_UnknownTokenException:
    std::cerr <<": Unknown token";
    break;
  case ExprEval::Exception::Type_InvalidArgumentCountException:
    std::cerr <<": Invalid argument count";
    break;
  case ExprEval::Exception::Type_ConstantAssignException:
    std::cerr <<": Constant assign";
    break;
  case ExprEval::Exception::Type_ConstantReferenceException:
    std::cerr <<": Constant reference";
    break;
  case ExprEval::Exception::Type_SyntaxException:
    std::cerr <<": Syntax error";
    break;
  case ExprEval::Exception::Type_UnmatchedParenthesisException:
    std::cerr <<": Unmatched parenthesis";
    break;
  default:
    std::cerr <<": Unknown exception";
  }
  std::cerr << std::endl;
  ExprEval::numError++;
}


/*!
  \brief Static helper that splits a function expression into components.
*/

std::vector<std::string> splitComps (const std::string& functions,
                                     const std::string& variables)
{
  std::vector<std::string> comps;
  size_t pos1 = functions.find("|");
  size_t pos2 = 0;
  while (pos2 < functions.size())
  {
    std::string func(variables);
    if (!func.empty() && func[func.size()-1] != ';')
      func += ';';
    if (pos1 == std::string::npos)
      func += functions.substr(pos2);
    else
      func += functions.substr(pos2,pos1-pos2);
    comps.push_back(func);
    pos2 = pos1 > 0 && pos1 < std::string::npos ? pos1+1 : pos1;
    pos1 = functions.find("|",pos1+1);
  }

  return comps;
}


/*!
  \brief Helper template to get size and dimension of a return type.
*/

template<class ArgType>
std::pair<size_t,size_t> getNoDims (size_t psize);


/*!
  \brief Template specialization for Vec3.
*/

template<>
std::pair<size_t,size_t> getNoDims<Vec3> (size_t psize)
{
  return {psize, psize};
}


/*!
  \brief Template specialization for Tensor.
*/

template<>
std::pair<size_t,size_t> getNoDims<Tensor> (size_t psize)
{
  size_t nsd = 0;
  if (psize > 8)
    nsd = 3;
  else if (psize > 3)
    nsd = 2;
  else if (psize > 0)
    nsd = 1;

  return {nsd, nsd*nsd};
}


/*!
  \brief Template specialization for SymmTensor.
*/

template<>
std::pair<size_t,size_t> getNoDims<SymmTensor> (size_t psize)
{
  size_t nsd = 0;
  if (psize > 5)
    nsd = 3;
  else if (psize > 2)
    nsd = 2;
  else if (psize > 0)
    nsd = 1;

  return {nsd, psize == 4 ? 4 : (nsd+1)*nsd/2};
}


/*!
  \brief Helper to obtain Voigt index.
*/

int voigtIdx (int d1, int d2)
{
  if (d1 > d2)
    std::swap(d1,d2); // Assuming symmetry

  if (d1 < 1 || d2 > 3)
    return -1; // Out-of-range
  else if (d2-d1 == 0) // diagonal term, 11, 22 and 33
    return d1-1;
  else if (d2-d1 == 1) // off-diagonal term, 12 and 23
    return d2+1;
  else // off-diagonal term, 13
    return 5;
}

}


template<class Scalar>
ExpressionHolder<Scalar>::ExpressionHolder (const char* function)
{
  try {
#ifdef USE_OPENMP
    const size_t nalloc = omp_get_max_threads();
#else
    const size_t nalloc = 1;
#endif
    expr.resize(nalloc);
    f.resize(nalloc);
    v.resize(nalloc);
    for (size_t i = 0; i < nalloc; i++)
    {
      expr[i] = std::make_unique<Expression>();
      f[i] = std::make_unique<FunctionList>();
      v[i] = std::make_unique<ValueList>();
      f[i]->AddDefaultFunctions();
      v[i]->AddDefaultValues();
      expr[i]->SetFunctionList(f[i].get());
      expr[i]->SetValueList(v[i].get());
      expr[i]->Parse(function);
    }
  }
  catch (ExprEval::Exception& e) {
    ExprException(e,"parsing",function);
  }
}


template<class Scalar>
ExpressionHolder<Scalar>::~ExpressionHolder () = default;


template<class Scalar>
void ExpressionHolder<Scalar>::setParameter (const std::string& name, Real val)
{
  auto setVal = [&name,val](ValueList& v)
  {
    Scalar* address = v.GetAddress(name);
    if (!address)
      v.Add(name,val,false);
    else
      *address = val;
  };

#ifdef USE_OPENMP
  if (omp_in_parallel())
    setVal(*v[omp_get_thread_num()]);
  else
#endif
    for (std::unique_ptr<ValueList>& v1 : v)
      setVal(*v1);
}


template<class Scalar>
Real ExpressionHolder<Scalar>::evaluateExpression (size_t i) const
{
  Real result = Real(0);
  if (i >= expr.size() || !expr[i].get())
    return result;

  try {
    if constexpr (std::is_same_v<Scalar,Real>)
      result = expr[i]->Evaluate();
    else
      result = expr[i]->Evaluate().expr->val;
  }
  catch (ExprEval::Exception& e) {
    ExprException(e,"evaluating expression");
  }

  return result;
}


template<class Scalar>
EvalFuncScalar<Scalar>::EvalFuncScalar (const char* function,
                                        const char* x, Real eps)
  : ExpressionHolder<Scalar>(function), dx(eps)
{
  if (ExprEval::numError > 0)
    return; // Faulty expression

  arg.resize(this->expr.size());
  for (size_t i = 0; i < arg.size(); i++)
    arg[i] = this->v[i]->GetAddress(x);
}


template<class Scalar>
EvalFuncScalar<Scalar>::~EvalFuncScalar () = default;


template<class Scalar>
void EvalFuncScalar<Scalar>::addDerivative (const std::string& function,
                                            const char* x)
{
  if (!gradient)
    gradient = std::make_unique<FuncType>(function.c_str(),x);
}


template<class Scalar>
Real EvalFuncScalar<Scalar>::evaluate (const Real& x) const
{
#ifdef USE_OPENMP
  const size_t i = omp_get_thread_num();
#else
  const size_t i = 0;
#endif
  if (i < arg.size() && arg[i])
    *arg[i] = x;

  return this->evaluateExpression(i);
}


template<>
Real EvalFuncScalar<Real>::deriv (Real x) const
{
  if (gradient)
    return gradient->evaluate(x);

  // Evaluate derivative using central difference
  return (this->evaluate(x+0.5*dx) - this->evaluate(x-0.5*dx)) / dx;
}


template<>
Real EvalFuncScalar<autodiff::var>::deriv (Real x) const
{
  if (gradient)
    return gradient->evaluate(x);

#ifdef USE_OPENMP
  const size_t i = omp_get_thread_num();
#else
  const size_t i = 0;
#endif
  if (i >= arg.size() || !arg[i])
    return Real(0);

  try {
    *arg[i] = x;
    return derivativesx(this->expr[i]->Evaluate(),
                        autodiff::wrt(*arg[i]))[0].expr->val;
  }
  catch (ExprEval::Exception& e) {
    ExprException(e,"evaluating expression");
  }

  return Real(0);
}


template<class Scalar>
EvalFuncSpatial<Scalar>::EvalFuncSpatial (const char* function,
                                          Real epsX, Real epsT)
  : ExpressionHolder<Scalar>(function), dx(epsX), dt(epsT)
{
  if (ExprEval::numError > 0)
    return; // Faulty expression

  // Check if the expression is time-dependent
  const bool isTimeDependent = utl::isTimeExpression(function);

  arg.resize(this->expr.size());
  for (size_t i = 0; i < arg.size(); i++)
  {
    arg[i].x = this->v[i]->GetAddress("x");
    arg[i].y = this->v[i]->GetAddress("y");
    arg[i].z = this->v[i]->GetAddress("z");
    if (isTimeDependent)
      arg[i].t = this->v[i]->GetAddress("t");
  }
}


template<class Scalar>
EvalFuncSpatial<Scalar>::~EvalFuncSpatial () = default;


template<class Scalar>
void EvalFuncSpatial<Scalar>::addDerivative (const std::string& function,
                                             const std::string& variables,
                                             int d1, int d2)
{
  if (d1 > 0 && d1 <= 4 && d2 < 1)
    --d1; // A first order derivative is specified
  else if ((d1 = voigtIdx(d1,d2)) >= 0)
    d1 += 4; // A second order derivative is specified
  else
    return;

  if (!derivative[d1])
    derivative[d1] = std::make_unique<FuncType>((variables+function).c_str());
}


template<class Scalar>
Real EvalFuncSpatial<Scalar>::evaluate (const Vec3& X) const
{
#ifdef USE_OPENMP
  const size_t i = omp_get_thread_num();
#else
  const size_t i = 0;
#endif
  if (i < arg.size())
    arg[i] = X;

  return this->evaluateExpression(i);
}


template<class Scalar>
Vec3 EvalFuncSpatial<Scalar>::gradient (const Vec3& X) const
{
  return this->RealFunc::gradient(X);
}


template<class Scalar>
SymmTensor EvalFuncSpatial<Scalar>::hessian (const Vec3& X) const
{
  return this->RealFunc::hessian(X);
}


template<>
Real EvalFuncSpatial<Real>::deriv (const Vec3& X, int dir) const
{
  if (dir < 1 || (dir > 3 && this->isConstant()))
    return Real(0);
  else if (dir < 4)
  {
    if (derivative[--dir])
      return derivative[dir]->evaluate(X);

    // Evaluate spatial derivative using central difference
    Vec4 X0, X1;
    X0.assign(X); X0[dir] -= 0.5*dx;
    X1.assign(X); X1[dir] += 0.5*dx;
    return (this->evaluate(X1) - this->evaluate(X0)) / dx;
  }
  else
  {
    if (derivative[3])
      return derivative[3]->evaluate(X);

    // Evaluate time-derivative using central difference
    Vec4 X0, X1;
    X0.assign(X); X0.t -= 0.5*dt;
    X1.assign(X); X1.t += 0.5*dt;
    return (this->evaluate(X1) - this->evaluate(X0)) / dt;
  }
}


template<>
Real EvalFuncSpatial<autodiff::var>::deriv (const Vec3& X, int dir) const
{
#ifdef USE_OPENMP
  const size_t i = omp_get_thread_num();
#else
  const size_t i = 0;
#endif
  if (i >= arg.size() || !arg[i].validComp(dir))
    return Real(0);

  arg[i] = X;

  // Evaluate spatial derivative using auto-diff
  return derivativesx(this->expr[i]->Evaluate(),
                      autodiff::wrt(arg[i](dir)))[0].expr->val;
}


template<>
Real EvalFuncSpatial<Real>::dderiv (const Vec3& X, int i, int j) const
{
  if ((i = voigtIdx(i,j)) < 0)
    return Real(0);
  else
    i += 4;

  return derivative[i] ? derivative[i]->evaluate(X) : Real(0);
}


template<>
Real EvalFuncSpatial<autodiff::var>::dderiv (const Vec3& X, int i, int j) const
{
#ifdef USE_OPENMP
  const size_t t = omp_get_thread_num();
#else
  const size_t t = 0;
#endif
  if (t >= arg.size() || !arg[t].validComp(i) || !arg[t].validComp(j))
    return Real(0);

  arg[t] = X;
  return derivativesx(derivativesx(this->expr[t]->Evaluate(),
                                   autodiff::wrt(arg[t](i)))[0],
                                   autodiff::wrt(arg[t](j)))[0].expr->val;
}


template<>
Vec3 EvalFuncSpatial<autodiff::var>::gradient (const Vec3& X) const
{
#ifdef USE_OPENMP
  const size_t i = omp_get_thread_num();
#else
  const size_t i = 0;
#endif
  if (i >= arg.size())
    return Vec3();

  arg[i] = X;

  const auto dx = derivativesx(this->expr[i]->Evaluate(),
                               autodiff::wrt(arg[i](1), arg[i](2), arg[i](3)));

  return Vec3(dx[0].expr->val, dx[1].expr->val, dx[2].expr->val);
}


template<>
SymmTensor EvalFuncSpatial<autodiff::var>::hessian (const Vec3& X) const
{
#ifdef USE_OPENMP
  const size_t i = omp_get_thread_num();
#else
  const size_t i = 0;
#endif
  if (i >= arg.size())
    return SymmTensor(3);

  arg[i] = X;

  const auto dx = derivativesx(this->expr[i]->Evaluate(),
                               autodiff::wrt(arg[i](1), arg[i](2), arg[i](3)));

  const auto [uxx, uxy, uxz] =
    derivativesx(dx[0], autodiff::wrt(arg[i](1), arg[i](2), arg[i](3)));

  const auto [uyy, uyz] =
    derivativesx(dx[1], autodiff::wrt(arg[i](2), arg[i](3)));

  const auto [uzz] =
    derivativesx(dx[2], autodiff::wrt(arg[i](3)));

  return SymmTensor({uxx.expr->val, uyy.expr->val, uzz.expr->val,
                     uxy.expr->val, uyz.expr->val, uxz.expr->val});
}


template<class Scalar>
EvalFunctions<Scalar>::EvalFunctions (const std::string& functions,
                                      const std::string& variables,
                                      const Real epsX, const Real epsT)
{
  std::vector<std::string> components = splitComps(functions,variables);
  for (const std::string& comp : components)
    p.emplace_back(std::make_unique<FuncType>(comp.c_str(),epsX,epsT));
}


template<class Scalar>
EvalFunctions<Scalar>::~EvalFunctions () = default;


template<class Scalar>
void EvalFunctions<Scalar>::addDerivative (const std::string& functions,
                                           const std::string& variables,
                                           int d1, int d2)
{
  std::vector<std::string> components = splitComps(functions,variables);
  for (size_t i = 0; i < p.size() && i < components.size(); i++)
    p[i]->addDerivative(components[i],variables,d1,d2);
}


template <class ParentFunc, class Ret, class Scalar>
EvalMultiFunction<ParentFunc,Ret,Scalar>::
EvalMultiFunction (const std::string& functions,
                   const std::string& variables,
                   const Real epsX, const Real epsT)
  : EvalFunctions<Scalar>(functions,variables,epsX,epsT)
{
  std::tie(this->nsd, this->ncmp) = getNoDims<Ret>(this->p.size());
}


template <class ParentFunc, class Ret, class Scalar>
Ret EvalMultiFunction<ParentFunc,Ret,Scalar>::
evaluate (const Vec3& X) const
{
  std::vector<Real> tmp;
  tmp.reserve(this->p.size());
  for (const std::unique_ptr<FuncType>& f : this->p)
    tmp.push_back((*f)(X));

  return Ret(tmp);
}


template<class ParentFunc, class Ret, class Scalar>
Ret EvalMultiFunction<ParentFunc,Ret,Scalar>::
deriv (const Vec3& X, int dir) const
{
  std::vector<Real> tmp;
  tmp.reserve(this->p.size());
  for (const std::unique_ptr<FuncType>& f : this->p)
    tmp.push_back(f->deriv(X,dir));

  return Ret(tmp);
}


template<class ParentFunc, class Ret, class Scalar>
Ret EvalMultiFunction<ParentFunc,Ret,Scalar>::
dderiv (const Vec3& X, int i, int j) const
{
  std::vector<Real> tmp;
  tmp.reserve(this->p.size());
  for (const std::unique_ptr<FuncType>& f : this->p)
    tmp.push_back(f->dderiv(X,i,j));

  return Ret(tmp);
}


template <class ParentFunc, class Ret, class Scalar>
std::vector<Real>
EvalMultiFunction<ParentFunc,Ret,Scalar>::
evalGradient (const Vec3& X) const
{
  std::vector<Vec3> dx;
  dx.reserve(this->p.size());
  for (const std::unique_ptr<FuncType>& f : this->p)
    dx.push_back(f->gradient(X));

  std::vector<Real> result;
  result.reserve(this->ncmp*this->nsd);
  for (size_t d = 1; d <= this->nsd; ++d)
    for (size_t i = 0; i < this->ncmp; ++i)
      result.push_back(dx[i](d));

  return result;
}


template <class ParentFunc, class Ret, class Scalar>
std::vector<Real>
EvalMultiFunction<ParentFunc,Ret,Scalar>::
evalHessian (const Vec3& X) const
{
  std::vector<SymmTensor> dx;
  dx.reserve(this->p.size());
  for (const std::unique_ptr<FuncType>& f : this->p)
    dx.push_back(f->hessian(X));

  std::vector<Real> result;
  result.reserve(this->p.size()*this->nsd*this->nsd);
  for (size_t d2 = 1; d2 <= this->nsd; ++d2)
    for (size_t d1 = 1; d1 <= this->nsd; ++d1)
      for (size_t i = 0; i < this->p.size(); ++i)
        result.push_back(dx[i](d1,d2));

  return result;
}


template <class ParentFunc, class Ret, class Scalar>
std::vector<Real>
EvalMultiFunction<ParentFunc,Ret,Scalar>::
evalTimeDerivative (const Vec3& X) const
{
  std::vector<Real> result;
  result.reserve(this->ncmp);
  for (const std::unique_ptr<FuncType>& f : this->p)
    result.push_back(f->timeDerivative(X));

  return result;
}


RealFunc* utl::parseExprRealFunc (const std::string& function, bool autodiff)
{
  ExprEval::numError = 0;

  RealFunc* f;
  if (autodiff)
    f = new EvalFuncSpatial<autodiff::var>(function.c_str());
  else
    f = new EvalFunction(function.c_str());

  if (ExprEval::numError == 0)
    return f;

  delete f;
  return nullptr;
}


VecFunc* utl::parseExprVecFunc (const std::string& function, bool autodiff)
{
  ExprEval::numError = 0;

  VecFunc* f;
  if (autodiff)
    f = new EvalMultiFunction<VecFunc,Vec3,autodiff::var>(function, "");
  else
    f = new EvalMultiFunction<VecFunc,Vec3,Real>(function, "");

  if (ExprEval::numError == 0)
    return f;

  delete f;
  return nullptr;
}


template class ExpressionHolder<Real>;
template class ExpressionHolder<autodiff::var>;
template class EvalFuncScalar<Real>;
template class EvalFuncScalar<autodiff::var>;
template class EvalFuncSpatial<Real>;
template class EvalFuncSpatial<autodiff::var>;
template class EvalFunctions<Real>;
template class EvalFunctions<autodiff::var>;
template class EvalMultiFunction<VecFunc,Vec3,Real>;
template class EvalMultiFunction<VecFunc,Vec3,autodiff::var>;
template class EvalMultiFunction<TensorFunc,Tensor,Real>;
template class EvalMultiFunction<TensorFunc,Tensor,autodiff::var>;
template class EvalMultiFunction<STensorFunc,SymmTensor,Real>;
template class EvalMultiFunction<STensorFunc,SymmTensor,autodiff::var>;
