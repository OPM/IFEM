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
#include "Vec3.h"
#include "Tensor.h"
#include "expreval.h"
#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <autodiff/reverse/var.hpp>


template<> int EvalFuncScalar<Real>::numError = 0; //!< Explicit instantiation


namespace {

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
  EvalFuncScalar<Real>::numError++;
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
std::pair<size_t,size_t> getNoDims(size_t psize);


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
EvalFuncScalar<Scalar>::EvalFuncScalar (const char* function, const char* x, Real eps)
  : dx(eps)
{
  try {
#ifdef USE_OPENMP
    size_t nalloc = omp_get_max_threads();
#else
    size_t nalloc = 1;
#endif
    expr.resize(nalloc);
    f.resize(nalloc);
    v.resize(nalloc);
    expr.resize(nalloc);
    arg.resize(nalloc);
    for (size_t i = 0; i < nalloc; ++i) {
      expr[i] = std::make_unique<Expression>();
      f[i] = std::make_unique<FunctionList>();
      v[i] = std::make_unique<ValueList>();
      f[i]->AddDefaultFunctions();
      v[i]->AddDefaultValues();
      v[i]->Add(x,0.0,false);
      expr[i]->SetFunctionList(f[i].get());
      expr[i]->SetValueList(v[i].get());
      expr[i]->Parse(function);
      arg[i] = v[i]->GetAddress(x);
    }
  }
  catch (ExprEval::Exception& e) {
    ExprException(e,"parsing",function);
  }
}


template<class Scalar>
EvalFuncScalar<Scalar>::~EvalFuncScalar () = default;


template<class Scalar>
void EvalFuncScalar<Scalar>::addDerivative (const std::string& function, const char* x)
{
  if (!gradient)
    gradient = std::make_unique<FuncType>(function.c_str(),x);
}


template<class Scalar>
Real EvalFuncScalar<Scalar>::evaluate (const Real& x) const
{
  Real result = Real(0);
  size_t i = 0;
#ifdef USE_OPENMP
  i = omp_get_thread_num();
#endif
  if (i >= arg.size())
    return result;
  try {
    *arg[i] = x;
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


template<>
Real EvalFunc::deriv (Real x) const
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

  size_t i = 0;
#ifdef USE_OPENMP
  i = omp_get_thread_num();
#endif
  try {
    *arg[i] = x;
    return derivativesx(expr[i]->Evaluate(),
                        autodiff::wrt(*this->arg[i]))[0].expr->val;
  }
  catch (ExprEval::Exception& e) {
    ExprException(e,"evaluating expression");
  }

  return 0.0;
}


template<class Scalar>
EvalFuncSpatial<Scalar>::
EvalFuncSpatial (const char* function, Real epsX, Real epsT)
  : dx(epsX), dt(epsT)
{
  try {
#ifdef USE_OPENMP
    size_t nalloc = omp_get_max_threads();
#else
    size_t nalloc = 1;
#endif
    expr.resize(nalloc);
    f.resize(nalloc);
    v.resize(nalloc);
    expr.resize(nalloc);
    arg.resize(nalloc);
    for (size_t i = 0; i < nalloc; ++i) {
      expr[i] = std::make_unique<Expression>();
      f[i] = std::make_unique<FunctionList>();
      v[i] = std::make_unique<ValueList>();
      f[i]->AddDefaultFunctions();
      v[i]->AddDefaultValues();
      v[i]->Add("x",0.0,false);
      v[i]->Add("y",0.0,false);
      v[i]->Add("z",0.0,false);
      v[i]->Add("t",0.0,false);
      expr[i]->SetFunctionList(f[i].get());
      expr[i]->SetValueList(v[i].get());
      expr[i]->Parse(function);
      arg[i].x = v[i]->GetAddress("x");
      arg[i].y = v[i]->GetAddress("y");
      arg[i].z = v[i]->GetAddress("z");
      arg[i].t = v[i]->GetAddress("t");
    }
  }
  catch (ExprEval::Exception& e) {
    ExprException(e,"parsing",function);
  }

  // Checking if the expression is time-independent
  // by searching for the occurance of 't' where the next character is not a letter
  IAmConstant = true;
  std::string expr(function);
  for (size_t i = expr.find_first_of('t'); IAmConstant; i = expr.find_first_of('t',i+1))
    if (i >= expr.size())
      return;
    else if (i+1 == expr.size() || !isalpha(expr[i+1]))
      IAmConstant = false;
}


template<class Scalar>
EvalFuncSpatial<Scalar>::~EvalFuncSpatial () = default;


template<class Scalar>
void EvalFuncSpatial<Scalar>::
addDerivative (const std::string& function,
               const std::string& variables,
               int d1, int d2)
{
  if (d1 > 0 && d1 <= 4 && d2 < 1) // A first order derivative is specified
  {
    if (!derivative1[--d1])
      derivative1[d1] = std::make_unique<FuncType>((variables+function).c_str());
  }
  else if ((d1 = voigtIdx(d1,d2)) >= 0) // A second order derivative is specified
  {
    if (!derivative2[d1])
      derivative2[d1] = std::make_unique<FuncType>((variables+function).c_str());
  }
}


template<class Scalar>
Real EvalFuncSpatial<Scalar>::evaluate (const Vec3& X) const
{
  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  Real result = Real(0);
  try {
    size_t i = 0;
#ifdef USE_OPENMP
    i = omp_get_thread_num();
#endif
    if (i >= arg.size())
      return result;

    *arg[i].x = X.x;
    *arg[i].y = X.y;
    *arg[i].z = X.z;
    *arg[i].t = Xt ? Xt->t : Real(0);
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


template<>
Real EvalFuncSpatial<Real>::deriv (const Vec3& X, int dir) const
{
  if (dir < 1)
    return Real(0);
  else if (dir < 4)
  {
    if (derivative1[--dir])
      return derivative1[dir]->evaluate(X);

    // Evaluate spatial derivative using central difference
    Vec4 X0, X1;
    X0.assign(X); X0[dir] -= 0.5*dx;
    X1.assign(X); X1[dir] += 0.5*dx;
    return (this->evaluate(X1) - this->evaluate(X0)) / dx;
  }
  else if (!IAmConstant)
  {
    if (derivative1[3])
      return derivative1[3]->evaluate(X);

    // Evaluate time-derivative using central difference
    Vec4 X0, X1;
    X0.assign(X); X0.t -= 0.5*dt;
    X1.assign(X); X1.t += 0.5*dt;
    return (this->evaluate(X1) - this->evaluate(X0)) / dt;
  }
  else
    return Real(0);
}


template<>
Real EvalFuncSpatial<autodiff::var>::deriv (const Vec3& X, int dir) const
{
  if (dir < 1 || dir > 4)
    return Real(0);

  size_t i = 0;
#ifdef USE_OPENMP
  i = omp_get_thread_num();
#endif

  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  *arg[i].x = X.x;
  *arg[i].y = X.y;
  *arg[i].z = X.z;
  *arg[i].t = Xt ? Xt->t : Real(0);

  // Evaluate spatial derivative using auto-diff
  return derivativesx(expr[i]->Evaluate(),
                      autodiff::wrt(arg[i].get(dir)))[0].expr->val;
}


template<>
Real EvalFuncSpatial<Real>::dderiv (const Vec3& X, int d1, int d2) const
{
  if ((d1 = voigtIdx(d1,d2)) < 0)
    return Real(0);

  return derivative2[d1] ? derivative2[d1]->evaluate(X) : Real(0);
}


template<>
Real EvalFuncSpatial<autodiff::var>::dderiv (const Vec3& X, int d1, int d2) const
{
  if (d1 < 1 || d1 > 3 ||
      d2 < 1 || d2 > 3)
    return Real(0);

  size_t i = 0;
#ifdef USE_OPENMP
  i = omp_get_thread_num();
#endif

  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  *arg[i].x = X.x;
  *arg[i].y = X.y;
  *arg[i].z = X.z;
  *arg[i].t = Xt ? Xt->t : Real(0);

  return derivativesx(derivativesx(expr[i]->Evaluate(),
                                   autodiff::wrt(arg[i].get(d1)))[0],
                                   autodiff::wrt(arg[i].get(d2)))[0].expr->val;
}


template<>
Vec3 EvalFuncSpatial<autodiff::var>::gradient (const Vec3& X) const
{
  size_t i = 0;
#ifdef USE_OPENMP
  i = omp_get_thread_num();
#endif

  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  double t;
  *arg[i].x = X.x;
  *arg[i].y = X.y;
  *arg[i].z = X.z;
  *arg[i].t = t = Xt ? Xt->t : Real(0);

  const auto dx = derivativesx(expr[i]->Evaluate(),
                               autodiff::wrt(*arg[i].x, *arg[i].y, *arg[i].z));

  return Vec4(dx[0].expr->val, dx[1].expr->val, dx[2].expr->val, t);
}


template<>
SymmTensor EvalFuncSpatial<autodiff::var>::hessian (const Vec3& X) const
{
  size_t i = 0;
#ifdef USE_OPENMP
  i = omp_get_thread_num();
#endif

  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  *arg[i].x = X.x;
  *arg[i].y = X.y;
  *arg[i].z = X.z;
  *arg[i].t = Xt ? Xt->t : Real(0);

  const auto dx =
    derivativesx(expr[i]->Evaluate(), autodiff::wrt(*arg[i].x, *arg[i].y, *arg[i].z));

  const auto [uxx, uxy, uxz] =
    derivativesx(dx[0], autodiff::wrt(*arg[i].x, *arg[i].y, *arg[i].z));

  const auto [uyy, uyz] =
    derivativesx(dx[1], autodiff::wrt(*arg[i].y, *arg[i].z));

  const auto [uzz] =
    derivativesx(dx[2], autodiff::wrt(*arg[i].z));

  return SymmTensor({uxx.expr->val, uyy.expr->val, uzz.expr->val,
                     uxy.expr->val, uyz.expr->val, uxz.expr->val});
}


template<class Scalar>
void EvalFuncSpatial<Scalar>::setParam (const std::string& name, double value)
{
  for (std::unique_ptr<ValueList>& v1 : v) {
    Scalar* address = v1->GetAddress(name);
    if (!address)
      v1->Add(name,value,false);
    else
      *address = value;
  }
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
                                           const std::string& variables, int d1, int d2)
{
  std::vector<std::string> components = splitComps(functions,variables);
  for (size_t i = 0; i < p.size() && i < components.size(); i++)
    p[i]->addDerivative(components[i],variables,d1,d2);
}


template <class ParentFunc, class Ret, class Scalar>
Ret EvalMultiFunction<ParentFunc,Ret,Scalar>::
evaluate (const Vec3& X) const
{
  std::vector<Real> res_array(this->p.size());
  for (size_t i = 0; i < this->p.size(); ++i)
    res_array[i] = (*this->p[i])(X);

  return Ret(res_array);
}


template <class ParentFunc, class Ret, class Scalar>
void EvalMultiFunction<ParentFunc,Ret,Scalar>::setNoDims ()
{
  std::tie(this->nsd, this->ncmp) = getNoDims<Ret>(this->p.size());
}


template<class ParentFunc, class Ret, class Scalar>
Ret EvalMultiFunction<ParentFunc,Ret,Scalar>::
deriv (const Vec3& X, int dir) const
{
  std::vector<Real> tmp(this->p.size());
  for (size_t i = 0; i < this->p.size(); ++i)
    tmp[i] = this->p[i]->deriv(X,dir);

  return Ret(tmp);
}


template<class ParentFunc, class Ret, class Scalar>
Ret EvalMultiFunction<ParentFunc,Ret,Scalar>::
dderiv (const Vec3& X, int d1, int d2) const
{
  std::vector<Real> tmp(this->p.size());
  for (size_t i = 0; i < this->p.size(); ++i)
    tmp[i] = this->p[i]->dderiv(X,d1,d2);

  return Ret(tmp);
}


template <class ParentFunc, class Ret, class Scalar>
std::vector<Real>
EvalMultiFunction<ParentFunc,Ret,Scalar>::
evalGradient (const Vec3& X) const
{
  std::vector<Real> result;
  result.reserve(this->ncmp*this->nsd);
  std::vector<Vec3> dx;
  dx.reserve(this->p.size());
  for (const std::unique_ptr<FuncType>& f : this->p)
    dx.push_back(f->gradient(X));

  for (size_t d = 1; d <= this->nsd; ++d)
    for (size_t i = 1; i <= this->ncmp; ++i)
      result.push_back(dx[i-1][d-1]);

  return result;
}


template <class ParentFunc, class Ret, class Scalar>
std::vector<Real>
EvalMultiFunction<ParentFunc,Ret,Scalar>::
evalHessian (const Vec3& X) const
{
  std::vector<Real> result;
  result.reserve(this->p.size()*this->nsd*this->nsd);
  std::vector<SymmTensor> dx;
  dx.reserve(this->p.size());
  for (const std::unique_ptr<FuncType>& f : this->p)
    dx.push_back(f->hessian(X));

  for (size_t d2 = 1; d2 <= this->nsd; ++d2)
    for (size_t d1 = 1; d1 <= this->nsd; ++d1)
      for (size_t i = 0; i < this->p.size(); ++i)
        result.push_back(dx[i](d1,d2));

  return result;
}


template <class ParentFunc, class Ret, class Scalar>
std::vector<Real>
EvalMultiFunction<ParentFunc,Ret,Scalar>::evalTimeDerivative (const Vec3& X) const
{
  std::vector<Real> result;
  result.reserve(this->ncmp);
  for (const std::unique_ptr<FuncType>& f : this->p)
    result.push_back(f->timeDerivative(X));

  return result;
}


RealFunc* utl::parseExprRealFunc (const std::string& function, bool autodiff)
{
  if (autodiff)
    return new EvalFuncSpatial<autodiff::var>(function.c_str());
  else
    return new EvalFunction(function.c_str());
}


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
